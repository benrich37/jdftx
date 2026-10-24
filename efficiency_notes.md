# Efficiency testing notes

Branch `efficiency_testing`, based on upstream e305ab5d. Two independent commits:

1. `bd608a26` Reuse ultrasoft pseudopotential projections in ElecMinimizer
2. `f047b9b0` Use symmetric (Lowdin) orthonormalization when fillings vary per band

All testing so far is CPU-only (arm64 Mac, Accelerate BLAS, Release build with `EnableProfiling`).
The changes have not yet been run on GPU.

## Motivation

The reference production run was `Patty/PBED3BJ_20Ha_addU/surfs/310_fixed_1x3_trim_Ni_B/p07V/ion_opt/jdftx_run/out`:
- MnO2 slab with 71 Mn atoms; nBands 1528; average nbasis 287k; about 2500 GBRV ultrasoft projectors.
- Settings: `target-mu`, LinearPCM CANDLE, `add-U`, `z-spin`, k-points 2x1x1.
- Hardware: 4 A100s, with build 52548a96.

What that run shows:
- About 16 s per ElecMinimize iteration and 20-80 iterations per ionic step, for 3.7 h in total.
- The LinearPCM solve is cheap: time per energy evaluation is the same with 0 or 27 fluid iterations. The wavefunction work dominates.
- About 37% of the wall time came after each electronic minimization was already within 1e-6 Ha of its final energy.
- Most ionic steps ended with "Probably at roundoff error limit" rather than the 1e-8 criterion.
- In that tail, "Undoing step" does not restore the energy. G comes back about 2e-7 Ha higher, and "Step increased G" grows as alpha shrinks. So the noise floor is not double-precision roundoff.

## Change 1: reuse of pseudopotential projections

For ultrasoft pseudopotentials, O(C) = C detR + V Q (V^C). The minimizer recomputed the projections V^C, and formed O(C) as a full plane-wave array, in places where the projections were already known or only the subspace overlap was needed. Each such product costs nProj x nBands x nbasis, which is larger than the nBands^2 x nbasis products because nProj (~2500) > nBands (1528).

| Location | Before | After |
|---|---|---|
| `ElecVars::orthonormalize` | `C^O(C)` forms O(C) | `dagO()`: detR C^C + VdagC^dag Q VdagC |
| gradient, `HC -= O(C) Hsub` | recomputes V^C, forms O(C) | Hsub from projections; `-Q VdagC Hsub` folded into `projectGrad` (one V*M product instead of three) |
| `ElecMinimizer::constrain` | `C^O(dir)` forms O(dir) | `dagO()` with V^dir |
| `ElecMinimizer::step` | `orthonormalize` recomputes V^C | V^C updated linearly from V^dir, cached in `ElecGradient::VdagC` by `constrain` and kept consistent through `axpy` / `*=` |

The net effect is about 9 projector-sized products per state per iteration reduced to about 4. The changes are algebraically exact: trajectories match baseline to about 1e-10 Ha over 25 iterations. When `dump End HC` is requested, the original path is used so that the dumped HC is unchanged.

## Change 2: Lowdin instead of Cholesky orthonormalization

`orthoMatrix()` is Cholesky-based, so it mixes wavefunction columns asymmetrically. With Haux (smearing) fillings, each column carries its own occupation, so Cholesky orthonormalization moves charge between occupied and partially empty states.

The line minimization takes a trial step and then steps back to the predicted alpha. With a symmetric orthonormalizer, that path is exactly reversible. I checked this algebraically for a direction with C^O d = 0, including the Haux rotations. With Cholesky it is not reversible.

On the CPU analog (below), converged to 1e-11:
- At iteration 38, every step raised G by 7.76e-3 Ha, even at alpha = 8e-8.
- "Undoing step" left G 7.76e-3 Ha higher, with a changed density: EH and Eloc moved by about 0.6 Ha at fixed fillings.
- With `invsqrt()`, the largest "Step increased G" event is 8e-10 Ha, which is genuine roundoff for |G| = 526 Ha.

Lowdin is now used whenever fillings are not constant and uniform. That includes the `orthonormalize` after ionic steps (`IonicMinimizer::step`) and lattice steps.

The extra cost is one nBands^3 eigendecomposition per step. This is the same size as the Haux diagonalization that each step already does.

## Testing

- **Regression suite:** `ctest` passes 11/11 with both changes.
- **Test system:** a 12-atom MnO2(110) slab analog of the production run.
  - Settings: GBRV 20/100 Ha, PBE + D3BJ, z-spin, AFM initial moments, `add-U Mn d 0.06`, LinearPCM CANDLE with H2O / F- / K+, k-points 2x4x1, Fermi smearing 0.001.
  - Starting state: a converged neutral run.
- **A/B method:** the baseline is unmodified e305ab5d built the same way. Baseline and patched runs were run side by side with 9 threads each, so timings are comparable.

| Case | Baseline | Patched | Notes |
|---|---|---|---|
| MnO2, `target-mu -0.25` | 45 iters, 5.2 s/iter, 226 s | 45 iters, 3.5 s/iter, 163 s | final G agrees to 5e-9 Ha |
| MnO2, `target-mu -0.19`, 4 ionic steps | 446 iters, 5.1 s/iter | 441 iters, 4.1 s/iter (change 1 only) | see note below |
| Pt(111) 3 layers, `target-mu -0.17`, CANDLE, from LCAO | 66 iters, 29.2 s | 65 iters, 26.1 s | small cell, few projectors |
| MnO2, 25 fixed iterations (change 1 v1 only) | 110 s | 101 s | before the step/gradient parts were added |
| MnO2 converged to 1e-11 | 1e-6 Ha reached at iter 110 | 1e-6 Ha reached at iter 94 | change 2 avoids the 7.8e-3 Ha jump |

### Effect of change 2 on iteration counts

- It removes the hysteresis jumps, and lowers the energy excess at the start of each ionic step by about 15%.
- It did not reliably reduce iteration counts in these tests: 441 vs 446 in total over 4 ionic steps.
- The MnO2 analog has a slow electronic mode, a plateau at -525.97086 Ha vs a minimum at -525.97233 Ha.
  - Whether a run escapes the plateau or stops on it under the |dG| < 1e-8 criterion depends on details of the path.
  - So the ionic-step comparison between builds is noisy.

### Estimate for the GPU runs (not yet measured)

Projector products dominate the per-iteration cost there as well, so I'd expect roughly 15-25% less time per iteration from change 1. Whether change 2 removes the ~1e-7 Ha noise tail in the GPU runs still needs testing; that tail may also come from the GPU-MPI race (see below).

## Other findings (no code changes)

- **GPU-MPI race:** build 52548a96 predates the upstream GPU-MPI sync fix (544d0e78). np = 4 GPU runs should use a newer build with `JDFTX_GPU_SYNC_MPI=1`.
- **cuSolver:** `EnableCuSolver` defaults to OFF. Without it, the per-step nBands^2 Haux diagonalization (and the new `invsqrt`) run in CPU LAPACK in GPU builds. Check the build's `CMakeCache.txt`.
- **Initial vacuum solve:** on fresh wavefunctions with a fluid, the vacuum pre-solve ran 100 iterations without converging. That took about 1650 s, about 12% of the first job.
- **Convergence threshold:** `energyDiffThreshold 1e-8` on a ~1e4 Ha system sits at or below the observed noise floor.

## Possible next steps

- A/B on Perlmutter: the same ionic relaxation with e305ab5d vs this branch, both with cuSolver and `JDFTX_GPU_SYNC_MPI=1`.
- Cap the threshold or iteration count of the initial vacuum pre-solve; it only provides a starting guess.
- Smaller items:
  - Use `zherk` for the Hermitian C^C overlap.
  - Cache the +U atomic orbitals per state (about 3% here).
  - Apply `dagO` / `augmentOverlapGrad` in BandMinimizer and BandDavidson too.

## Reproducing

Local files live in `~/Desktop/jdftx_hold/`:
- `perf_test/`: inputs (`common.in`, `bench.in`, `tight.in`, `caseB/C/D.in`), outputs, and `run_ab.sh`, which runs the baseline and patched builds side by side.
- `perf_patches/`: the two patches, plus `all_perf_changes.patch` with both combined.
- Builds:
  - `build_base_prof`: baseline e305ab5d.
  - `build_perf`: this branch.
  - All builds are Release builds with `EnableProfiling`.

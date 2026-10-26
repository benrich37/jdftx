# Efficiency testing notes

Branch `efficiency_testing`, based on upstream e305ab5d. Two independent commits:

1. `bd608a26` Reuse ultrasoft pseudopotential projections in ElecMinimizer
2. `f047b9b0` Use symmetric (Lowdin) orthonormalization when fillings vary per band

CPU testing was on an arm64 Mac (Accelerate BLAS, Release build with `EnableProfiling`).
The GPU results are in "GPU results" below (Perlmutter, 4 A100-80GB, 2026-10-24/25).
Summary: change 1 is 27% faster per iteration on GPU, with identical results, and uses 6.5 GiB less memory per GPU.
Change 2 neither helps nor removes the noise floor on GPU.

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

### Estimate for the GPU runs (made before they were measured)

Projector products dominate the per-iteration cost there as well, so I'd expect roughly 15-25% less time per iteration from change 1. Whether change 2 removes the ~1e-7 Ha noise tail in the GPU runs still needs testing; that tail may also come from the GPU-MPI race (see below).
The measurements are in the next section: 27% less per iteration, and the noise tail remains with both the sync fix and change 2.

## GPU results (Perlmutter, 2026-10-24/25)

### Setup

- **Builds:** Perlmutter recipe from CLAUDE.md: CUDA 13.2, `EnableCuSolver=yes`, `CudaAwareMPI=yes`, `PinnedHostMemory=yes`,
  `EnableProfiling=yes`.
  - base: `54315e95` (fork master; contains sync fix 544d0e78).
  - p1: `bd608a26` (change 1).
  - p12: `fe9df38d` (changes 1 + 2).
- **Environment:** every run had `JDFTX_GPU_SYNC_MPI=1` (checked via `gpuInit: sync before MPI is on`) and
  `JDFTX_CACHE_SIZE=36000`. The binaries were called by full path, and `ldd` confirmed each build's own `libjdftx_gpu.so`.
- **System:** the production run `310_fixed_1x3_trim_Ni_B/p07V/ion_opt/jdftx_run` (218 atoms, nBands 1528, 4 states, one per GPU).
  - Restarted from its final `wfns`/`fillings`/`eigenvals`/`fluidState` with the `in` positions, which are the start of that
    relaxation and up to 0.16 bohr from where the wfns were converged. The first G is 2.6 Ha above the minimum.
  - The State, BandProjections and ElecDensity dumps are removed.
- **Layout:** 1 node, `-c 32 --ntasks-per-node=4 --gpus-per-task=1 --gpu-bind=none`, 80GB nodes.
  Each build ran on a different node.

### Regression suite on GPU

`JDFTX_LAUNCH="" JDFTX_SUFFIX=_gpu make test`: p12 11/11 and base 11/11 (47 checks each). p12 and base energies
differ by at most 1.1e-7 Ha (ionSolvation CANDLE).

### Test A: cost per iteration (20 iterations, `energyDiffThreshold 0`)

| | base | p1 | p12 |
|---|---|---|---|
| s/iter | 16.42 | 11.87 (−27.7%) | 11.96 (−27.2%) |
| max \|G − G_base\| | — | 2.5e-10 Ha | 0.14 Ha at iter 1 → 7e-5 at iter 20 (different path) |
| `Y1^Y2` calls / s | 569 / 85.0 | 325 / 47.2 | 325 / 47.3 |
| `Y*M` calls / s | 625 / 112.8 | 297 / 62.4 | 297 / 62.9 |
| `augmentOverlap` (ColumnBundle) calls / s | 336 / 106.7 | 8 / 3.3 | 8 / 2.5 |
| `matrix::diagonalize` calls / s | 62 / 2.5 | 62 / 2.6 | 102 / 4.1 |
| `orthoMatrix(matrix)` calls / s | 41 / 0.74 | 41 / 0.59 | 1 / 0.07 |
| `Idag_DiagV_I` / `diagouterI` / `EdensityAndVscloc` s | 30.0 / 30.0 / 20.1 | 29.9 / 30.0 / 20.3 | 30.0 / 30.0 / 20.2 |
| MEMUSAGE total (ColumnBundle) | 42.43 (40.72) GiB | 35.81 (34.18) GiB | 35.81 (34.18) GiB |

- The new matrix-level `augmentOverlap` and `augmentOverlapGrad` have no StopWatch. Their small (nProj × nBands) products
  appear under `Y*M`/matrix multiplies.
- p1 saves 6.5 GiB of peak memory, exactly one wavefunction-sized array (nbasis 287k × 1528 × 16 B).
  That is the full O(C) array that is no longer formed.
- Change 2's extra `invsqrt` costs about 0.08 s/iter on GPU (cuSolver), under 1%.
- Per-iteration profile for p1/p12 (about 12 s):
  - `Y*M` + `Y1^Y2` ≈ 5.5 s
  - density and local potential (`diagouterI`, `Idag_DiagV_I`) ≈ 3 s
  - augmentation ≈ 1.75 s
  - XC ≈ 1.35 s (`ExCorrCommunication` 0.67 s; this probably includes the sync wait for the preceding kernels)
  - `VanDerWaalsD3::energyAndGrad` costs 18 s per call, but runs only once per ionic evaluation.

### Test B: convergence tail (default thresholds, `nIterations 150`)

| | base | p1 | p12 |
|---|---|---|---|
| iterations, stop | 52, roundoff limit | 52, roundoff limit | 60, roundoff limit |
| s/iter mean (median) | 18.39 (16.39) | 13.29 (11.81) | 13.31 (11.91) |
| "Step increased G": n, max | 9, 2.1e-7 | 9, 2.0e-7 | 12, 4.4e-7 |
| iterations after \|G − Gfinal\| < 1e-6 | 9 (17%) | 9 (17%) | 12 (20%) |
| G final | −9575.920310639 | −9575.920310672 | −9575.920310883 |
| job wall time | 18.2 min | 13.7 min | 15.8 min |

- With the sync on, base still has the ~1e-7 Ha floor, so the floor is not the MPI race. Change 2 does not remove it.
- The signature in p12, iterations 49 → 50:
  - After one trial evaluation, each later evaluation sits about 1.2e-7 Ha higher, even at alpha = 8e-6.
  - "Undoing step" restores C, but G stays 1.17e-7 higher. Some history-dependent state changed.
  - Energy evaluation itself is reproducible: base and p1 agree to ≤ 2.5e-10 Ha for 25 iterations.
- The prime suspect, not yet proven, is the warm-started LinearPCM solve.
  - It is the only iterative solve inside an energy evaluation that starts from the previous state
    (`fluid-minimize knormThreshold 1e-11` on sqrt(|r.z|)). Near convergence most solves "Complete after 0 iterations".
  - Against it: `LinearPCM::get_Adiel_and_grad_internal` uses a first-order-correct energy, so a fluid residual should change
    G only at second order.
- The 7.8e-3 Ha Cholesky hysteresis of the CPU MnO2 analog does not occur in this system, even in base.

### Test C: ionic relaxation (`ionic-minimize nIterations 4` from the `in` geometry)

| | base | p1 | p12 |
|---|---|---|---|
| electronic minimizations | 6 (step 3 needed a cubic line-search step) | 5 | 5 |
| ElecMinimize iterations per minimization | 52, 43, 37, 47, 26, 11 | 52, 33, 42, 46, 18 | 60, 34, 52, 35, 50 |
| total iterations | 216 | 191 | 231 |
| stop | 5× roundoff, 1× converged | 5× roundoff | 5× roundoff |
| "Step increased G" events (max) | 45 (3.9e-7) | 38 (5.3e-7) | 82 (1.2e-6) |
| iterations after \|G − Gfinal\| < 1e-6 | 34 (16%) | 38 (20%) | 73 (32%) |
| wall time to ionic step 4 | 76.5 min | 50.6 min (−34%) | 64.6 min (−16%) |
| wall per iteration incl. forces/D3 | 21.2 s | 15.9 s | 16.8 s |
| G after 4 ionic steps, \|grad\|_K | −9575.925637, 4.1e-4 | −9575.925655, 5.1e-4 | −9575.925627, 5.2e-4 |

- The ionic paths diverge after step 2 because tiny differences pass through the noisy stopping point.
  So iteration counts are single noisy samples, as on CPU.
- p12 was worse than p1 in both B (60 vs 52 iterations) and C (231 vs 191, with twice the "Step increased G" events).

### Conclusions

- **Change 1 (projection reuse): keep.**
  - Results are identical to base at GPU roundoff.
  - It is 27-28% faster per ElecMinimize iteration and about 25% faster overall including forces.
  - Peak memory drops by one wavefunction-sized array: −6.5 GiB/GPU here.
- **Change 2 (Löwdin): no benefit on GPU for this system.**
  - It costs nothing per iteration.
  - It does not remove the ~1e-7 Ha floor.
  - Both B and C needed more iterations, with more failed steps.
  - Keep it out of production builds unless a case like the CPU 7.8e-3 Ha hysteresis turns up.
- **The GPU noise floor is not Cholesky and not the MPI race.** About 16-20% of iterations (p1/base) are spent within 1e-6 Ha of
  the final energy, and most minimizations end with "Probably at roundoff error limit".
- No further code changes were made:
  - cuBLAS `zherk` was not faster than `zgemm` at 1528 × 72k on a (busy) login A100.
  - Caching the +U orbitals would save about 3% for about 1.6 GiB/GPU. That is the wrong trade while 40GB nodes are memory-bound.

### Next steps (GPU)

1. Find the source of the noise floor. Rerun Test B with p1 and `fluid-minimize knormThreshold 1e-14` (and 1e-9).
   If the floor scales with the threshold, it is the fluid warm start. Each run is about 15 min, so it fits in the debug queue.
2. Once the floor is understood, use a looser electronic stop for these 1e4-Ha systems: `energyDiffThreshold` of about 1e-7,
   or a `knormThreshold`. That would save the roughly 16-20% of iterations spent in the tail.
3. Try p1 on 40GB nodes for the 255-atom Ni_B adsorbates. They peaked at 45.4 GiB on base and should be about 39 GiB with change 1,
   which is marginal but may avoid the 2-2.5x slower managed-memory build.
4. Production builds: use change 1 (`bd608a26`) on top of master, and leave out change 2.

## Other findings (no code changes)

- **GPU-MPI race:** build 52548a96 predates the upstream GPU-MPI sync fix (544d0e78). np = 4 GPU runs should use a newer build with `JDFTX_GPU_SYNC_MPI=1`.
- **cuSolver:** `EnableCuSolver` defaults to OFF. Without it, the per-step nBands^2 Haux diagonalization (and the new `invsqrt`) run in CPU LAPACK in GPU builds. Check the build's `CMakeCache.txt`.
- **Initial vacuum solve:** on fresh wavefunctions with a fluid, the vacuum pre-solve ran 100 iterations without converging. That took about 1650 s, about 12% of the first job.
- **Convergence threshold:** `energyDiffThreshold 1e-8` on a ~1e4 Ha system sits at or below the observed noise floor.

## Possible next steps

- ~~A/B on Perlmutter~~: done (see "GPU results"; GPU-specific next steps are listed there).
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

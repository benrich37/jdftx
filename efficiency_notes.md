# Efficiency testing notes

Branch `efficiency_testing`, originally based on upstream e305ab5d, merged with upstream through b3c54e1d on 2026-10-26.
Two changes were tested:

1. `bd608a26` Reuse ultrasoft pseudopotential projections in ElecMinimizer. **Kept**; this is the only code change on the branch
   (true until `a37d055f`, item 3, was added on 2026-10-28).
2. `f047b9b0` Use symmetric (Lowdin) orthonormalization when fillings vary per band. **Tested and found to show no
   improvement; removed** by the revert `5c09112b`. See "Change 2" for how to recover it.
3. `a37d055f` Speed up DFT-D3 pair sums in `VanDerWaalsD3.cpp` (2026-10-28). **Kept**; results unchanged (Mac: 2e-12;
   Perlmutter GPU: 5e-15 Ha). Perlmutter timings are in "DFT-D3 speedup `a37d055f`: Perlmutter GPU results" below.

CPU testing was on an arm64 Mac (Accelerate BLAS, Release build with `EnableProfiling`).
The GPU results are in "GPU results" below (Perlmutter, 4 A100-80GB, 2026-10-24/25).
Perlmutter results for the D3 speedup, and the end-to-end gain of the branch tip over the current production build on
40GB nodes, are in "DFT-D3 speedup `a37d055f`: Perlmutter GPU results" (2026-10-29) at the end.
Summary: change 1 is 27% faster per iteration on GPU, with identical results, and uses 6.5 GiB less memory per GPU.
Change 2 neither helped nor removed the noise floor on GPU, so it was removed.

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

The net effect, per state, species and iteration (one trial step and one accepted step), is 9 projector-sized products reduced to 2: orthonormalization 4 → 0, gradient 3 → 1, `constrain` 2 → 1. The GPU profiler counts match this (see Test A). The changes are algebraically exact: trajectories match baseline to about 1e-10 Ha over 25 iterations. When `dump End HC` is requested, the original path is used so that the dumped HC is unchanged.

## Change 2: Lowdin instead of Cholesky orthonormalization (tested; no improvement; removed)

**Status:** tested on CPU and GPU and found to show no performance improvement (see "Effect of change 2 on iteration
counts" and "GPU results"), so it was removed from the branch by the revert `5c09112b`. It can be recovered from any of:
- the original commit `f047b9b0` (still in the branch history);
- the git stash entry "efficiency_testing change 2: Lowdin orthonormalization when fillings vary ..." (`git stash list`);
- `efficiency_patches/not_adopted/02_lowdin_orthonormalize.patch`, which applies on top of change 1.

The rest of this section describes the change as it was tested.

`orthoMatrix()` is Cholesky-based, so it mixes wavefunction columns asymmetrically. With Haux (smearing) fillings, each column carries its own occupation, so Cholesky orthonormalization moves charge between occupied and partially empty states.

The line minimization takes a trial step and then steps back to the predicted alpha. With a symmetric orthonormalizer, that path is exactly reversible. I checked this algebraically for a direction with C^O d = 0, including the Haux rotations. With Cholesky it is not reversible.

On the CPU analog (below), converged to 1e-11:
- At iteration 38, every step raised G by 7.76e-3 Ha, even at alpha = 8e-8.
- "Undoing step" left G 7.76e-3 Ha higher, with a changed density: EH and Eloc moved by about 0.6 Ha at fixed fillings.
- With `invsqrt()`, the largest "Step increased G" event is 8e-10 Ha, which is genuine roundoff for |G| = 526 Ha.

As tested, Lowdin was used whenever fillings are not constant and uniform. That includes the `orthonormalize` after ionic steps (`IonicMinimizer::step`) and lattice steps.

The extra cost is one nBands^3 eigendecomposition per step. This is the same size as the Haux diagonalization that each step already does.

## Testing

- **Regression suite:** `ctest` passed 11/11 with both changes, and again with change 1 alone after the upstream merge.
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
   - Perlmutter, 2026-10-26/29: with the *managed* build, `74dd6931` peaked at 38.4 GB on the 255-atom input (45.4 on
     `54315e95`) and `4fe71c77` at 36.3 GB on the 224-atom production input, both on 40GB nodes. The non-managed build has
     not been tried on 40GB nodes yet (see the Perlmutter D3 section at the end).
4. Production builds: use change 1 (`bd608a26`) on top of master, and leave out change 2.
   - Perlmutter, 2026-10-26: done as `~/jdftx_compile_test/ben/build-gpu-{,managed-}eff-74dd6931` (validated against the
     `54315e95` builds, job 60345288). 2026-10-29: `$JDFTx_GPU` in `~/.jdft.bashrc` still points at
     `build-gpu-managed-54315e95`; builds of the tip `4fe71c77` exist only on scratch (`claude_scratch/eff_d3/`).

## Merge with upstream (2026-10-26)

Upstream e305ab5d..b3c54e1d (the `vib_chk` vibrations checkpointing work) changes only `commands/dump.cpp`,
`electronic/Vibrations.{cpp,h}`, and one comment in `electronic/SpeciesInfo.h`. None of it touches the code that change 1
modifies. The new vibrations code moves atoms only through `IonicMinimizer::step`, which calls `orthonormalize(q)` with
the default `projected=false`, so the cached projections are recomputed there as before.

Checks after the merge (CPU):
- `01_uspp_projection_reuse.patch` (and the not-adopted `02` on top of it) applies cleanly to b3c54e1d, and the result is
  identical, file for file, to the branch code.
- After removing change 2, the branch's code differs from upstream b3c54e1d by exactly `01_uspp_projection_reuse.patch`.
- `ctest`: 11/11, including `vibrations`.
- MnO2 benchmark (`bench.in`, 25 iterations, baseline and branch side by side at 9 threads): energies agree with the
  baseline to 6.7e-11 Ha at every iteration, and time per iteration is 3.59 vs 4.50 s (−20%).

## Other findings (no code changes)

- **GPU-MPI race:** build 52548a96 predates the upstream GPU-MPI sync fix (544d0e78). np = 4 GPU runs should use a newer build with `JDFTX_GPU_SYNC_MPI=1`.
- **cuSolver:** `EnableCuSolver` defaults to OFF. Without it, the per-step nBands^2 Haux diagonalization (and change 2's `invsqrt`, now removed) run in CPU LAPACK in GPU builds. The Perlmutter recipe in CLAUDE.md already enables it. Whether the shared m4025 build used for the analysed run has it is unknown; check its `CMakeCache.txt`.
- **Initial vacuum solve:** on fresh wavefunctions with a fluid, the vacuum pre-solve ran 100 iterations without converging. That took about 1650 s, about 12% of the first job.
- **Convergence threshold:** `energyDiffThreshold 1e-8` on a ~1e4 Ha system sits at or below the observed noise floor.

## Possible next steps

- ~~A/B on Perlmutter~~: done (see "GPU results"; GPU-specific next steps are listed there).
- Cap the threshold or iteration count of the initial vacuum pre-solve; it only provides a starting guess.
- Smaller items:
  - `zherk` for the Hermitian C^C overlap: tested on GPU, not faster than `zgemm`; not pursued.
  - Caching the +U atomic orbitals per state: estimated at about 3% for about 1.6 GiB/GPU; not pursued.
  - Apply `dagO` / `augmentOverlapGrad` in BandMinimizer and BandDavidson too (untested; same algebra as change 1).

## Reproducing

Local files live in `~/Desktop/jdftx_hold/`:
- `perf_test/`: inputs (`common.in`, `bench.in`, `tight.in`, `caseB/C/D.in`), outputs, and `run_ab.sh`, which runs the baseline and patched builds side by side.
- Patch files (in this branch, `efficiency_patches/`): `01_uspp_projection_reuse.patch` (change 1) and
  `not_adopted/02_lowdin_orthonormalize.patch` (change 2, applies on top of change 1).
- Builds:
  - `build_base_prof`: baseline e305ab5d.
  - `build_perf`: changes 1 + 2 as tested before the merge.
  - `build_merged`: this branch after the upstream merge (change 1 only).
  - All builds are Release builds with `EnableProfiling`.

## Survey of further efficiency candidates (2026-10-28)

Question: what else can be made faster while keeping JDFTx functionally equivalent, and would a documented
minimization algorithm beat the current one? Sources: the code on this branch (merged with upstream through b3c54e1d;
upstream has nothing newer), the GPU profile of Test A (local copy `~/Desktop/jdftx_hold/efficiency_gpu_outputs/runs/A_p1/out`), and CPU experiments on
the Mac with a profiling build of a37d055f (`~/Desktop/jdftx_hold/build_eff`, source clone `jdftx_eff`). Inputs, outputs,
scripts, prototype patches and the toy model are in `~/Desktop/jdftx_hold/efficiency_survey/`.

### Where a GPU ElecMinimize iteration goes (A_p1: 218 atoms, 1 state per A100, 11.9 s per iteration)

One iteration = one energy-only (trial) evaluation + one energy-and-gradient evaluation + `constrain(d)`.
StopWatch totals divided by 20 iterations; call counts per iteration. The totals include the iteration-0 evaluation
(21 gradient and 41 energy evaluations) and the force calculation, so the rows sum to ~12.3 s against the 11.9 s
measured between iterations 1 and 20. The 3 D3 calls in this profile cost 18.4 s each; a37d055f removes that.

| Part | s/iter | share | calls | what it is |
|---|---|---|---|---|
| `Y*M` + `Y1^Y2` | 5.5 | 46% | 15 + 16 | ~13 full nbasis x nBands x nBands products, see below |
| `Idag_DiagV_I` + `diagouterI` | 3.0 | 25% | 1 + 2 | 1528 bands x 4 single FFTs (2 + 1 + 1) of the 70x240x300 wfns grid |
| `augmentDensityGridGrad` | 1.25 | 10% | 4 (species) | GPU kernel launched once per atom, see item 1 |
| `augmentDensityGrid` | 0.50 | 4% | 8 | per atom too, but a plain forward kernel |
| `ExCorrCommunication` | 0.69 | 6% | ~12 timed blocks | ~20 grid-sized allReduces (sigma, E_n, E_sigma, E_nTilde); 99% of `ExCorrTotal` |
| rest of `EdensityAndVscloc` (fluid, Hartree, local) | 0.33 | 3% | 2 | up to 44 LinearPCM CG iterations per evaluation |
| `rhoAtom_calc` + `rhoAtom_grad` (+U) | 0.40 | 3% | 8 + 4 (per species; only Mn has U) | orbitals rebuilt and re-projected each call |
| `EnlAndGrad`, `augmentDensitySpherical(Grad)`, dense nBands^3 (diagonalize, inv, orthoMatrix) | 0.7 | 5% | | |

The ~13 full-size products per iteration (0.4-0.5 s each): 4 for frame bookkeeping in `step` (`dir.C*rotPrevC` and
the orthonormalizing `C*rot`, in both the trial and the actual step), 2 for rotating `grad` and `Kgrad` back by
`rotPrevCinv` in `compute`, 1 for `C^HC`, 1 for `HC -= C*Hsub`, ~1.6 for `projectGrad`, ~1.6 for `V^d` in
`constrain`, 1 for `d -= C*(C^Od)`, ~0.7 for the +U orbital projections. The trial evaluation costs about 2.7 s (22%).

### Functionally equivalent candidates, ranked by expected GPU gain

1. **`nAugmentGrad_gpu`: per-launch overhead and one launch per atom (~6% of the iteration, up to 9%).**
   `SpeciesInfo::augmentDensityGridGrad` loops over spins and atoms and calls `nAugmentGrad` for each (436 launches per
   gradient evaluation here). The GPU version (`SpeciesInfo.cu`) costs 0.28 s per species call, 4.6x the forward
   `nAugment` kernel (0.061 s), while on CPU the ratio is 2.0x (0.092 vs 0.045 s per call, `b10.in` at 8 threads; the
   gradient kernel does more work per point: spline derivatives and 6 accumulators). The excess over 2x is overhead,
   with two visible sources. (a) Every call runs `cudaGetDevice`, `cudaGetDeviceProperties` and `cudaFuncGetAttributes`
   on the host to size the launch; `cudaGetDeviceProperties` alone can take milliseconds, and the whole launch costs
   about 2 ms for Mn (0.28 s / 142 launches). `GpuUtil.cpp` already stores the device properties at init, and the
   other kernels' `GpuLaunchConfig` only calls `cudaFuncGetAttributes`. (b) The kernel uses `nCoeff` blocks (the |G|
   bins of `nagIndex`) of `warpSize` threads, about 16k threads for a 7M-point grid, and ends with
   `cudaDeviceSynchronize()` because its scratch buffer is a local `GpuBuffer` (pool-allocated, so the allocation is
   cheap; the sync protects its lifetime). Step 1, a few lines and strictly equivalent: cache the launch configuration
   in statics and keep a persistent scratch buffer so the sync can go. Step 2, if the kernel is still slow: make the
   atom a second grid dimension (`blockIdx.y`), index `atpos`, `nRadial` and the `E_nRadialTemp` copies by atom, then
   one collect kernel and one sync per species call. Same arithmetic per grid point, only the reduction grouping
   changes (roundoff-level, like the existing 8-copy split). Expected 1.25 s -> ~0.5 s per iteration if the GPU reaches
   the CPU's 2x ratio (~6%), ~0.15 s if it matches the forward kernel (~9%). Needs a GPU build to measure; the CPU
   path is untouched. `augmentDensityGrid` could be batched the same way (one kernel per species reading `n` once per
   point); ~2%.
2. **XC: skip the MPI grid split on GPUs (~5%).** `ExCorr::operator()` divides the three gradient directions
   (`TaskDivision(3, mpiWorld)`) and the functional evaluation (`gInfo.irStart/irStop`) over processes and then
   allReduces sigma (3), E_n (2), E_sigma (3) and E_nTilde (2) grid arrays (10 per call, timed as 6 blocks): 0.3 s per
   XC call on 4 GPUs, 99% of the XC time, while the
   functional itself takes 1 ms (`ExCorrFunctional`) and the FFTs the split saves take ~10 ms. Prototype:
   `efficiency_survey/patches/xc_nosplit_gpu.prototype.patch` (`splitMPI = not isGpuEnabled()`, plus an env override
   `JDFTX_XC_NOSPLIT` used only for the CPU check). CPU check, np=2, 10 iterations of `b10.in`: energies agree with the
   split path to 1e-10 at every iteration (for scale: the same input at np=1 with 8 threads differs from the np=2 run
   by 3e-10, and the thread-threshold variants below by 4e-11), `ExCorrCommunication` disappears and
   `ExCorrFunctional` doubles per process as expected (so on CPU the split stays).
3. **Batched FFTs in `Idag_DiagV_I` / `diagouterI` (2-5%, uncertain).** Both transform one band at a time
   (`getColumn`: memset + scatter, cuFFT, multiply, cuFFT, gather): 2 transforms per band in `Idag_DiagV_I` and 1 in
   each of the two `diagouterI` calls, about 6100 single 5M-point transforms per iteration. At ~0.5 ms each they are
   already near A100 bandwidth, so a `cufftPlanMany` batch would mostly remove the launch and
   allocation overhead of the ~8 small kernels per band. I would not expect more than 10-20% of the 3 s. Not tested
   (no CUDA on the Mac).
4. **Reuse the +U projections within an evaluation (~1%).** `rhoAtom_grad` rebuilds the orbitals and recomputes
   `Opsi^Cq` that `rhoAtom_calc` just computed for the same `C[q]`. Caching `psiOCdag` per state (355 x 1528) saves one
   0.23-size product and one orbital build per iteration. Caching the orbitals themselves is the 3% / 1.6 GiB item above.
5. **Frame bookkeeping (up to 4%, not exactly equivalent).** One of the 4 rotation products per iteration could go by
   keeping `d` and `gPrev` in the current frame and rotating them in `step`, but only by dropping the re-unitarization
   snap in `report()`, which perturbs the CG direction at O(alpha^2). Not pursued.
6. **CPU only: a thread launch per small zgemm.** Without MKL, `eblas_zgemm` (`BlasExtra.cpp`) spawns
   `nProcsAvailable-1` threads for every product, including the 18x84x18 ones in `augmentDensitySpherical` and all
   `matrix*matrix` products. Prototype threshold (`JDFTX_ZGEMM_THREAD_MIN`, `patches/zgemm_thread_threshold.prototype.patch`),
   `b10.in` (10 iterations, 12-atom slab, 8 threads, runs one at a time, `efficiency_survey/zgemm_threshold/`):
   wall 40.1 s (HEAD) and 45.8 s (prototype with threshold 1, i.e. always threaded; the spread shows the timing noise)
   vs 35.5 / 34.8 / 34.9 s with thresholds 1e5 / 1e6 / 1e7; `augmentDensitySpherical` 0.41 -> 0.09 s, its gradient
   0.23 -> 0.06 s, `EnlAndGrad` 0.75 -> 0.16 s, and the nBands-sized `Y*M`/`Y1^Y2` products 12.5 -> 10.9 s; energies
   agree to 4e-11 at every iteration. So 10-25% on a small
   CPU system; the gain shrinks with system size as the big products dominate. A threshold of 1e6 complex MACs
   (about 0.3 ms single-threaded) is safe; the final patch should hard-code it and drop the env variable.
   Irrelevant on GPU (cuBLAS) and with MKL (`THREADED_BLAS`).

### Algorithms (not functionally equivalent): caseB on the 12-atom MnO2 slab

Restart from the neutral state at `target-mu -0.25`, converge to 1e-8; 8 threads, two runs side by side
(`efficiency_survey/algo/`).

| Scheme | iterations | evaluations (fluid solves) | wall | final G |
|---|---|---|---|---|
| CG Polak-Ribiere + quadratic linmin (default) | 45 | 92 (46 with gradient) | 2:33 | -519.497837474793 |
| CG Fletcher-Reeves | 45 | 92 | 2:42 | -519.497837474786 |
| CG + cubic Wolfe linmin (`linminMethod CubicWolfe`) | 71 | 73 (all with gradient) | 3:16 | -519.497837252730 |
| L-BFGS, default Wolfe linmin | 6, stopped at "roundoff" | 12 | 0:37 | -519.4927 |
| L-BFGS, `alphaTstart 0.1 nAlphaAdjustMax 6` (default history 15 / `history 5`) | 87 / 85 | 156 / 177 | 10:33 / 12:41 | diverged, aborted / -519.458 (roundoff stop) |
| L-BFGS, `linminMethod Quad` | 106 | 247 | 14:38 | diverged (-161) |
| SCF, Pulay density mixing + Davidson (`electronic-scf`, defaults) | 72 cycles | 144 | 5:30 | -519.497837438154 |

- **L-BFGS** is unusable as implemented: the first tentative step is `alphaTstart` along the preconditioned gradient
  and overshoots by ~100x ((E-E0)/|g.d| = 47; g.d changes sign and grows 134x), the cubic search gives up after
  `nAlphaAdjustMax` tries without satisfying Wolfe, and `Minimize_lBFGS.h` then stores `rho = 1/(y.s)` without checking
  the curvature condition y.s > 0 (the standard skip rule; Nocedal and Wright, Numerical Optimization, ch. 6-7), so
  the inverse-Hessian update goes indefinite and the next direction is not a descent direction ("Bad step direction:
  g.d > 0"). Smaller initial steps and the quadratic linmin did not help: `alphaTstart 0.1` diverged to -1.3e5 Ha and
  aborted on the hermiticity check in `diagonalize`, its `history 5` variant stopped at -519.458 ("roundoff"), and the
  Quad variant diverged to -161 Ha. Its large steps also show the irreversibility below (Quad variant: G up by 0.15
  and 0.17 Ha after the undos at iterations 2 and 28). Making it work would need the curvature guard, a
  Hessian-scaled first step, and care with the frame: `report()` re-snaps `rotPrevC` to the unitary `rotPrev` each
  iteration, so s = alpha*d and y = g - gPrev stored in different iterations are consistent only to O(alpha^2). The
  prize is bounded by the trial evaluation (22% on GPU) plus any iteration-count difference.
- Not tested: `linminMethod Relax` (constant step `alphaTstart`, no line search, only checks that G is finite; the
  code marks it "not recommended for CG"). It would remove the trial evaluation but gives up the step-size control
  these runs need (alpha ranges from 1e-4 to 1.5 within one run here).
- **CG with the Wolfe search** does reach one evaluation per iteration, but each is a full gradient evaluation (9.3 s on
  GPU vs 2.7 s for the trial energy) and the inexact line search costs 60% more iterations. Net loss.
- **SCF + Davidson** converges but needs 72 cycles with a wild start (cycle 1 at -426 Ha) and 2x the wall time; mixing
  was not tuned. For `target-mu` + fluid + metallic slab the direct minimizer stays the right default. (Davidson's own
  cost could still drop with `dagO`/cached projections, as in change 1.)
- Fletcher-Reeves follows the same path as Polak-Ribiere here (G agrees to 1e-14; both runs had one CG reset, "Bad
  step direction", in 45 iterations).

### Tight convergence: the step is not reversible with Cholesky orthonormalization (CPU, reproducible)

`caseA_tight.in` (target-mu -0.19 restart, `energyDiffThreshold 1e-11`, 200 iterations) on three builds:

| build | result in 200 iterations | undo events | alpha-independent energy jumps | final G | min G |
|---|---|---|---|---|---|
| HEAD a37d055f (Cholesky, projection reuse) | not converged, grad 3e-7 | 2 | 1: +9.9e-4 Ha at iter 38 | -525.972332414 | same |
| base_prof e305ab5d (Cholesky, no projection reuse) | not converged, grad 3.6e-7 | 1 | 1: +7.9e-4 Ha at iter 38 | -525.972337164 | same |
| perf (Lowdin `invsqrt` + projection reuse, = patch 02) | roundoff stop at iter 137 (\|dG\| 2e-10, grad 8e-9) | 8 (17 failed steps, 7e-12 to 5.5e-10) | none | -525.972332134 | same |

What the jump looks like (HEAD log, iteration 38): the test step at alpha 0.61 (the previous step size, via
`updateTestStepSize`) is followed by "Step increased G by 9.9178e-04 / 9.9175e-04 / 9.9175e-04" at alpha 8e-5, 8e-6,
8e-7, then "Undoing step". The fluid solver needs 21 CG iterations at the first stepped-back point and 0 at the next
two and after the undo, so those four states coincide to the fluid tolerance although alpha spans 8e-5 to 0, and all
sit 1 mHa above G at iteration 38. The next `FillingsUpdate` has the same nElectrons (108.264800) and total moment
(+0.23970) but |m| 14.32447 instead of 14.32590. So fillings and Haux come back exactly and the wavefunctions do not.
Both Cholesky runs then leave the plateau at -525.97086 (where they had crept at 1e-7 to 3e-7 per iteration since
iteration ~25) and descend to -525.97233 by iteration ~200. The Lowdin run meets the same hard spot at iterations
38-41 (trial steps overshoot; the fitted alpha is 4e-5, then -3e-2), but every step back is exact, G keeps
decreasing, and it reaches -525.97233 by iteration 100.

Mechanism (toy model `efficiency_survey/reversibility/step_reversibility.py`): a trial step C1 = ortho(C + a d) rot1
followed by the step back C2 = ortho(C1 - a d L1 rot1) rot1^H gives C2 = C U with U = L1 rot1 L2 rot1^H. With a
triangular (Cholesky) L this U is unitary but not the identity whenever the auxiliary-Hamiltonian rotation rot1 is
not trivial: |U - I| = 1.3e-6, 1.3e-4, 1.2e-2 for a/|C| = 3e-3, 3e-2, 0.3 (energy error up to 1e-3 for distinct
fillings). With symmetric L = (C^OC)^-1/2 it is the identity to 1e-15 at every step size, and without the Haux
rotation both are exactly reversible (a unitary triangular matrix is diagonal). In Haux mode rot1 is large whenever
Haux eigenvalues are close, so this is the normal case for metals. This is the mechanism patch 02 was written for; the
new evidence is that with it the CPU run reaches the double-precision floor (|dG| 2e-10 on a 526 Ha quantity; the
1e-11 criterion is below roundoff) with step failures only at the 1e-10 level, while without it the run cannot be
converged below ~1e-7 in 200 iterations and projection reuse plays no role (base_prof shows the same jump).

Caveat on the final energies: the three runs end in different electronic states. Final nElectrons 108.3129 / 108.3179
/ 108.2702 and total moment +0.063 / +0.048 / +0.213 (HEAD / base_prof / Lowdin; the surface-O moment is 0.80 / 0.78 /
0.92), with G within 5e-6 Ha of each other. The Lowdin run's stationary point (|grad| 8e-9) lies 4.7e-6 Ha above
base_prof's iteration-200 value, which was still descending at 1e-7 per iteration. The descent from the plateau is a
slow collective mode (charge exchange with the reservoir plus moment) along which G is flat to ~1e-6 Ha, so neither
|dG| nor |grad| discriminates between these states and a "roundoff" stop can land anywhere along it. The
reversibility conclusion rests on the iteration-38 signature, not on the final energies; which of the three states
is the minimum is open, and this flat mode, not per-iteration cost, is what the production tails are made of.

Consequence for the GPU floor: the GPU runs show the same signature (undo does not restore G, failed steps at alpha
8e-6) but patch 02 did not remove it there (test B: 12 vs 9 failed steps, max increase 4.4e-7 vs 2.0e-7), so there is
an additional GPU/MPI-specific source. Suggested GPU tests, cheapest first: (a) `caseA_tight`-style long run (1e-11,
200 iterations) with p1 and p12 at np=1 (all states on one GPU) and np=4; (b) the same without the fluid;
(c) `fluid-minimize knormThreshold 1e-14`. If (a) at np=1 converges with p12, the floor is in the replicated or
distributed real-space work (the XC split above is one such rank-dependent reduction).

### Cheap non-equivalent wins already identified

- Cap the initial vacuum pre-solve in `elecFluidMinimize` (it only provides a starting guess but runs with the
  production `energyDiffThreshold`; 100 iterations and 12% of the first job in the production log).
- A looser `energyDiffThreshold` (~1e-7) for 1e4-Ha systems once the floor is understood (16-20% of iterations).

## Local CPU builds on a small solvated molecule (2026-10-29)

Written by a research-side session choosing the local build for the Patty reference-molecule opts and vibs
(`~/vs/research/Patty/gen/run_local_mol_{opts,vibs}.py`). No code was changed.

### Setup

- **System:** glycinate NH2CH2COO⁻ (9 atoms, charge −1), single point from its relaxed geometry, Patty settings.
  - GBRV 20/100 Ha, PBE-D3BJ, z-spin, Fermi 0.001, `coulomb-interaction Isolated` in a 37.79-bohr (20 Å) cube.
  - LinearPCM CANDLE with 1 M K+/F-, `energyDiffThreshold 1e-9`, `fluid-minimize knormThreshold 1e-13`.
  - nbasis 230k, nBands 24, 2 states.
  - Fresh LCAO start, so each run is a vacuum pre-solve followed by the solvated minimization.
- **Builds** (all Release with `EnableProfiling`, Accelerate BLAS, no MKL):
  - `build_merged`: upstream b3c54e1d + projection reuse.
  - `build_eff`: branch HEAD `a37d055f`, i.e. `build_merged` + the D3 pair-sum speedup.
  - `build_eff2`: `build_eff` + the uncommitted survey prototypes in the `jdftx_eff` clone (zgemm thread threshold
    1e6 MACs; XC no-split, which does nothing on CPU at np = 1). Its banner says `a37d055f` anyway.
- **Method:** one build after another at `-c 6` on the M5 Pro (6 + 12 cores, 64 GB), in normal interactive use, so
  differences of a few percent are noise.

### Results

| Build | Wall (s) | Vacuum: iterations, s/iter | Solvated: iterations, s/iter | Final F (Ha) |
|---|---|---|---|---|
| `build_merged` | 658.6 | 32, 5.75 | 17, 20.38 | −56.102161447 |
| `build_eff` | 676.8 | 32, 6.06 | 17, 20.57 | −56.102161447 |
| `build_eff2` | 664.4 | 32, 5.90 | 17, 20.32 | −56.102161447 |

- All three are the same within the noise (2.7% spread), with identical iteration counts. The final F agrees between
  builds to 1e-9 Ha, and with the original GPU ion opt (−56.102161448).
- Peak memory is 11.7 GB (`MEMUSAGE Total`) in every build. Nearly all of it is grid arrays: `ScalarField` 7.95 GB and
  `ScalarFieldTilde` 4.94 GB, against `ColumnBundle` 0.75 GB.
- `user/real` is 4.8 at `-c 6` (about 80% thread efficiency). Other thread counts were not measured.
- Profile of `build_merged` (the other builds are within ±4% per entry), out of 659 s wall:

  | Entry | Total (s) | Share |
  |---|---|---|
  | `EdensityAndVscloc` (includes the LinearPCM solves) | 347 | 53% |
  | `augmentDensityGridGrad` + `augmentDensityGrid` | 95 + 73 | 25% |
  | `Idag_DiagV_I` + `diagouterI` | 45 + 40 | 13% |
  | `ExCorrTotal` | 27 | 4% |
  | `Y1^Y2` + `Y*M` | 21 + 13 | 5% |
  | `VanDerWaalsD3::energyAndGrad` | 0.07 | — |

### Conclusions

- For isolated molecules in large boxes, neither the kept changes (projection reuse, D3 speedup) nor the
  `build_eff2` prototypes matter. The wavefunction products they target are about 5% of the run, and D3 costs
  nothing at 9 atoms. Local production uses `build_eff`: it is the committed HEAD, and its banner matches its code.
- The cost is in grid work, so different targets matter here than for the slabs:
  - **Fluid solves.** A solvated iteration costs 3.5x a vacuum one (20.4 vs 5.8 s). The difference is mostly the
    LinearPCM solve inside `EdensityAndVscloc`, which takes 14-34 CG iterations per call at `knormThreshold 1e-13`.
    Whether 1e-13 is needed for molecule free energies (and vibrations) is untested. It overlaps with the
    GPU next step of varying `knormThreshold`.
  - **Vacuum pre-solve.** It took 32 iterations, about 185 s, or 28% of this single point. That makes the existing
    "cap the vacuum pre-solve" item much bigger for small molecules than for the slab (12%). Opts and vibs restart
    from saved wavefunctions, so the cap would mainly help first runs.
  - **USPP grid augmentation.** It took 25% for only 9 atoms (0.16-0.38 s per call). Worth checking whether
    `augmentDensityGrid(Grad)` scales with the box (full grid) or with the atoms (spheres around them).
  - **Memory.** It scales with the box volume, not with the molecule. In production it caps local concurrency at
    2 solvated molecules, not 3-4 (see the follow-up below).

### Follow-up: concurrency in production (2026-10-29)

These numbers come from the Patty production runs (`build_eff`, runner logs in
`/Volumes/Ext4TB/scratch_backup/perl/Patty/PBED3BJ_20Ha_addU/mols/run_local_mol_{opts,vibs}.log`).
- **Footprint.** A solvated ion opt reaches ~15 GB per process (macOS `top` MEM, including compressed pages), against
  11.7 GB `MEMUSAGE` in the single point above. Gas-phase calcs and the vacuum pre-solve take ~4 GB.
- **Three at once does not fit.** With three solvated opts at `-c 5` plus normal desktop use, the 64 GB machine reached
  25 GB of compression and a full 6 GB swap. The solvated SCF slowed from ~20 s/iter (when gas jobs were alongside)
  to 45-55 s/iter. That is worse in total than one job at `-c 12`, which ran 11.6 s/iter (a glycinate vib on
  `build_merged`, 2026-10-26).
- **Two at once fits.** After one job was stopped, the other two (still `-c 5`) dropped to 20-28 s/iter within minutes,
  with 10 GB of compression and 7.7 GB free. Production now runs 2 x `-c 8`. In five solvated vibs (3-8 atoms), the
  median was 13.8-16.6 s/iter, about 3.5 min per displacement. That is 0.13 iter/s in total, against 0.086 for one
  job at `-c 12`. The `run_threads.sh` test below would settle the best split.
- **Profiling builds and output checks.** `EnableProfiling` builds print the PROFILER and MEMUSAGE tables after
  `Done!`. A script that tests whether `out` ends with `Done!` will flag every finished run as failed. The Patty vib
  runner did this, and it now searches for a `Done!` line instead.
- **Memory is the local throughput limit for molecules.** Shrinking the grid arrays (`ScalarField` 7.95 GB,
  `ScalarFieldTilde` 4.94 GB in a 20 Å box) would let a third job run.
- **Vacuum pre-solve, worst case.** For CO3²⁻ (charge −2, `Isolated`, 20 Å box), the vacuum pre-solve is an unbound
  dianion (mu = +0.43 Ha). It crept through 441 iterations to `energyDiffThreshold 1e-9` (~1 h under 3-way load),
  only to serve as the solvated start. This strengthens the "cap the vacuum pre-solve" item for charged molecules.
- **Vib cost.** perlstuff `vib2.py` does one-sided displacements, 3N + 1 SCFs. Each one restarts from the reference
  wavefunctions and takes 11-12 iterations to `energyDiffThreshold 1e-9` when solvated (8-12 in gas). Fluid
  solve cost and the convergence threshold therefore dominate the vibs as well.

### Reproducing

`~/Desktop/jdftx_hold/patty_mol_build_benchmark/`:
- `bench.in` and `run_bench.sh [threads] [builds...]` run the builds one after another.
- `summarize.py` prints the tables from `out_<build>_c<n>` and `time_<build>_c<n>`.
- `bench_short.in` and `run_threads.sh <build> <threads...>` (12 vacuum + 12 solvated iterations) are a ready
  thread-scaling test, not yet run.

## DFT-D3 speedup `a37d055f`: Perlmutter GPU results (2026-10-29)

Written by the Perlmutter session (NERSC, A100 nodes). No code was changed; this measures `a37d055f` as pushed. The
timings in its commit message (250-atom slab, one call: 16.8 s → 0.017 s at `-c 8`, 0.34 s → 0.098 s at `-c 1`) are Mac CPU
numbers and stand as they are; below is what the same change does on Perlmutter GPUs in the production layout.

### Setup

- **Builds.** Old = the `74dd6931` production builds in `~/jdftx_compile_test/ben`: `build-gpu-eff-74dd6931` (CUDA 13.2) and
  `build-gpu-managed-eff-74dd6931` (`CudaManagedMemory=yes`, CUDA 12.9 + cray-mpich 9.0.1, RPATH-pinned). New = `4fe71c77`
  built with the same two recipes into `/pscratch/sd/b/beri9208/claude_scratch/eff_d3/build-gpu-{,managed-}4fe71c77`
  (CMake caches identical apart from the loaded conda module version; `ldd` NEEDED lists identical). All builds have
  `EnableProfiling`; the timed quantity is `PROFILER: VanDerWaalsD3::energyAndGrad`, which covers only the vdW routine.
- **Layout and environment,** as in production: 1 node, 4 ranks × 32 threads, one A100 per rank, `--gpu-bind=none`,
  `JDFTX_GPU_SYNC_MPI=1`, `JDFTX_MEMPOOL_SIZE=36000`, `jdftx/gpu` module loaded, binaries called by full path (`ldd` showed
  each build's own `libjdftx_gpu.so`). Old and new ran back to back on the same node. Debug QOS, 40GB nodes for the managed
  builds, an 80GB node for the default builds.
- **Method:** lj-override, as in CLAUDE.md "Testing and benchmarking the van der Waals module": the production `in` reduced to
  lattice, species, ions, `elec-ex-corr`, `van-der-waals`, plus `elec-cutoff 5`, `lj-override 105.83544` (= the production
  200 bohr D3 cutoff) and `ionic-minimize nIterations 1` (energy and forces; the number of D3 calls depends on the line
  search), `lattice-minimize nIterations 1` for stress. Inputs: the 224-atom `310_fixed_1x3_trim_Ni_B/CH2COO/p07V/Mn_Mn-O_O-bi`
  production slab (D3BJ; also run with `van-der-waals D3` for zero damping), the 255-atom `NH2CH2COO Ni_Mn-O_O-bi`
  adsorbate (D3BJ), and the NO3H molecule (D3). Generator: `eff_d3/inputs/gen_vdw_input.py`.

### Time per `VanDerWaalsD3::energyAndGrad` call

| case (D3BJ unless noted) | managed build, 40GB node: old → new | default build, 80GB node: old → new |
|---|---|---|
| 224-atom slab, `ionic-minimize 1` (13 calls) | 41.3 → 0.26 s | 19.3 → 0.20 s |
| 255-atom slab, `ionic-minimize 1` (5 calls) | 53.4 → 0.29 s | 24.9 → 0.25 s (3.1 s at 1 thread per rank) |
| 224-atom slab, `lattice-minimize 1` (stress, 10 calls) | | 19.6 → 0.23 s |
| 224-atom slab, D3 zero damping (3 calls) | | 19.1 → 0.22 s |
| NO3H molecule, D3 (3 calls) | | 0.028 → 0.014 s |

- Old and new agree to ≤ 5e-15 Ha in the total energy and EvdW, ≤ 1e-14 Ha/bohr in every force component; the printed
  stress tensors and coordination numbers are identical. Checked for D3BJ, D3 zero damping and the molecule
  (`eff_d3/compare.py`).
- The old code was 2.1x slower per call on the managed build than on the default build (41 vs 19 s). The new code costs
  about the same on both (0.20-0.29 s), so on the managed 40GB production setup the saving is ~40-53 s per call.
- Per call, the new code is about 10x slower here than the Mac's `-c 8` figure (0.20-0.29 s vs 0.017 s), and 3.1 s at 1 thread
  per rank (Mac `-c 1`: 0.098 s). Not investigated; the profiled time includes the MPI reductions and the per-call
  `report()` output, and these slab cells (16 × 59 × 70 bohr) need many neighbor cells for the 200 bohr cutoff. At 0.3 s
  per call this no longer matters for production.

### End-to-end on the production configuration

Same node type and environment as production (40GB node, managed builds), the production `in` of the 224-atom
`Ni_B/CH2COO/p07V/Mn_Mn-O_O-bi` run restarted from its state files (wfns, fillings, eigenvals, fluidState of
2026-10-29 10:00, copied to `eff_d3/state224/`), 10 ElecMinimize iterations with `energyDiffThreshold 0`, then forces.
Old = `build-gpu-managed-54315e95` (the current `$JDFTx_GPU`), new = managed `4fe71c77`. Job 60503566.

| | old: managed `54315e95` | new: managed `4fe71c77` |
|---|---|---|
| s/iter, steady iterations | 29.8-30.2 | 19.1-19.2 (−36%) |
| s/iter mean over 10 (incl. one step-back iteration) | 34.1 | 22.3 |
| D3 per call (3 calls) | 42.3 s | 0.27 s |
| job wall (setup + 10 iterations + forces) | 521 s | 285 s |
| MEMUSAGE total per GPU | 42.96 GB | 36.26 GB |

Both runs took the same path (identical alpha, "Predicted alpha/alphaT" and "Step increased G" events). G differs by 5e-11 Ha at
iteration 1 and 5e-8 Ha at iteration 10 (floating-point order of projection reuse, while G is still falling by ~1e-2 Ha per
iteration); final forces differ by 4e-7 Ha/bohr.

### Estimate for the current production runs (Perlmutter)

Four recent 5 h jobs (224-228-atom slabs, PBE+U+D3BJ+CANDLE, `ionic-minimize`, managed `54315e95` on 40GB nodes): ElecMinimize
iterations are 81-86% of wall time at 33.5-36.9 s/iter (414-437 iterations, 6-9 ionic steps per job); D3 is 14-20 calls at ~40 s
(~4%); the rest (LCAO, forces, the 28 GB `wfns` dump per ionic step) ~13%. With −36% per iteration and D3 → 0, switching
`$JDFTx_GPU` from `managed-54315e95` to a `4fe71c77` build should cut wall time per ionic step by about a third (~1.5x
throughput, ~10 → ~14-15 ionic steps per 5 h job). About 5 points of that are the D3 change; the rest is projection reuse,
which the `74dd6931` production builds already have.

### Open (Perlmutter)

- The peak of 36.3 GB with projection reuse suggests the *non-managed* build might now fit these slabs on 40GB nodes. In Test A
  (80GB, 218 atoms) the default build did 11.9 s/iter, vs 19.2 s/iter for managed on 40GB here, so it could be much faster if the
  allocation cache does not push it over 40 GB. Untested.
- Production builds of `4fe71c77` (both flavors) in `~/jdftx_compile_test/ben`, then the `$JDFTx_GPU` switch by the user.

### Reproducing (Perlmutter)

`/pscratch/sd/b/beri9208/claude_scratch/eff_d3/`: `README.md` (results), `build.sh default|managed`, `inputs/` with
`gen_vdw_input.py` and `gen_e2e_input.py`, `job_vdw.sh` (old vs new per case, same node), `job_e2e.sh`, `compare.py vdw
default|managed` / `compare.py e2e`, outputs in `runs/`, Slurm logs in `logs/`. Jobs 60503563 (vdw, managed, 40GB),
60503565 (vdw, default, 80GB), 60503566 (e2e, 40GB).

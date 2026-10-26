# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> **This file lives only on the orphan branch `claude-md`.** Never commit it to `master` or any code branch, so it
> can't end up in PRs from this fork. Each clone ignores it locally (`/CLAUDE.md` in `.git/info/exclude`, not
> `.gitignore`), and you must never add it with `git add -f`. Claude Code agents on several machines share it
> through `origin/claude-md`, so keep it in sync:
> - **Before reading it, update from the remote:** `git fetch origin claude-md:claude-md && git load-claude-md`
>   (without the alias: `git show claude-md:CLAUDE.md > CLAUDE.md`). If the fetch is refused as non-fast-forward,
>   the local branch has unpushed saves: push or reconcile them first. Loading overwrites the working copy, so
>   save any local edits before you load.
> - **After every edit, save and push:** `git save-claude-md "what changed" && git push origin claude-md`.
>   `save-claude-md` commits the working-tree file to `claude-md` without switching branches. Both aliases are
>   local to each clone's `.git/config`; if `git save-claude-md` says "not a git command", run the setup below.
>
> **First-time setup in a new clone** (run once from the repository root; safe to re-run). Before setup you can
> read this file with `git fetch origin claude-md && git show origin/claude-md:CLAUDE.md`.
> ```bash
> git fetch origin claude-md:claude-md
> grep -qxF /CLAUDE.md .git/info/exclude || echo /CLAUDE.md >> .git/info/exclude
> git config alias.load-claude-md '!git show claude-md:CLAUDE.md > CLAUDE.md && echo "CLAUDE.md restored from claude-md"'
> git config alias.save-claude-md '!f() { parent=$(git rev-parse --verify -q refs/heads/claude-md) || { echo "no local claude-md branch (git fetch origin claude-md:claude-md)"; return 1; }; blob=$(git hash-object -w CLAUDE.md) || return 1; if [ "$blob" = "$(git rev-parse claude-md:CLAUDE.md)" ]; then echo "CLAUDE.md unchanged; nothing to save"; return 0; fi; tree=$(printf "100644 blob %s\tCLAUDE.md\n" "$blob" | git mktree) && commit=$(git commit-tree "$tree" -p "$parent" -m "${1:-Update CLAUDE.md}") && git update-ref -m "save-claude-md" refs/heads/claude-md "$commit" "$parent" && echo "claude-md -> $(git rev-parse --short "$commit"); push with: git push origin claude-md"; }; f'
> git load-claude-md
> ```
> Check it: `git status` shows no `CLAUDE.md`, and `git save-claude-md` prints "unchanged". Saving makes a commit, so
> it needs a git identity. If `git config user.email` prints nothing, set `user.name` and `user.email` first
> (`benrich37`, `benjamin.rich@colorado.edu`, as on earlier `claude-md` commits).

## Repository layout

The CMake source root is `jdftx/`, not the repository root. The repo root also contains `fluid1D/`
(a separate 1D classical-DFT fluid code with its own CMakeLists), a `Dockerfile` (Ubuntu build of upstream
master), and `generate_tarball.sh`. This is a fork: `origin` = benrich37/jdftx, `upstream` = shankar1729/jdftx.

## Build

Out-of-source builds only (`build/` is gitignored):
```bash
mkdir -p build && cd build
cmake ../jdftx [-D EnableCUDA=yes] [-D EnableMKL=yes] [-D CMAKE_BUILD_TYPE=Release]
make -j
```
- Library sources are collected with `FILE(GLOB ...)` from `core/ fluid/ electronic/ perturb/ commands/`.
  After adding or removing a `.cpp`/`.cu` file, rerun `cmake` or it won't be compiled.
- The build produces `libjdftx` plus thin front-end executables `jdftx`, `phonon` and `wannier`
  (`jdftx.cpp` is under 20 lines). With `EnableCUDA`, the same sources are compiled a second time into
  `libjdftx_gpu` with `GPU_ENABLED` defined, plus the `.cu` kernels, giving `jdftx_gpu`, `phonon_gpu` and
  `wannier_gpu`.
- All CMake options and their meaning are in `jdftx/doc/compiling/Customization.dox`. GPU-relevant ones:
  `CudaAwareMPI`, `PinnedHostMemory`, `CudaManagedMemory`, `EnableCuSolver`.
- `make doc` builds the Doxygen manual (command reference pages are generated from command source; see below).

## Cluster builds (HPC; e.g. Perlmutter, Anvil)

Read this before compiling on any cluster. The Perlmutter recipe below is known-good (verified 2026-10-20:
fork master `54315e95`, PrgEnv-gnu, CUDA 13.2, A100, `make -j16` on a login node takes ~9 min). Other clusters
need their own modules and paths, but the same structure.
```bash
module use --append /global/cfs/cdirs/m4025/Software/Perlmutter/modules
module load jdftx/gpu   # depends_on: gpu PrgEnv-gnu cmake craype-x86-milan cray-fftw cray-hdf5-parallel libxc
mkdir build-gpu && cd build-gpu            # sibling of the git clone
CC=cc CXX=CC CUDAARCHS=80 cmake \
    -D EnableProfiling=yes -D EnableHDF5=yes -D EnableLibXC=yes \
    -D EnableLibSci=yes -D EnableScaLAPACK=yes \
    -D FFTW3_PATH=${FFTW_ROOT} -D GSL_PATH=${GSL_ROOT} \
    -D EnableCUDA=yes -D EnableCuSolver=yes -D CudaAwareMPI=yes -D PinnedHostMemory=yes \
    -D CMAKE_LIBRARY_PATH="${LD_LIBRARY_PATH//:/;}" \
    ../jdftx-git/jdftx
make -j16
```
Translating to another cluster:
- `CC`/`CXX` must be the MPI compiler wrappers (`cc`/`CC` on Cray; `mpicc`/`mpicxx` elsewhere). If several MPI
  stacks are installed, pass the one FFTW was built against by full path (`-D MPI_CXX_COMPILER=...`,
  `-D MPI_C_COMPILER=...`). The user's Mac build needed this because Open MPI and MPICH are both installed.
- `CUDAARCHS` is the GPU compute capability (80 = A100, 90 = H100). For a CPU-only build, drop the CUDA options.
- `EnableLibSci` is Cray-only. Elsewhere, take BLAS/LAPACK from `-D EnableMKL=yes` (Intel stacks) or let CMake
  find system/OpenBLAS libraries. `EnableScaLAPACK`, `EnableHDF5` and `EnableLibXC` are optional; turn on only what
  the cluster's modules provide.
- `CMAKE_LIBRARY_PATH="${LD_LIBRARY_PATH//:/;}"` lets CMake find libraries that modules expose only through
  `LD_LIBRARY_PATH`. Run cmake, make *and jobs* with the same modules loaded. Perlmutter job scripts load no
  modules and rely on `sbatch` exporting the submitting shell's environment (the default `--export=ALL`).
- Set `-D CudaAwareMPI=yes` only if the MPI library really is CUDA-aware (Cray MPICH with
  `MPICH_GPU_SUPPORT_ENABLED=1`; Open MPI built with CUDA). Otherwise leave it OFF.
- Never modify group or shared installs (e.g. `/global/cfs/cdirs/m4025/**` on Perlmutter). Build in your own
  scratch directory.

**CUDA-aware MPI race (any GPU build with `CudaAwareMPI=yes` run with more than one MPI process):** MPI can read
device buffers before the asynchronous kernels or cuBLAS calls writing them finish. The symptom is silently wrong
energies and failed convergence: "Bad step direction", "Step increased F" scaling with alpha, and "Probably at
roundoff error limit" at |grad| ~1e-3. With 1 k-point and z-spin, any np>=2 is affected because the two spins
sit on different ranks. Confirmed on Perlmutter: NO3H in CANDLE LinearPCM gives first `LCAOMinimize: Iter: 0`
F = -57.9042 Ha at np=1 but -56.9602 at np=4. With the fix, np=4 converges to the np=1 energy within 1e-10 Ha.
- Upstream fix `544d0e78` (Sundararaman, 2026-10-02) syncs in `ManagedMemory::dataMPI()`, but **only if env
  `JDFTX_GPU_SYNC_MPI` is set**. Any value enables it, even `0`. Fork master contains it, so export
  `JDFTX_GPU_SYNC_MPI=1` in job scripts. Branch `gpu-sync-mpi-default` makes it default-on (`=0`/`no`/`off`
  disables it); that branch compiles but has not yet been run as a full build.
- Always check that `out` says `gpuInit: sync before MPI is on`. Builds older than `544d0e78` print no such line
  and ignore the variable; those, including the Perlmutter m4025 shared build as of 2026-10, need rebuilding.
- Quick check for a new GPU build: run the same small spin-polarized molecule at np=1 and np=4 (same node,
  `--gpus-per-task=1`) and compare the first `LCAOMinimize: Iter:   0` energy. With the sync on, the
  two agreed to 1.6e-5 Ha. A difference near 1 Ha means the race is present.

Running a freshly built binary:
- Call it by full path. A loaded JDFTx module may put another `jdftx_gpu` on `PATH` and its library directory
  first in `LD_LIBRARY_PATH`. Check with `ldd <build>/jdftx_gpu | grep libjdftx` that the library resolves into
  your build directory. That holds on Perlmutter because the binary carries `DT_RPATH`; a toolchain that links
  with `RUNPATH` (`--enable-new-dtags`) lets `LD_LIBRARY_PATH` win instead.
- `JDFTX_MEMPOOL_SIZE` is deprecated; use `JDFTX_CACHE_SIZE` (MB per GPU).
- Use the debug queue with short limits for tests, and ask the user before submitting any job. Perlmutter tests
  used `-C gpu -A m4025_g -q debug -N 1 -c 32 --ntasks-per-node=4 --gpus-per-task=1 --gpu-bind=none`. The
  `--gpu-bind=none` causes a harmless Cray GTL warning ("not using the same device it did during gtl_init").

### Anvil (Purdue, A100)

A working GPU build lives in `/anvil/projects/x-chm260113/Software/JDFTx/`: `anvil-env.sh` (modules),
`run-cmake-gpu.sh` (cmake), `build-gpu/`, and `examples/gold_sp/asubmit.sh` (a job script known to work). Its
`jdftx-git` clones upstream `shankar1729/jdftx` (not the fork) at `52548a96` (2026-07-13), with `CudaAwareMPI=OFF`.
Those files are the ground truth; the notes below explain why they look the way they do. This is a shared group
directory, so ask before rebuilding in place.

Environment:
- `module load modtree/gpu` comes first. Anvil's GPU software is a separate module tree, and without it the
  CUDA builds of OpenMPI are invisible. Then load, in this order: `gcc/13.3.0 cuda/12.6.1 openmpi/4.1.6-cu12.6.1
  openblas/0.3.27 fftw/3.3.10 gsl/2.7`. OpenMPI and OpenBLAS appear only after `gcc/13.3.0` is loaded. Pair the
  `-cu12.6.1` OpenMPI with `cuda/12.6.1`, not with another CUDA version.
- Anvil has no MKL under any name. Use `-D EnableMKL=no -D BLA_VENDOR=OpenBLAS`, and don't copy MKL settings
  from Bridges.
- The GPUs are A100s: `CUDAARCHS=80`, `-D CUDA_ARCH=compute_80 -D CUDA_CODE=sm_80` (Bridges used 70).
- Login nodes have the CUDA toolkit but only a stub driver. Compiling works there; running `jdftx_gpu` doesn't.

Build: `source anvil-env.sh`, then run `bash ../run-cmake-gpu.sh && make -j16` from inside `build-gpu/`.
- Source `anvil-env.sh`; don't execute it. The `.cu` targets find GSL headers only through `CPATH` from the gsl
  module, while the C++ targets get them from CMake. `run-cmake-gpu.sh` loads the same modules, but they are
  unloaded when its subshell exits. `make` then fails with `gsl/gsl_integration.h: No such file or directory` on
  `Coulomb.cu`/`BlasExtra.cu` while every `.cpp` file compiles.
- `run-cmake-gpu.sh` also sets `CMAKE_CUDA_HOST_COMPILER="$(which g++)"` (otherwise nvcc can fall back to the
  system gcc 8) and passes explicit GSL and FFTW `-I` flags in `CMAKE_CUDA_FLAGS`. Keep both.
- A changed CUDA host compiler does not take effect in an existing CMake cache. After any toolchain change,
  `rm -rf build-gpu` and reconfigure.

Running:
- Every job script must `source anvil-env.sh`. The binary's RPATH covers gcc and OpenMPI but not OpenBLAS, FFTW, GSL
  or the CUDA libraries. Without the modules, `ldd` reports them as "not found".
- Launch with `mpirun -np N <path>/jdftx_gpu`, as `asubmit.sh` does. A bare `srun` fails at `MPI_Init` with
  "OPAL ERROR: Unreachable in ext3x_client.c". This OpenMPI was built `--without-pmi` against PMIx 5.0.3, while
  Slurm offers only `pmi2` and `pmix_v2` (`srun --mpi=list`). Don't "simplify" the launch to `srun`.
- Allocations use the `-gpu` account suffix: `-A chm260113-gpu`. The test partition is `gpu-debug`.
- Thread count: the binary reads `-c`, then `SLURM_CPUS_PER_TASK`, then `OMP_NUM_THREADS`. With none of these set,
  it splits the node's cores among its ranks (`core/Util.cpp`), so set `--cpus-per-task`. It does not read
  `JDFTX_NUM_PROCS` in either capitalization: the `JDFTx_NUM_PROCS=1` in `asubmit.sh` did nothing, and gold_sp
  ran 2 ranks × 64 threads.
- `MPICH_GPU_SUPPORT_ENABLED` is an MPICH/Cray variable and does nothing under OpenMPI.
- This build predates `544d0e78`, so it ignores `JDFTX_GPU_SYNC_MPI`. It is also built without `CudaAwareMPI`,
  which the race described above requires. Pull `544d0e78` or later before rebuilding with `CudaAwareMPI=yes`.

## Tests

Regression tests are ctest-driven bash harnesses in `jdftx/test/` (see `jdftx/test/README`):
```bash
cd build
make test                                         # all 11 tests
ctest -R vibrations --output-on-failure           # single test
make testresults                                  # print summary of previous runs
make testclean                                    # REQUIRED before rerunning: runs whose .out ends in "Done!" are skipped
JDFTX_LAUNCH="mpirun -n %d" make test             # MPI (%d -> nProcs from each test's sequence.sh)
JDFTX_LAUNCH="" JDFTX_SUFFIX="_gpu" make test     # GPU binaries
```
Per-test results are in `build/test/<name>/results` and `summary`. Each test directory has `*.in` files
(which include shared files via `${SRCDIR}/...`), `sequence.sh` (run order and nProcs) and
`checkResults.sh` (prints `nChecks`, then `obtained expected tolerance description` lines).
There is no linter or unit-test framework.

## Running

`jdftx -i in -o out` (`-c N` threads, `-n` dry run: parse and initialize only, `-t` prints every command's
documentation and syntax). **`-o` appends to an existing output file** unless `-d`/`--no-append` is given.
Inputs are command files; `$VAR` in filenames expands per dumped variable (`dump-name`).
Internal units are atomic (Hartree, bohr); see `core/Units.h`.

## Architecture (big picture)

- **`Everything`** (`electronic/Everything.h`) is the container for a whole calculation: `gInfo` (grid/FFT),
  `iInfo` (ions and species), `eInfo` (electronic states and fillings), `eVars` (wavefunctions, density,
  potentials), `coulomb`, `symm`, `dump`, fluid solver, `vibrations`, control parameters. Commands write into
  it, `Everything::setup()` initializes it, and `electronic/run.cpp::run()` chooses the driver:
  dry run → vibrations (bypasses ionic minimization) → perturbation solver → lattice/ionic/electronic minimize.
- **Commands** (`commands/*.cpp`) register themselves: each is a `struct CommandX : public Command`
  with a file-scope instance, so there's no central list to edit. `format`/`comments` strings *are* the user
  documentation, shown by `jdftx -t` and in the Doxygen manual; keep them accurate when changing behavior.
  `process()` parses via `ParamList::get(...)` and rejects bad input by throwing `string("...")`.
  `printStatus()` must print a command line that reproduces the parsed state (it is echoed into every output file).
- **CPU/GPU dual implementation:** the numerical kernel for one element is an inline `__hostanddev__`
  `foo_calc(...)` in a `*_internal.h` header. It is looped over with `threadLaunch` in `foo()` (CPU, `.cpp`)
  and launched as a CUDA kernel in `foo_gpu()` (`.cu`). Callers use `callPref(foo)(..., x.dataPref())`,
  which resolves to `foo_gpu` in GPU builds (`core/BlasExtra.h`).
- **Memory:** `ManagedMemory` (`core/ManagedMemory.*`) owns all arrays (`ScalarField`, `ColumnBundle`,
  `matrix`, ...) and moves data between host and device on demand (`data()` vs `dataGpu()`/`dataPref()`).
  Allocations go through a cache (`JDFTX_CACHE_SIZE`), so `cudaFree` no longer provides implicit syncs.
  With `CudaAwareMPI`, MPI calls take device pointers via `dataMPI()`; GPU work must be finished before MPI
  reads the buffer (the `JDFTX_GPU_SYNC_MPI` env var in upstream, on by default in branch `gpu-sync-mpi-default`).
- **Parallelization:** pthreads within a process; MPI distributes electronic *states* (k-points × spins)
  over processes (`eInfo.qStart/qStop`, `eInfo.isMine(q)`), while real-space grids and fluid state are
  replicated on every process and combined with `allReduce`. Code after a reduction runs redundantly on all
  ranks, so every rank must make the same collective calls. The templated CG/L-BFGS minimizers in
  `core/Minimize*.h` call `sync()` on scalars for this reason. `die()` must be reached by all ranks;
  use `die_alone()` where that can't be guaranteed. Use `logPrintf` for output (only the head process logs).
- **Physics modules:** `electronic/` (Kohn-Sham DFT, pseudopotentials in `SpeciesInfo`, minimizers, `Dump`,
  `Vibrations`), `fluid/` (solvation hierarchy under `FluidSolver`: `PCM` → `LinearPCM`, `NonlinearPCM`, `SaLSA`,
  plus classical-DFT `FluidMixture`), `perturb/` (variational perturbation theory), `phonon/` and `wannier/`
  (extra executables and commands), `core/` (grids, FFT, Coulomb kernels incl. truncation and embedding, BLAS
  wrappers, minimizers).
- Developer overview with an architecture diagram: `jdftx/doc/devguide/DevGuide.dox`.

## Branch notes: vibrations checkpointing (`vib_chk`)

On `vib_chk`, `vibrations` saves each configuration's forces and dipole as `grad<i>`/`Pel<i>` (named via
`dump-name`) and **reads any of these files that already exist instead of recomputing them**. Options
`iConfiguration`/`nConfigurations` compute a range without analysis. When testing, always use a fresh
directory: stale files from an earlier build or input are silently reused, and `make testclean` does not
delete them.

## Branch notes: electronic minimizer efficiency (`efficiency_testing`)

Rationale, tests and A/B timings are in `efficiency_notes.md` on that branch (patches in `efficiency_patches/`).
Merged with upstream through b3c54e1d (2026-10-26); the branch's only code change is projection reuse
(ctest 11/11 on CPU and GPU; upstream merge checks in `efficiency_notes.md`):
- **Projection reuse** (`bd608a26`): with ultrasoft psps, O(C) = C·detR + V·Q·(V^C), and V^C products
  (nProj×nBands×nbasis) dominate an ElecMinimize iteration because nProj > nBands. `ElecVars::VdagC[q]` must
  always hold the projections of the current `C[q]`: `dagO()` builds overlaps from cached projections,
  `applyHamiltonian(..., orthoGrad=true)` folds `-O(C)·Hsub` into `projectGrad`, and `ElecGradient::VdagC`
  (set in `ElecMinimizer::constrain`, propagated by `axpy`/`*=`, cleared by `init`/`randomize`) lets `step`
  update projections linearly via `orthonormalize(q, rot, projected=true)`. Code that modifies `C` must keep
  `VdagC` consistent or clear the direction's cached projections.
- **Löwdin orthonormalization** (`f047b9b0`): tested, no performance gain, **removed** (revert `5c09112b`).
  `orthoMatrix()` is Cholesky (non-symmetric); with per-band fillings it mixes occupied and empty states, so linmin's
  trial step plus step back is not reversible (jumps up to ~1e-2 Ha seen in a CPU test system). `invsqrt()` fixed that
  but needed more iterations on GPU. Recoverable from `git stash list`, the commit, or
  `efficiency_patches/not_adopted/02_lowdin_orthonormalize.patch`. `ElecVars::orthonormalize` uses Cholesky again.
- Profile with `-D EnableProfiling=yes` (per-function `PROFILER:` lines at the end of the output). Compare
  builds side by side at equal thread counts; iteration counts, not wall time, are the fair metric when runs
  share the CPU.

#### GPU results (Perlmutter, 2026-10-24/25)

Full tables and discussion: "GPU results" in `efficiency_notes.md` on `efficiency_testing` (commit `f4466321`).
- **Builds** (Perlmutter recipe, cuSolver on, CudaAwareMPI, profiling): base `54315e95` (master), p1 `bd608a26`,
  p12 `fe9df38d`.
  - They live in `/pscratch/sd/b/beri9208/claude_scratch/efficiency_gpu/build-*`. Raw outputs are in `.../efficiency_gpu/runs/{A,B,C}_{base,p1,p12}/`,
    together with `analyze.py`, `analyzeC.py` and `job.sh` in the parent directory.
  - The jobs unload `jdftx/gpu` (the user's profile loads it and it puts the shared m4025 build on the paths) and load only its
    dependencies. Every run had `SYNC_OK` and `ldd` showed each build's own library.
- **System and layout:** the production Ni_B p07V slab (218 atoms, nBands 1528, 4 states), restarted from its final state
  files at the `in` positions. 1 node, 4 A100-80GB, `-c 32 --ntasks-per-node=4 --gpus-per-task=1 --gpu-bind=none`, account
  `m5414_g`.
- **GPU regression suite:** 11/11 for both p12 and base.

| Test | base | p1 | p12 |
|---|---|---|---|
| A: s/iter (20 iters) | 16.42 | 11.87 (−27.7%) | 11.96 (−27.2%) |
| A: max \|G−G_base\| | — | 2.5e-10 Ha | 0.14 → 7e-5 Ha (different path) |
| A: Y1^Y2 / Y*M / augmentOverlap calls | 569 / 625 / 336 | 325 / 297 / 8 | 325 / 297 / 8 |
| A: MEMUSAGE total | 42.43 GiB | 35.81 GiB | 35.81 GiB |
| B: iters, stop | 52, roundoff | 52, roundoff | 60, roundoff |
| B: "Step increased G" n, max | 9, 2.1e-7 | 9, 2.0e-7 | 12, 4.4e-7 |
| B: iters after \|G−Gf\|<1e-6 | 9 | 9 | 12 |
| C: ElecMin iters (4 ionic steps) | 216 (6 minimizations) | 191 (5) | 231 (5) |
| C: wall time | 76.5 min | 50.6 min | 64.6 min |
| C: "Step increased G" events | 45 | 38 | 82 |

**Conclusions:**
- **Projection reuse: keep.** It reproduces base at GPU roundoff and is 27% faster per iteration (25% faster per
  iteration including forces/D3). It also uses 6.5 GiB/GPU less memory, one wavefunction-sized array.
- **Löwdin: no benefit on GPU.** Its per-iteration cost is negligible, but it does not remove the noise floor, and it needed
  more iterations and more failed steps in both B and C. These are single runs, but they agree. Leave it out of production builds.
- **The ~1e-7 Ha floor is still there with the sync on, so it is not the MPI race and not Cholesky.** After one trial
  evaluation, every later evaluation is about 1.2e-7 Ha higher, even at alpha 8e-6, and "Undoing step" does not remove it. Leading
  suspect (unproven): the warm-started LinearPCM solve. About 16-20% of iterations are spent within 1e-6 Ha of the final G.
- No new code commits. cuBLAS zherk was no faster than zgemm. +U orbital caching would save about 3% for about 1.6 GiB, the wrong
  trade on 40GB nodes.

**Next steps:**
1. Rerun B_p1 with `fluid-minimize knormThreshold 1e-14` (and 1e-9). If the floor scales with the threshold, it is the fluid.
   About 15 min each in debug.
2. Then use a looser stop for 1e4-Ha systems (`energyDiffThreshold` about 1e-7) to cut the tail.
3. Try p1 for the 255-atom Ni_B adsorbates on 40GB nodes (45.4 GiB on base, so about 39 GiB expected).
4. Production builds: change 1 on top of master, without change 2.
5. The 9 `wfns` copies were deleted on 2026-10-26 (with the user's OK). The run directories keep `in`, `out`, the small state files and `run_info.txt` (502 MB). A rerun needs `wfns` copied again from the production directory.

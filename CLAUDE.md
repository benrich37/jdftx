# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

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
> - **STANDING PERMISSION: push Claude notes without asking** (granted by the user 2026-10-29; an exception to the
>   usual "ask before pushing" rule). Agents on different machines communicate through these notes, so pull before
>   editing and push right after saving, as part of normal work. Don't stop to ask, and don't wait on a permission prompt.
>   - Always for `claude-md` (this file).
>   - Also for other Claude note-taking markdown on code branches (e.g. `efficiency_notes.md` on
>     `efficiency_testing`), but only if every commit in the push touches nothing except such notes files. Check with
>     `git log --stat origin/<branch>..<branch>` first. If the push would also carry code commits, or the branch is
>     `master`, ask.
>   - Everything else follows the usual rules: push only Claude-owned branches, never `master`, and ask otherwise.
> - **Notes stay machine-aware** (user, 2026-10-29). Observations in these files belong to the machine they were made on
>   (Mac CPU, Perlmutter GPU, Anvil, ...). If a new measurement disagrees with an existing claim, do not change the
>   claim: leave it in place and add a note next to it saying what you saw, where and when (e.g. "Perlmutter,
>   2026-10-29: 10% here, not 30%"). Label every new result with its machine and setup. Plain typo fixes in your
>   own machine's notes are fine.
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
*Note (mbook2019, 2026-11-05; not verified on Perlmutter):* `jdftx/CMakeLists.txt` sets no default build type and no `-O` level
(checked on `b3c54e1d` and `4839c10f`), and the recipe above passes no `CMAKE_BUILD_TYPE`. On mbook2019 (gcc 15.2, CMake 4.2.3) a
plain `cmake` gives `CXX_FLAGS = -Wall -Wno-unused-result -ftemplate-depth-512`, i.e. gcc's default -O0; `-D CMAKE_BUILD_TYPE=Release`
adds -O3. Check `grep CMAKE_BUILD_TYPE <build>/CMakeCache.txt` for the Perlmutter builds: if it is empty, host-side code (e.g.
`VanDerWaalsD3`) is unoptimised, which may be related to the unexplained ~10x D3 per-call gap vs the Mac in `efficiency_notes.md`.

Perlmutter production builds since 2026-10-29 are instead compiled with the user's own module `~/modules/jdftx/gpu.lua`
through a self-contained driver, `~/jdftx_compile_test/ben/build-eff-4fe71c77.sh` (copy and edit it for a new commit).
It writes `run_env.sh` (the build's library environment) and `BUILD_INFO` into the build directory; both are required
by the channel system below. The group module directory `/global/cfs/cdirs/m4025/Software/Perlmutter/modules` is
read-only for agents.

Translating to another cluster:
- `CC`/`CXX` must be the MPI compiler wrappers (`cc`/`CC` on Cray; `mpicc`/`mpicxx` elsewhere). If several MPI
  stacks are installed, pass the one FFTW was built against by full path (`-D MPI_CXX_COMPILER=...`,
  `-D MPI_C_COMPILER=...`). The user's Mac build needed this because Open MPI and MPICH are both installed.
- `CUDAARCHS` is the GPU compute capability (80 = A100, 90 = H100). For a CPU-only build, drop the CUDA options.
- `EnableLibSci` is Cray-only. Elsewhere, take BLAS/LAPACK from `-D EnableMKL=yes` (Intel stacks) or let CMake
  find system/OpenBLAS libraries. `EnableScaLAPACK`, `EnableHDF5` and `EnableLibXC` are optional; turn on only what
  the cluster's modules provide.
- `CMAKE_LIBRARY_PATH="${LD_LIBRARY_PATH//:/;}"` lets CMake find libraries that modules expose only through
  `LD_LIBRARY_PATH`. Run cmake, make *and jobs* with the same modules loaded, unless the build carries its own
  `run_env.sh`: for Perlmutter production jobs since 2026-10-29, perlStuff sources the chosen build's `run_env.sh` for
  the `srun` launch only, so the modules loaded in the submitting shell no longer matter (channel system below). The
  legacy path (`$JDFTx_GPU` frozen at submission, resolution step 4) still loads no modules in the job script and
  relies on `sbatch` exporting the submitting shell's environment (the default `--export=ALL`).
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

### Perlmutter production builds and channels (since 2026-10-29)

Verified on Perlmutter 2026-10-29 to 2026-11-03. The calculation-side rules (how to write `psubmit.sh`, which node type,
the `jdftx_build` override file) are in the research repo's root `CLAUDE.md` (`~/vs/research`, "Submitting to
Perlmutter: JDFTx build channels"); this section is the build side.
- A production job no longer uses the executable in `$JDFTx_GPU` at `sbatch` time. The build is chosen when the job
  starts by `~/perlStuff/perlStuff/helpers/jdftx_channel.py`, called from `helpers/calc_helpers.py::get_exe_cmd` (opt4.py
  and ~40 other perlStuff scripts; `opt6*.py` and `ase_opt_self_contained/` run JDFTx in-process through `pyjdftx` and
  are NOT covered). On 2026-11-05, with the user's permission, the hook was committed as perlStuff `3121f8b`, merged
  with the user's 2026-10-21 local-run commits (`6e4a92b`, `0a64e90`) as `f15b6fd`, and pushed to `origin/master`.
  The Perlmutter working copy was fast-forwarded after a check that the opt4 production path is unaffected. As of
  `1fb5533`, perlStuff also ignores `__pycache__` and its `CLAUDE.md` describes the resolver. That file also holds the
  user's keep-synced rule for the Perlmutter working copy and the Anvil to-do list. Pushes from Perlmutter go over SSH
  (`pushurl`); the https remote has no credentials there. Every clone that pulls master (laptop, Anvil) gets the
  resolver too. Without `~/jdftx_prod`, it falls back to `$JDFTx_GPU`/`$JDFTx` unchanged.
- Resolution order, first hit wins: (1) `<job submit dir>/jdftx_build`, first non-comment word = channel name or build
  directory; (2) `$JDFTx_GPU_CHANNEL` from `psubmit.sh`; (3) `~/jdftx_prod/adopt`, mapping a frozen legacy `$JDFTx_GPU`
  path to a channel (`build-gpu-managed-54315e95 -> managed`, `build-gpu-54315e95 -> default`,
  `build-gpu-mpifix-10-20-26 -> managed`); (4) `$JDFTx_GPU` as frozen at submission (legacy, no `run_env.sh`). A failure
  in 1-3 logs `WARNING: JDFTx build resolution failed (...)` and falls through to 4, so the job runs on the OLD build.
- Channels are symlinks in `~/jdftx_prod/channels/`: `managed` ->
  `~/jdftx_compile_test/ben/build-gpu-managed-eff-4fe71c77` (CUDA managed memory, CUDA 12.9 + cray-mpich 9.0.1) and
  `default` -> `~/jdftx_compile_test/ben/build-gpu-eff-4fe71c77` (CUDA 13.2, plain device memory); both commit
  4fe71c77 of benrich37/jdftx efficiency_testing (D3 rewrite + projection reuse), promoted 2026-10-29 14:40.
  `~/jdftx_prod/channels/HISTORY` records every change. A channel is resolved once per process: a running job keeps
  its build when the channel is repointed; every job that starts afterwards, across all projects, gets the new one.
  Pin a job with `jdftx_build <build dir>` if it must keep a specific build.
- Promote: `~/jdftx_prod/bin/jdftx-promote <channel> <build_dir>` (refuses a build without `run_env.sh`, with
  unresolved libraries, or on `/pscratch`; swaps the symlink atomically). Rollback:
  `~/jdftx_prod/bin/jdftx-promote managed ~/jdftx_compile_test/ben/build-gpu-managed-54315e95` and
  `... default ~/jdftx_compile_test/ben/build-gpu-54315e95`; no resubmission needed. Build directories referenced by
  channels, adopt or any `jdftx_build` file must stay in `~/jdftx_compile_test/ben/` (never `/pscratch`); the
  54315e95 builds stay as rollback targets. Retarget one queued job with `echo managed > jdftx_build` (or `default`,
  or a build dir) in its directory.
- Which build where (pairing rule, never broken: **a job on 40GB nodes runs the managed build**): the ~220-266-atom
  slabs (nBands ~1600, nbasis ~287k) peak at ~45 GiB/GPU with `default` and OOM on 40GB A100s; with `managed` they run
  at 18.3 s per ElecMinimize iteration (220-228 atoms) and ~25 s (242-266 atoms), ~1.55x slower than `default` (was
  1.8x with 54315e95). Use `default` for anything that fits in 40 GB: `managed` is ~5x slower on small systems (NO3H
  single point 3m07 vs 38 s). 80GB + `default` is reserved for short expected waits (restart runs of ~220 atoms peak
  42.6 GiB; fresh LCAO starts of the ~255-atom systems OOM'd even on 80GB).
- Verify in `opt.iolog`: `grep "JDFTx build:\|WARNING\|exe_cmd" opt.iolog` shows the resolved executable with the
  winning step in brackets, the `BUILD_INFO` line (`commit=4fe71c77 | flavor=managed | built=...`) and the launch
  (`. <build>/run_env.sh && srun <build>/jdftx_gpu`). `run_env.sh` also exports `MPICH_GPU_SUPPORT_ENABLED=1` and
  `JDFTX_GPU_SYNC_MPI=1`. In `out` expect `Run totals: 4 processes, 128 threads, 4 GPUs` and `gpuInit: sync before MPI
  is on`. Rollout tests, the template `psubmit.sh` and the original files: `/pscratch/sd/b/beri9208/claude_scratch/build_channels/`.

### Porting the build channels to another cluster, and choosing node types (for the Anvil agent)

Written by the Perlmutter session, 2026-11-05. Everything below was built and measured on Perlmutter only. The
resolver code is on perlStuff `origin/master` (`f15b6fd`; see the next paragraph and "Where Anvil differs"), but
`~/jdftx_prod`, the promote script, the build driver and the research-repo `CLAUDE.md` are not on Anvil. So those
parts are a spec to reimplement, not files to copy. Ask the user before you patch the
live job-launch code (every queued job picks the change up when it starts), before you write into a shared group
directory, and before you submit any job.

**Why it exists.** The user rebuilds JDFTx faster than queued jobs start, and resubmitting a job loses its queue
position. `sbatch` copies the job script and (with the default `--export=ALL`) the submitting shell's environment at
submission, and Slurm can't change either for a pending job. What the job *reads when it starts* (an imported Python
module, a sourced file, a symlink) can still change. So the choice of build moved from submission to job start.

**Parts** (~150 lines in total on Perlmutter; resolution order and fallback are in the section above):
- *Immutable, versioned build directories* (`build-gpu-[managed-]<tag>-<commit>`, ~430 MB each plus a ~60 MB source
  worktree) on persistent storage, never purgeable scratch. Never rebuild into an existing directory. The build driver
  (see `build-eff-4fe71c77.sh` above) checks that the source worktree is at the intended commit and refuses an existing
  build dir. In the same shell that ran cmake and make, with the modules still loaded, it writes two files:
  `run_env.sh` (exports that shell's `LD_LIBRARY_PATH`, plus `MPICH_GPU_SUPPORT_ENABLED=1` and `JDFTX_GPU_SYNC_MPI=1`)
  and `BUILD_INFO` (`commit=`, `flavor=`, `built=`, `modulefile=`, `modules=`).
- *Channels:* symlinks `<prod root>/channels/{managed,default}` pointing at build dirs. `channels/HISTORY` is
  append-only (`<ISO time> <channel>: <old> -> <new>`) and doubles as the rollback list.
- *`jdftx-promote <channel> <build dir>`* refuses a build that has no `run_env.sh`, shows any `not found` in
  `ldd` run under its `run_env.sh`, or sits on scratch. An optional GPU dry run (`jdftx_gpu -n -i in`) checks a small
  input. The script then swaps the link atomically (`ln -sfn <build> .<ch>.new && mv -T .<ch>.new <ch>`; a bare
  `ln -sfn` unlinks before it creates, so a job starting in between would see no channel) and appends to HISTORY.
- *Resolver* `resolve_build(gpu, log_fn) -> (exe, run_env or None, reason)`. It is called by the single function that
  builds the launch command, so ~40 perlStuff scripts change with one 8-line edit there. The override file is read from
  `$SLURM_SUBMIT_DIR`. A target with no `/` is a channel name; anything else is a build dir or executable path. The
  resolver follows symlinks to a real path, caches the result per process, and logs the executable, the rule that won
  and the `BUILD_INFO` line. Any exception means: log it and fall back to the frozen legacy variable, never crash.
- *Launch string:* `. <build>/run_env.sh && srun <build>/jdftx_gpu`. The environment applies to that launch only.
- *Adopt table* `<prod root>/adopt` (`<frozen legacy exe path> <channel>` per line). This is what moved jobs that were
  already queued: each one had frozen a build path at submission. Without the table, they all stay on the old build.
- *Job scripts:* `export JDFTx_GPU_CHANNEL=<channel>`, with no `export JDFTx_GPU=...` and no JDFTx module loads. To
  retarget one queued job (also for its later resubmissions), run `echo <channel or build dir> > jdftx_build` in its
  directory.

**Rollout that worked:**
1. Mock-test the resolver: override, channel variable, adopt, and a channel-name typo that must fall back.
2. Make patched *copies* next to the originals (`*_rollout_test.py`).
3. Run a debug job: NO3H single point at np=1. Every resolution path must reproduce the reference F to ~1e-12 Ha.
4. Run a 4-rank check with the production SBATCH block (Pt fcc, 4x4x4 k). Expect `Run totals: 4 processes, 128
   threads, 4 GPUs` and equal peak memory on every GPU.
5. Get the user's OK, then rename the copies over the originals.
6. Confirm in the first production job's log that it resolved through the adopt table.

**Where Anvil differs** (cross-check with the Anvil section below):
- The launcher is `mpirun -np N`, not `srun`. `run_env.sh` must therefore also put the matching OpenMPI on `PATH`;
  `LD_LIBRARY_PATH` alone is not enough. The binary's RPATH covers gcc and OpenMPI but not OpenBLAS, FFTW, GSL or CUDA,
  so `run_env.sh` (a frozen copy of the `anvil-env.sh` module loads, or a `PATH`+`LD_LIBRARY_PATH` snapshot) is what
  makes runs independent of the submit shell. Test it from a clean shell:
  `env -i HOME=$HOME bash -l -c '. <build>/run_env.sh && ldd <build>/jdftx_gpu | grep "not found"; which mpirun'`.
- The launch path on Anvil is the same perlStuff `get_exe_cmd` (Anvil session's findings in the Anvil section). Since
  2026-11-05 the resolver is on perlStuff `origin/master` (`f15b6fd`), so **the next `git pull` in Anvil's `~/perlstuff`
  puts it live for every pending Anvil job**. Mock-tested on Perlmutter (2026-11-05, merged code, job Python, Anvil's
  frozen string `mpirun -np 2 .../build-gpu/jdftx_gpu`, `use_srun=False`). With no `~/jdftx_prod` on Anvil, the launch
  string comes out byte-identical to before, so the pull is inert for jobs. The log gains two lines: `JDFTx build: ...
  [$JDFTx_GPU as frozen at submission]` and `No run_env.sh for this build`. To use channels on Anvil, the resolver needs
  two changes:
  (a) **The adopt table can't match Anvil's frozen string.** The parser keeps only lines with exactly two
  whitespace-separated fields and compares build dirs. The full string with `mpirun` has too many fields, and a
  path-only key is compared against the whole frozen string, launcher included. Both fail silently, with no `WARNING`.
  Fix: match on the string's last token (the executable).
  (b) **A resolved channel loses the launcher.** With `use_srun=False`, `get_exe_cmd` returns `. <build>/run_env.sh &&
  <build>/jdftx_gpu`, which has no `mpirun`. Fix: keep the frozen string's prefix before the executable (`mpirun -np 2`),
  or build the launcher from `$SLURM_NTASKS`.
  Both fixes are untested. `calc_helpers.py` and `jdftx_channel.py` are shared with Perlmutter production, so any change
  there must leave the `srun` path and the two-field adopt lines working exactly as they do now. Test both machines'
  cases (the four in "Rollout that worked", plus the Anvil string) before pushing, and ask the user before pulling it
  on either cluster. The user wants these Anvil changes to come from Anvil sessions (2026-11-05). perlStuff's
  `CLAUDE.md` ("Anvil: what the build channels will need") has the same to-do list next to the code, plus the
  check list for changes to `get_exe_cmd`.
- The current build is in the shared group dir `/anvil/projects/x-chm260113/Software/JDFTx/`. Ask the user where the
  production root and the new versioned build dirs should live. Never repoint or rebuild what others run.
  *Answered (Anvil, 2026-11-05):* builds go in that projects Software dir; see "Where Anvil builds live" in the Anvil section.
- Node types: partition `ai` has 21 nodes with 4x H100, besides the 16 nodes with 4x A100 in `gpu`. That comes from
  the Anvil session's `sinfo` (2026-11-05); the user's H100 access is unconfirmed (see the Anvil section). The JDFTx CMake
  ignores `CUDA_ARCH`/`CUDA_CODE` (they have been gone since b202b5f5, also at 52548a96); only `CUDAARCHS` matters. The
  Perlmutter `CUDAARCHS=80` library holds sm_80 SASS plus compute_80 PTX (`cuobjdump`), so an sm_80 build would
  probably JIT-compile on H100 rather than fail. For H100 nodes, build natively with `CUDAARCHS="80;90"` or `90`, and
  give each node type its channel the way Perlmutter pairs `managed` with 40GB.
- A managed-memory flavor (`-D CudaManagedMemory=yes`) needs CUDA 12.x, because the prefetch code does not compile
  with 13.x; Anvil's `cuda/12.6.1` qualifies. Use it only for systems that don't fit in GPU memory. On Perlmutter A100s
  it was ~1.55x slower than the plain build on 220-228-atom slabs and ~5x slower on a small molecule.
- Untested idea for mixed GPU types: an `auto` channel that picks the flavor from the GPU the job lands on
  (`nvidia-smi --query-gpu=name,memory.total --format=csv,noheader` at resolve time). Jobs could then accept either
  node type and draw from the larger pool. Perlmutter doesn't do this; its node constraint and channel are set by hand.

**Using the allocation well: the method transfers, the numbers don't** (Perlmutter, 2026-11-01 to 11-03):
- Optimize work per calendar day, (ElecMinimize iterations per job) / (queue wait + wall time), not speed per
  iteration. 80GB + `default` was 1.55x faster per iteration, but 80GB 3-5 h jobs waited a median ~214 h, against
  ~24 h on 40GB (all users' jobs that started 2026-11-02/03). Estimated throughput: 40GB + `managed` 5 h ~450-600
  iterations/day, 80GB 5 h ~120. Data sources: `sacct` Submit/Start/Elapsed (use Elapsed, not End-Start, across DST)
  and iteration counts in the `out` files.
- Read the scheduler limits before choosing a strategy: `scontrol show config | grep -i "prior\|sched"` and
  `sacctmgr show qos <qos> format=name,maxwall,maxjobspu,maxsubmitpu,maxjobsaccruepu`. Perlmutter `gpu_regular` has
  `MaxJobsAccruePU=2`: only the user's 2 oldest pending jobs age, so jobs start in pairs, in submit order. A deeper queue
  then adds nothing, and a resubmitted segment goes to the back. `PriorityWeightFairShare=0` there.
- Each restart costs ~10-20 min (re-convergence plus the lost partial ionic step). Fewer, longer jobs should raise
  throughput as long as the wait grows more slowly than the wall time (longer than 5 h untested on Perlmutter).
- You can move your own pending jobs to another node type without losing queue position (NERSC):
  `scontrol update JobId=<id> Features="gpu&a100&hbm40g"`. Write the matching `jdftx_build` file first, because moving
  from 80GB to 40GB with `default` OOMs. On Anvil, test whether users may change `Partition=` or `Features=` this way.
- Compare charge rates per node type. On Perlmutter, hbm40g and hbm80g both bill 128 per node-hour. Anvil's GPU
  charging was not in the public docs on 2026-11-05; check the Anvil user guide or `mybalance`-style tools there.

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
  *Note (Perlmutter session, 2026-11-05, from the source, not run on Anvil):* `jdftx/CMakeLists.txt` at `52548a96` no
  longer reads `CUDA_ARCH` or `CUDA_CODE` (they were removed in b202b5f5), so `CUDAARCHS=80` alone sets the targets.
  Check `CMakeCache.txt` for `CMAKE_CUDA_ARCHITECTURES` and cmake.log for "Manually-specified variables were not used".
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

**Anvil state and job-launch path (Anvil session, 2026-11-05; read from files and Slurm, nothing run):**
- Builds in the group dir: `build-gpu/` (2026-09-21, the `52548a96` build above) and `build-gpu-cur/` (2026-10-22). The
  shared `jdftx-git` is now at upstream `e305ab5d`, which contains `544d0e78`. Both builds have `CudaAwareMPI=OFF`,
  `CMAKE_CUDA_ARCHITECTURES=80` and profiling on. **Both have an empty `CMAKE_BUILD_TYPE`**, and `flags.make` shows
  `CXX_FLAGS = -Wall -Wno-unused-result -ftemplate-depth-512 -DGPU_ENABLED` with no `-O`. Host code is therefore built at
  gcc's default -O0, which confirms the mbook2019 pitfall on Anvil. Add `-D CMAKE_BUILD_TYPE=Release` to the next Anvil build.
- Production `$JDFTx_GPU` comes from `~/.bashrc_export` (sourced by `~/.bashrc`):
  `"mpirun -np 2 /anvil/projects/x-chm260113/Software/JDFTx/build-gpu/jdftx_gpu"`, which is the old build and includes the
  launcher. `sbatch` freezes it at submission.
- Job scripts (`anvil_submit.sh` in each calc dir: partition `gpu`, `--ntasks=2 --gpus-per-node=2`, 1 h) source
  `anvil-env.sh` and then run `python /home/x-brich/perlstuff/opt4.py`, **both read at job start**. `opt4.py` sets
  `use_srun=False` when the work dir contains `anvil`, and `get_exe_cmd` returns `$JDFTx_GPU` as it is. A resolver in
  `~/perlstuff/helpers/calc_helpers.py` would therefore reach jobs that are already queued. On Anvil the adopt key is the
  whole frozen string, launcher included. It also means **every `git pull` in `~/perlstuff` changes the code that all
  pending jobs will run**. The 2026-11-05 fast-forward to `0a64e90` was checked: `apply_local_protections` does nothing on
  Anvil (`SLURM_JOB_ID` is set), and `oid["gpu"]`/`oid["ddec6"]` have defaults. Later on 2026-11-05, with the user's
  OK, it was fast-forwarded again to `1fb5533` (resolver included). Before the pull, `jdftx_channel.py` from
  `origin/master` was mock-run on Anvil with the real `$JDFTx_GPU` and a pending job's submit dir. There is no
  `~/jdftx_prod`, no `JDFTX_PROD_ROOT`/`JDFTx_GPU_CHANNEL`, and no `jdftx_build` file in any pending job's dir, so the
  launch string came out byte-identical. The resolver is now live on Anvil but changes nothing until a prod root exists.
- **Where Anvil builds live (user, 2026-11-05):** all JDFTx builds go in the group projects Software directory,
  `/anvil/projects/x-chm260113/Software/JDFTx/`, as new versioned dirs next to `build-gpu`/`build-gpu-cur` (never
  rebuild into or repoint an existing one). The prod root (channels, adopt, HISTORY) will also live there, as
  `.../Software/JDFTx/prod` set through `JDFTX_PROD_ROOT`. The user agreed on 2026-11-05 but asked that it is not created
  yet. **User rule (2026-11-05):** the user is the only Anvil user in the group for now, but others will join. The channels
  will switch every group member's jobs, so production must never be switched to an experimental build before it has
  been tested and verified. The group dir's `README.md` states this, defines "verified" (3 checks), and requires each
  build's `STATUS` file to say `status=verified` with a test record before a production channel points at it. A future
  `jdftx-promote` on Anvil must enforce that (refuse without `status=verified`).
- **Versioned builds on Anvil** (Anvil session, 2026-11-05): driver `.../Software/JDFTx/build-versioned.sh <full SHA> [tag]`.
  It makes a shallow source checkout `src-<tag>-<short>/` from the fork, refuses existing dirs, and runs the
  `run-cmake-gpu.sh` options plus `CMAKE_BUILD_TYPE=Release` (it drops the unused `CUDA_ARCH`/`CUDA_CODE`). It writes
  `run_env.sh`, `BUILD_INFO` and `STATUS=unverified`. `run_env.sh` prepends only the module `/apps/...` PATH entries to
  the caller's PATH, because a full snapshot leaked the builder's conda env and `~/bin`, and it sets the full module
  `LD_LIBRARY_PATH` and `JDFTX_GPU_SYNC_MPI=1`.
  First build: **`build-gpu-eff-4fe71c77`** (fork `efficiency_testing` 4fe71c77, the same code as tip 6a74e4e2 and as
  Perlmutter's channels; Release; `-O3 -DNDEBUG` in `flags.make`; sm_80; CudaAwareMPI OFF; default flavor; 563 MB plus a
  210 MB source tree; `make -j16` on login00 ~10 min). Login-node checks passed: from `env -i` shells, with and without a
  login profile, `ldd` shows 0 not found and its own `libjdftx_gpu.so`, and `mpirun` is OpenMPI 4.1.6 matching
  `libmpi.so.40`. The CPU `jdftx -n` dry run succeeded (git hash 4fe71c7). **Status: unverified, not used by any job.**
  Next: a `gpu-debug` reference comparison against `build-gpu` (README check 2), then a 40GB production-like run
  (check 3). The Anvil GPUs are A100-SXM4-40GB (from `out` files), and production jobs run 2 ranks × 64 threads.
- **Anvil scheduler (Anvil session, 2026-11-05; Slurm 26.05.4, read from config/sacct/sprio, nothing changed):**
  - `gpu` and `gpu-debug` share nodes g[000-004,006-015] and have the same PriorityTier/JobFactor, so debug jobs get no
    priority boost. QOS `part-gpu`: MaxWall 2 days, MaxTRESPU gres/gpu=12. QOS `part-gpu-debug`: MaxWall 30 min,
    MaxJobsPU 1, MaxSubmitPU 2.
  - Priority is age only: Age weight 20000 (PriorityMaxAge 7 days), FairShare 0, Partition 750, QOS 500, JobSize 500.
    **Only 4 of a user's pending jobs accrue age at a time.** That was observed through `AccrueTime`: the user's 4 oldest
    had it set and the other 62 showed N/A with age 0. The limit doesn't appear in the `sacctmgr` qos/assoc output. As on
    Perlmutter, queueing deeper than a few jobs gains no priority.
    *Correction (same session): the limit is per ACCOUNT, not per user.* `scontrol show assoc_mgr flags=qos qos=gpu`
    shows `MaxJobsAccruePA=4` on job QOS `gpu`, which `gpu-debug` jobs also run under. So all chm260113-gpu members
    (x-sgerits is also in the account) share 4 accrue slots, and pending `gpu` jobs keep a debug job from accruing.
  - **2026-11-05 ~13:50: all 66 of the user's pending jobs were put on user hold** (`scontrol hold`, with the user's OK)
    so the gpu-debug tests of `build-gpu-eff-4fe71c77` have free accrue slots and backfill. Verified: `JobHeldUser`,
    priority 0, AccrueTime reset, account accrue count `4(0)`, submit times unchanged. **Release them when testing ends:**
    `scontrol release <ids>`, with the ids in `queue_logs/hold_2026-11-05.log` (also in `pending_2026-11-05.tsv`).
  - `bf_max_job_assoc=50`: backfill tests at most 50 of a user's jobs, in priority order. With 50 or more of the user's
    own older jobs pending, a new debug job is never tested. Evidence: the earlier debug jobs, with no older jobs of
    the user pending, waited 0.8-26 h. Debug job 21591351 waited 13.7 days with 158 older jobs pending; that count fell
    below 50 on 11-03 21:20 and the job started 36 h later.
  - User's `gpu` jobs started since 2026-10-22: 124, median wait 274 h, median run 1.00 h. Since 10-01, 184 of 276 that
    ended hit TIMEOUT (1 h or 2 h limits). Users can only lower a pending job's time limit, so a longer walltime needs a
    resubmission.
  - Full pending-queue snapshot: `/anvil/scratch/x-brich/claude_scratch/queue_logs/pending_2026-11-05.tsv`.
- Node types (`sinfo`, 2026-11-05): partitions `gpu`/`gpu-debug` have 16 nodes with 4x A100. Partition `ai` has 21 nodes
  with 4x H100 (feature `H100`, QoS `part-ai`). The user's only association is `chm260113-gpu` on partition `gpu`, so H100
  access is unconfirmed.

### mbook2019 (Ben's home CPU server: 2019 MacBook Pro 16", Ubuntu 26.04)

Verified 2026-11-05 by the Claude session on that machine; its system log is `~/SERVER_NOTES.md` there. Intel i9-9880H (8 cores /
16 threads, AVX2), 64 GB, no GPU. CPU package power capped at 45 W (sustained ~3.0 GHz, no thermal throttling). Single-node Slurm
25.11 (partition `main`, `SLURM_HINT=nomultithread` in `~/.bashrc`, cgroup core/memory limits), OpenMPI 5.0.10 with external PMIx:
`srun` launches MPI directly (`MpiDefault=pmix`).
- **Builds** (identical CMake caches, compared option by option): `~/software/jdftx` = upstream `b3c54e1d`;
  `~/software/jdftx-fork` = fork `efficiency_testing` `4839c10f` (CLAUDE.md protocol set up in `src/`). **Since 2026-11-05 the fork is the
  default** (`jdftx` on PATH, job template); the upstream build is kept for A/B (`JDFTX=~/software/jdftx/build/jdftx sbatch ...`). Ubuntu packages: gcc 15.2,
  CMake 4.2.3, OpenBLAS 0.3.32 (pthreads), FFTW 3.3.10, GSL 2.8, LibXC 5.2.3. No CompileNative (more nodes may be added).
  ```bash
  cmake -D CMAKE_BUILD_TYPE=Release -D EnableLibXC=yes -D CBLAS_LIBRARY=/usr/lib/x86_64-linux-gnu/libopenblas.so ../src/jdftx
  srun -n 1 -c 16 --hint=multithread make -j16      # 72 s
  ```
- **Two CMake pitfalls seen here:** (1) no build type means no `-O` (see the note under the Perlmutter recipe); (2) FindCBLAS picks
  `libgslcblas` (GSL's reference CBLAS) ahead of OpenBLAS. `CBLAS_LIBRARY` fixes it; checked with `LD_BIND_NOW=1 LD_DEBUG=bindings`
  that `cblas_zgemm` then binds to libopenblas.
- Binaries link with `RUNPATH` (Ubuntu's default), so an `LD_LIBRARY_PATH` would override the build's own `libjdftx.so`. None is set here.
- **Running:** always `OPENBLAS_NUM_THREADS=1`. JDFTx threads itself; with OpenBLAS threads = `-c`, one process ran 22 threads on 8 cores
  and was 25-36% slower. **MPI over states beats threads here** (bulk Si, PBE, GBRV, 8x8x8 k = 29 states, electronic-SCF, all 8 cores:
  1x8 23.1 s, 2x4 6.8 s, 4x2 4.5 s, 8x1 3.6 s). Job template `~/slurm-tests/jdftx-job-template.sh` (default `-n 4 -c 2`; `-n 1 -c 8`
  for one Gamma-only system or molecule). **Several molecules: run them side by side** (2 x `-c 4` gave 1.24x the throughput of
  1 x `-c 8` on solvated glycinate, 2026-11-06). Memory sets the limit: ~12.6 GB per solvated 20 Å box, ~3.8 GB gas phase.
- **Tests:** ctest 11/11 (47 checks) for both builds, threads-only (`-n 1 -c 8`) and MPI (`JDFTX_LAUNCH="srun -n %d -c 2"`).
- **`efficiency_testing` A/B** (4839c10f vs b3c54e1d, back to back, alternating order, 2 repetitions, otherwise idle machine; inputs and
  scripts in `~/software/jdftx-fork/bench/`; details in `efficiency_notes.md`, "mbook2019 CPU A/B"):

  | case | layout n×c | upstream | fork | change |
  |---|---|---|---|---|
  | bulk Si, GBRV, 29 states, electronic-minimize (13 iters both) | 4×2 | 0.151 s/iter | 0.113 s/iter | −25% per iter (−14% wall) |
  | Pt(111) 2x2x3, 12 atoms, 144 bands, 216 projectors, 4 states, Fermi 0.01 (22 iters both) | 4×2 | 2.91 s/iter | 2.52 s/iter | −13% per iter (−9% wall) |
  | D3, 180-atom Pt(111) 6x6x5 slab, lj-override 105.83544 + ionic-minimize 1 (13 D3 calls) | 1×8 | 36.0 s/call | 0.038 s/call | ~950x |

  Energies agree to 8e-15 Ha (Si), 5e-11 Ha (Pt) and 3e-13 Ha (D3), with the same iteration counts. These systems are small and the
  hardware differs, so the per-iteration gains are not directly comparable with the 27% for the 218-atom slab on Perlmutter GPUs.
- **Molecule thread scaling** (2026-11-06, solvated glycinate in a 20 Å box, details in `efficiency_notes.md`, "mbook2019: thread
  scaling on a solvated molecule"): 8 threads = 2.6x one thread. Causes: the 45 W cap (one core ~4.2 GHz, eight ~2.9 GHz),
  17-19% thermal throttling at 8 threads, and memory bandwidth (solvated cycles are ~11x slower than the M5 Pro, gas ~4x).
  Open lead: ~160 M minor page faults per solvated SCF cycle, 15-20% of CPU time in the kernel.
- **Python:** conda env `dft` (Miniforge, conda-forge, no compilers/MPI in the env) with `benrich37/pymatgen-core@misc_fixes` installed
  editable from `~/software/pymatgen-core`. 2026-11-05: its `tests/io/jdftx` gave 390 passed and 1 tolerance-only failure (exact-0 compare at
  `atol=0`, see `~/envs/pymatgen-jdftx-issues.md` there).
- **pyjdftx** (`jdftx/opt/pyjdftx`, built 2026-11-06 against the fork build `4839c10f`; script `~/envs/build-pyjdftx.sh` there, 7 s).
  No env changes beyond pure-Python build tools: conda-forge `scikit-build-core` + `pybind11`; `mpi4py` pip-built from source against the
  system OpenMPI (`pip install --no-binary mpi4py mpi4py`; conda-forge's mpi4py would pull conda's own OpenMPI). Then
  `pip install --no-build-isolation -C cmake.define.CMAKE_CXX_COMPILER=/usr/bin/mpicxx -C cmake.define.JDFTX_LIBRARY=<build>/libjdftx.so
  -C cmake.define.CMAKE_INSTALL_RPATH='$ORIGIN/../../..;<build>' <src>/jdftx/opt/pyjdftx` from a shell with the env active and the
  system cmake. The `$ORIGIN/../../..` entry (the env's lib) matters: without it, `import pyjdftx` before numpy loads the OS libstdc++
  (older than conda's), so later conda C++ modules get the older one. Re-run the script after every rebuild of the JDFTx build.
  Results (Si, 29 states, electronic-SCF): `pyjdftx.run(input)` = CLI to 2e-14 Ha at 1x8 and 4x2; ASE calculator = CLI to 2e-13 Ha.
  - *pyjdftx log bug (all machines; found here 2026-11-06, source not changed):* `initialize(comm, commAll, log, False)` opens `log`
    with "w" on **every** rank, so a rank that gets there late truncates what rank 0 already wrote (NUL bytes at the top of the log;
    the energies are unaffected). The CLI is safe because `-o` opens with "a". Workaround: pass `os.devnull` on ranks other than 0.
    Fix: open the file only on the head rank in `_pyjdftx.cpp::initialize`. **User (2026-11-06): known bug, not to be fixed for now.**
    Use the workaround in every multi-rank pyjdftx script, and don't patch `_pyjdftx.cpp`.
  - *Unexplained, here only:* at 1x8 (one process, 8 threads) pyjdftx ran 18 s against the CLI's 23 s on the same input (3 reps each,
    0.71 vs 0.96 s per SCF cycle, 5x fewer minor page faults). At 4x2 the two are the same. Ruled out: library binding (LD_DEBUG), the
    env's libstdc++, the MPI thread level and the glibc mmap threshold. On this 2-atom cell 1x2 (9.4 s) is also faster than 1x8, so
    threading overhead dominates tiny systems. *2026-11-06:* not seen on molecules: on glycinate (gas and solvated) pyjdftx = CLI
    to <1% with the same page-fault counts, so it only matters for tiny cells.

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

### Testing and benchmarking the van der Waals module (default method)

To test or time `VanDerWaals*` code, skip electronic minimization with these built-in options. Don't write drivers or
patch `run.cpp` for this.
- **Repeated calls: `lj-override` with D3.** With `lj-override`, setup sets `skipWfnsInit`, and `IonInfo::update` and
  `ionicEnergyAndGrad` return right after `pairPotentialsAndGrad`. The minimizers then skip all electronic work, and
  Ewald runs with zero charges. For D3/D3BJ the only other change is the cutoff (`rCut = ljOverride` instead of 200 bohr
  in `VanDerWaalsD3::energyAndGrad`), so pass `lj-override 105.83544` (= 200 bohr) to reproduce production D3 (checked: energy matches a normal run to 6e-12 Ha).
  Add `elec-cutoff 5` (grids and basis are still set up), then `ionic-minimize nIterations N` (forces),
  `lattice-minimize` (adds stress) or `ionic-dynamics` (a fixed number of steps). Each ionic step calls vdW twice:
  energy in `update`, then forces in `ionicEnergyAndGrad`. **Not for D2:** `lj-override` replaces D2 with a
  Lennard-Jones potential.
- **One call: dry run `jdftx -n`.** `Everything::setup()` runs `iInfo.update()`, which evaluates vdW once (energy and
  forces, no stress). This works for D2 too.
- Timing: build with `-D EnableProfiling=yes` and read `PROFILER: VanDerWaalsD3::energyAndGrad` (or `...D2...`), which
  counts only the vdW routine.

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
Merged with upstream through b3c54e1d (2026-10-26). Code changes: projection reuse (below; ctest 11/11 on CPU and
GPU, upstream merge checks in `efficiency_notes.md`) and the DFT-D3 pair-sum rewrite `a37d055f` (2026-10-28; see
"DFT-D3 speedup" below). `efficiency_notes.md` also has a 2026-10-28 survey of further candidates (`nAugmentGrad_gpu`
launch overhead, skipping the XC MPI split on GPU, batched FFTs, Cholesky step irreversibility); none are in the code yet.
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
- **CPU A/B on mbook2019 (Intel i9, 2026-11-05):** −25% / −13% per ElecMinimize iteration (bulk Si / 12-atom Pt slab), D3 ~950x faster
  per call, identical energies and iteration counts. Table in "mbook2019" under Cluster builds; details in `efficiency_notes.md`.

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
4. Production builds: change 1 on top of master, without change 2. Done 2026-10-26 for `74dd6931`
   (`~/jdftx_compile_test/ben/build-gpu-{,managed-}eff-74dd6931`, validated against the `54315e95` builds); builds of the
   branch tip `4fe71c77` (adds the D3 speedup) are still pending there, see below.
5. The 9 `wfns` copies were deleted on 2026-10-26 (with the user's OK). The run directories keep `in`, `out`, the small state files and `run_info.txt` (502 MB). A rerun needs `wfns` copied again from the production directory.

#### DFT-D3 speedup `a37d055f` (GPU results, Perlmutter, 2026-10-29)

`a37d055f` (branch tip `4fe71c77`) rewrites the pair sums in `VanDerWaalsD3.cpp` only (header unchanged): C6
interpolation weights computed once per atom and contracted on plain arrays instead of three tiny matrix products per
ordered atom pair (on GPU each was a cuBLAS call with transfers), neighbor cells listed once per call, threads split the
first atom of each pair, integer powers instead of `std::pow`. The heteronuclear C6 indexing bug in `getPairParams` is
kept on purpose (marked TODO) so results are unchanged. Tested with the lj-override method ("Testing and benchmarking the
van der Waals module" above) in the production layout (4 ranks × 32 threads, 1 A100 each, `JDFTX_GPU_SYNC_MPI=1`,
`JDFTX_MEMPOOL_SIZE=36000`): old = the `74dd6931` production builds, new = `4fe71c77` built with the same recipes (CMake
cache identical), run back to back on the same node. Workspace `/pscratch/sd/b/beri9208/claude_scratch/eff_d3/`
(README with results, `build.sh`, `inputs/gen_vdw_input.py`, `job_vdw.sh`, `job_e2e.sh`, `compare.py`, outputs in
`runs/`, Slurm logs in `logs/`); debug jobs 60503563, 60503565, 60503566. Full write-up with setup details:
"DFT-D3 speedup `a37d055f`: Perlmutter GPU results" at the end of `efficiency_notes.md` (`efficiency_testing`).

| case (D3BJ unless noted) | managed build, 40GB node: old → new | default build, 80GB node: old → new |
|---|---|---|
| 224-atom slab, `ionic-minimize nIterations 1` (13 D3 calls) | 41.3 → 0.26 s/call | 19.3 → 0.20 s/call |
| 255-atom slab, `ionic-minimize nIterations 1` (5 calls) | 53.4 → 0.29 s | 24.9 → 0.25 s (3.1 s with 1 thread/rank) |
| 224-atom slab, `lattice-minimize nIterations 1` (stress) | | 19.6 → 0.23 s |
| 224-atom slab, D3 zero damping | | 19.1 → 0.22 s |
| NO3H molecule, D3 | | 0.028 → 0.014 s |

Old and new agree to ≤5e-15 Ha (energy, EvdW) and ≤1e-14 Ha/bohr (forces); stress and coordination numbers are identical.

End-to-end on the production configuration (40GB node, managed builds, production env and `in` of the 224-atom
`Ni_B/CH2COO/p07V/Mn_Mn-O_O-bi` run, 10 ElecMinimize iterations with `energyDiffThreshold 0` from its state files),
current `$JDFTx_GPU` `build-gpu-managed-54315e95` vs managed `4fe71c77`: 30.0 → 19.2 s/iter (−36%), D3 42.3 → 0.27 s per
call, MEMUSAGE total 42.96 → 36.26 GB per GPU, same trajectory (identical step events; G within 5e-8 Ha after 10
iterations, final forces within 4e-7 Ha/bohr). In the current production runs (224–228-atom slabs, PBE+U+D3BJ+CANDLE,
`ionic-minimize`, 5 h jobs) ElecMinimize is ~83% of wall time at ~35 s/iter and D3 (14–20 calls × ~40 s) ~4%, so
switching the primary build from `managed-54315e95` to a `4fe71c77` build should cut wall time per ionic step by about
a third (~1.5× throughput); about 5 points of that are the D3 change, the rest projection reuse.

Status 2026-10-29: `~/.jdft.bashrc` still points `$JDFTx_GPU` at `build-gpu-managed-54315e95` (only the user can edit
that file). The `4fe71c77` builds exist only on scratch (`eff_d3/build-gpu-{,managed-}4fe71c77`); production copies in
`~/jdftx_compile_test/ben` are not made yet. Open question: whether the non-managed build now fits these slabs on 40GB
nodes (peak 36.3 GB with projection reuse; the default build did 11.9 s/iter on the 218-atom slab on 80GB in Test A, vs
19.2 s/iter managed on 40GB here, so it could be much faster if the cache does not push it over 40 GB). Untested.

Status 2026-11-03 (Perlmutter): resolved. The `4fe71c77` builds were copied to `~/jdftx_compile_test/ben/build-gpu-{,managed-}eff-4fe71c77`
(each with `run_env.sh` and `BUILD_INFO`) and promoted to the `managed` and `default` channels on 2026-10-29 14:40; production jobs
pick them up through perlStuff's `jdftx_channel.py` ("Perlmutter production builds and channels" above), so `~/.jdft.bashrc` no
longer matters for new `psubmit.sh` files. The open question is answered: the non-managed build does NOT fit the 219-266-atom
adsorbate slabs on 40GB nodes (peak ~45 GiB/GPU, OOM); managed `4fe71c77` runs them there at 18.3 s/iter (220-228 atoms) and
~25 s/iter (242-266 atoms) in restart runs, ~1.55x slower than `default` on 80GB, and 24 such 5 h jobs ran on 2026-11-01/02
without an OOM. Queue wait decides: 40GB 3-5 h jobs ~24 h vs 80GB ~214 h median on 2026-11-02/03.

# HPC and Slurm Workflows

This file is the source of truth for preparing, reviewing, and running long
simulations or Slurm workflows in FLOWPanel.jl.

## When and Where to Run

Login nodes are for editing, inspection, and light checks only. Jobs expected
to take more than roughly 20–30 laptop minutes belong on HPC. Request only the
resources and hardware features the case needs, and configure account, QOS,
constraints, and modules for the user's allocation instead of assuming
project-specific values.

Repository Slurm scripts are submitted from the top level of the
`FLOWPanel.jl/` checkout. Keep project, example, data, and output paths relative
to that working directory, for example `--project=.`, `examples/case.jl`, and
`data/case/`. Do not derive or export a `REPO_ROOT` from `BASH_SOURCE`, and do
not `cd` inside a batch script unless the user explicitly requests a launcher
that works from arbitrary directories.

## Where to Submit (partition and QOS choice)

Do not guess a partition. Refresh the live picture first, with the thresholds
set to the job's real ask:

```bash
ssh orc -fN                                   # only if the ControlMaster socket is cold

SKILL=~/.claude/skills/slurm-availability/slurm_availability.py
python3 $SKILL --cpus 16 --mem-gb 96          # p018 hover chains (--mem=96G)
python3 $SKILL --cpus 64 --mem-gb 500         # p021 benchmarks (zen3 512 GiB nodes only)
python3 $SKILL --gpus 1 --print               # GPU work, full table
```

This is the system-wide `slurm-availability` skill
(`~/.claude/skills/slurm-availability/`), not a repo script — invoke the skill
by name and it will explain itself. It rewrites
`~/.claude/slurm/orc_availability.csv` in place, one row per partition
("cluster"), from a read-only `sinfo` / `scontrol show partition` /
`sacctmgr show assoc` / `squeue` probe piped over the existing ssh socket.
Nothing is deployed cluster-side.

Match `--cpus`/`--mem-gb` to what the job actually requests or the `nodes_fit`
column is meaningless: several ORC partitions have 28-core nodes, so a 32-core
threshold silently reports them as unusable. Probe two or three points across
the plausible ask to find the knee — the p021 500 G profile is restricted to
zen3 512 GiB nodes, while a 96 G ask opens most of the cluster.

Before committing to a partition, confirm with `--eta`:

```bash
python3 $SKILL --cpus 28 --mem-gb 96 --time 3-00:00:00 --eta
python3 $SKILL --qos-map        # which partitions each QOS reaches
```

`--eta` runs `sbatch --test-only` (allocates nothing) and reports estimated
start times plus rejections. It catches gates the node counts cannot: `knlg`
looks like 11 free nodes but needs `-C knl`; `qos=test` caps walltime at 60
minutes; GPU partitions reject CPU-only jobs. Treat a `REJECTED` row as
disqualifying.

A QOS is an entitlement, not a queue: `-p m9,m12 --qos=normal` queues across
both and takes the earliest free node. But Slurm **silently drops** partitions
whose `AllowQos` excludes the QOS, so multi-partition only works within one QOS
family.

The decisive column is `access`:

- `access=normal` — a non-preemptible QOS is available; `qos_normal` lists the
  alternatives, pass **one** as `--qos=`. Jobs here are not preempted.
- `access=preempt` — reachable only through `standby`/`gstandby`. The job runs
  on someone else's private hardware and **will** be requeued when the owner
  wants it (`preempt_mode=REQUEUE`). Only submit here with `--requeue` plus
  row-level resume, as in `benchmark/slurm/p2_table.sh`.
- `access=none` — not submittable; kept in the CSV for inventory only.

For GPU work always pass `--gpus N` (or `--gpu-report`): it adds a by-model
breakdown across all reachable partitions, including fully-busy ones that the
shortlist filters out. On ORC the H200s are two separate pools reached by
different QOS — 32 in `m13h` (`--qos=gpu`) and 8 in `eng` (`--qos=eng`), 40
cards total — so a per-partition view alone under-reports the model.

GPU columns (`gpus_total`, `gpus_alloc`, `gpus_idle`, `gpu_types`,
`gpu_free_types`, `nodes_gpu`, `nodes_gpu_free`) come from per-node
`Gres`/`GresUsed`. `gpus_total` is capacity including down/maint nodes;
`gpus_idle` and `gpu_free_types` count only GPUs on *schedulable* nodes. GPU
QOS follows the same split — `gpu` and `eng` are normal, `gstandby` is
preemptible.

Queue depth (`pend_jobs`, `run_jobs`) and the user's own footprint (`my_pend`,
`my_run`) are in the same row: a partition with idle nodes and a 500-deep
pending queue is idle for a reason, usually a reservation or a node feature the
queue is waiting on.

If the probe comes back empty or truncated the script prints `PROBE-FAIL`, exits
non-zero, and leaves the CSV untouched. An empty answer is never "nothing
available" — re-open the ssh socket and retry.

## Standard Single-Node Julia Launcher

Use the following pattern for an ordinary single-process, threaded Julia job:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=flowpanel-case
#SBATCH --nodes=1
#SBATCH --ntasks=<N>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<HH:MM:SS>
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Submit from the top level of the FLOWPanel.jl checkout.
set -euo pipefail

THREADS=<N>

export JULIA_NUM_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export BLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

echo "FLOWPanel case"
echo "  repo:    $(pwd)"
echo "  threads: $THREADS"

julia --project=. -t "$THREADS" examples/case.jl
```

Replace both occurrences of `<N>` with the CPU count appropriate for the case;
it may be lower or higher than the counts used by existing launchers. The
`--ntasks=<N>` convention is standard, but no particular value of `<N>` is a
repository-wide default. Slurm launchers must declare explicit memory, time,
output-log, and error-log directives. Use `set -euo pipefail` so a failed
workflow stage stops the job. Slurm opens log paths before the script runs, so
any requested log directory must already exist when the user submits the job.

If a run does not require GPU and should not exceed 1 hour of walltime, use
the test queue by appending a `SBATCH qos=test` line. This will dodge the waiting queue.

If a run should use Julia or Python, include the appropriate module with `module load julia python` etc. Note that Julia 1.12 is used on the HPC by default, but 1.11 can be requested during the `module load julia` command by specifying the julia version.

Set the single `THREADS=<N>` variable explicitly to the same CPU count requested
by `#SBATCH --ntasks=<N>`, and export it consistently through
`JULIA_NUM_THREADS`, `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`BLAS_NUM_THREADS`, and `MKL_NUM_THREADS`. Pass the same value to Julia with
`-t "$THREADS"` and print it near the start of the log for auditability.

The batch script already executes inside its allocation, so invoke Julia
directly. Avoid `srun` for ordinary single-process threaded Julia jobs. It is
appropriate only for a genuine multi-process launch step or when required by
the site, in which case its task and CPU request must match the allocation.

## Inputs and Checkpoints

Before submission, confirm that every input asset used by the script is present
on the cluster checkout. An untracked local CSV or mesh will not appear in a
fresh clone; copy it to the matching repository path or pass its cluster path
explicitly through the script's environment. Do not silently substitute a
different geometry or polar input merely to make a job start.

For multi-stage jobs, retain completed, validated stage artifacts in the run
directory and let the driver checkpoint-skip that stage on resubmission. The
rotor axial comparison treats both tagged CCBlade sectional CSVs, its polar
CSV, and its validation report as the complete XFOIL/CCBlade checkpoint. Set
`FORCE_CCBLADE=1` only when those polars must be deliberately rebuilt.

## Simulation Outputs

When launching a simulation from an example, especially a diagnostic or
reproduction run for an active investigation, leave VTK output enabled by
default. The I/O cost is usually small relative to the value of retaining
ParaView-ready state. Do not set `SAVE_VTK=false` unless the user asks for a
no-output run.

To avoid filling the disk across repeated iterations, ParaView output from prior
runs must be moved off `/home`. **Since Ryan's ruling of 2026-08-27 that is done
by ARCHIVING, not deleting**: a finished run is tarred whole into
`/nobackup/archive/usr/rander39/FLOWPanel_runs/`, and only then are its ParaView
files removed from `/home` — **except the newest 5 restartable steps, which stay
put** (Ryan, 2026-08-27). An archived run is therefore still warm-startable
directly off `/home`: continuing it needs no tarball unpacking. The `.csv` logs,
`*_case_metadata.toml`, `*.metadata.toml`, `*.pvd` and `monitors/` stay on
`/home` exactly as before — they are small, they are the scientific output, and
`/home` is backed up daily where the archive is not. Deletion without an archive copy is now the fallback
path, not the default, and applies only to runs that are still writing.

Don't let FLOWPanel take up more than 400G at any time. The cap is measured
across ALL of `/home/rander39` (raised from a 200G FLOWPanel-only budget by
Ryan, 2026-08-24; lowered 500G -> 400G, with every G tier dropped 100 G, by Ryan
2026-08-25), not just the FLOWPanel checkout.

When the user explicitly wants to preserve a previous run for comparison, ask
before overwriting and offer to move the old directory aside.

### The archive tier (Ryan, 2026-08-27)

Facts from [BYU RC storage](https://rc.byu.edu/wiki/?id=Storage), measured on the
cluster 2026-08-27. Get these wrong and the reclaim path fails in ways that are
not obvious from an error message:

| | `/home/rander39` | `/nobackup/archive/usr/rander39` |
|---|---|---|
| quota | 2 TiB / 2 M files | 20 TiB / **1 M files** |
| backed up | **yes, daily** | **no** |
| snapshots | daily | no |
| auto-delete | no | no |
| visible from a compute node | yes | **no — login nodes only** |
| baseline 2026-08-27 | 383 G used | **61.64 G / 92,998 files** |

Consequences that drive the design:

- **A Slurm job cannot read or write the archive.** Archiving happens from a
  login node, never from inside a batch script. A "permission denied" or
  "no such directory" on `/nobackup/archive` almost always means the wrong node.
- **The binding constraint is the file count, not the bytes** — 1 M files against
  20 TiB. Per-timestep VTK is exactly the wrong shape for that, which is why
  runs are archived as **one tarball each** rather than copied as trees. At
  ~900 k files free, tree-copying a few campaigns would exhaust the inode budget
  while using a rounding error of the space.
- **The archive is not backed up.** Once a run's VTK is removed from `/home`, the
  tarball is the only copy. This is why the archiver verifies by reading the
  tarball back before deleting anything, and why the CSV/monitor residue is
  deliberately left on the backed-up filesystem.
- **The archive is slow** and the wiki warns of migration to cheaper media.
  It is cold storage: fine for a tarball you may never open, wrong for anything
  a job or an analysis touches repeatedly.

### Retention ladder — the LIVE-run fallback (Ryan, 2026-08-24; tiers lowered 100 G 2026-08-25; demoted to fallback 2026-08-27)

**Archive everything finished first.** The ladder below applies only to runs
still being written, which cannot be tarred consistently and so cannot be
archived. Reach for it only when archiving every finished run still leaves the
cap breached. 36 remains the absolute floor.

| Home usage | Retention | Sweep style |
|---|---|---|
| below 300 G | keep **288** steps | normal whole-tree sweep once past the 150 G trigger |
| above 300 G | keep **72** steps | escalation: oldest runs first, one at a time, **stop as soon as usage drops below 300 G** |
| still breaching 400 G at 72 | keep **36** steps | last resort only; requires the sweeper's `--last-resort` flag |

The escalation tier is deliberately incremental — runs that were not needed to
get back under 300 G keep their full 288 steps. It is driven by a **dry-run
work-list**: one `--keep 72` dry run (~60 s) reports `delete=<n>files/<mb>MB`
per run, and only runs with non-zero yield are swept, in oldest-first order.
The between-run stop check uses `df -BM /home/rander39` (instant, and valid
because home is its own dedicated VAST mount — agrees with `du` to 0.8%), not
`du`, which takes ~14 min per call. One authoritative `du` per cycle still
provides the reported headline. Measured 2026-08-24: 8 of 189 runs held all
9953 MB of reclaimable VTK; the naive re-measure-everything loop would have
cost ~27 h of inode walking to free nothing. **36 steps is the absolute
floor**, and is a genuine last resort: it applies only if sweeping everything at
72 still leaves the cap breached (Ryan, 2026-08-24). Dropping a run to 36 cannot
be undone, so it is gated behind an explicit flag rather than reachable by a
typo.

`scripts/p018_vtk_sweeper.sh` enforces the ladder: a `--keep` value that is
non-integer, below 36, or below 72 without `--last-resort` exits 4 and deletes
nothing. It also refuses an `--only` that matches no run directory (exit 5) and
holds an flock so two sweeps cannot race on the same tree (exit 6). Verified
against a synthetic fixture on the cluster, 2026-08-24 — before that guard, a
`--keep 0` or an empty `--keep` silently deleted 100% of a run's VTK.

### Enforcing the 400 G cap

#### Scope: every checkout, not just this one

There are several FLOWPanel clones under `/home/rander39`, each with its own
`data/` tree, and **all of them are in scope** (Ryan, 2026-08-27). Measured
2026-08-27, the main checkout is only ~84 G of the ~398 G total; the bulk sits in
`FLOWPanel-018-gpu-h200` (95.7 G), `FLOWPanel-018-gpu-gh200` (80.2 G),
`FLOWPanel-052-h200` (64.6 G) and `FLOWPanel-052` (23.2 G). A reclaim pass that
only looks at `projects/FLOWPanel.jl` is looking in the wrong place.

`run_archiver.sh --root DIR` (repeatable) or `--all-checkouts` covers them.
`--all-checkouts` treats a directory as a checkout only if it has **both**
`Project.toml` and `data/` — `Project.toml` alone would drag in unrelated Julia
projects, `data/` alone would match scratch directories. One protect list governs
all of them, since run names are campaign-wide rather than per-clone.

#### Archive layout: which clone did this come from?

```
/nobackup/archive/usr/rander39/FLOWPanel_runs/
├── INDEX.tsv                                  <- one append-only row per archived run
├── projects_FLOWPanel.jl/<run>.tar.zst
├── FLOWPanel-018-gpu-h200/<run>.tar.zst
└── FLOWPanel-052/<run>.tar.zst
```

The slug is the checkout path relative to `$HOME` with `/` → `_`. **Basename
alone is not enough** — `projects/FLOWPanel.jl` and `flowpanel-021/FLOWPanel.jl`
would collide, and so would identically-named runs inside them.

Provenance is recorded three ways, so losing any one of them is survivable:

1. **the slug directory** the tarball sits in;
2. **`ARCHIVE_SOURCE.txt` inside the tarball** — absolute source path, checkout
   slug, git rev, host, UTC timestamp. A tarball that gets moved or renamed still
   says where it came from. It is staged in a temp dir and tarred with a second
   `-C`, never written into the run directory: touching `/home` would bump the
   run's mtime and make it read LIVE on the next pass;
3. **`INDEX.tsv` at the archive root** — `archived_utc, checkout_slug, run,
   tarball, src_bytes, tar_bytes, files, kept_steps, git_rev, source_checkout`.
   "Which clone did this come from?" is one `grep`.

`--restore` and `--resume-delete` act on exactly one checkout and **refuse to
guess**: with several roots in play they exit 2 and ask for `--root`. Restoring
into the wrong clone would be silent and wrong.

#### The two scripts

Two scripts, with a strict order of precedence:

| script | scope | action |
|---|---|---|
| `scripts/run_archiver.sh` | **finished** runs | tar whole run → archive → verify → remove ParaView files from `/home` outside the newest **5** restartable steps |
| `scripts/p018_vtk_sweeper.sh` | **live** runs | fallback only: delete VTK outside the newest **288** restartable steps |

Note the two keep windows are deliberately very different, and both are correct.
The sweeper's 288 is a *safety margin against irreversible loss* — what it
deletes is gone. The archiver's 5 is a *convenience window* — everything it
removes is already in a verified tarball, so the only question is how much stays
instantly warm-startable without a `--restore`.

`run_archiver.sh` runs first and does the bulk of the reclaim. It is dry-run by
default (`--apply` acts) and honours the same protect list. Locking (since
2026-08-28) lives on the SHARED filesystem under `$LOCKROOT` (default
`~/.cache/flowpanel/locks`, mkdir-based, so it holds across login nodes):
archivers publish reader markers and may run **in parallel, one worker per
checkout** (a `checkout.<slug>` lock makes two workers on one checkout
impossible); the sweeper takes `sweeper.excl` exclusively, so a sweep and an
archive can never touch `data/` at once. No lock auto-expires — a stale lock
after a crash is removed by a human, and the error message names the exact
`rm -rf`. Its exit codes:
`2` bad usage — including an `--only` run name that matches directories in more
than one checkout, which must be scoped with a single `--root` (run names are
not unique across clones), `3` no protect list, `5` unmatched
`--only`/`--restore`, `6` lock conflict (the sweeper holds the shared lock,
or a checkout was skipped because another archiver holds its lock), `7` archive
unreachable or unwritable
(**you are probably on a compute node**), `8` a tarball failed verification —
*those* runs deleted nothing — or a `--restore` that left members missing or
truncated on disk (`RESTORE-INCOMPLETE`; the tarball stays intact), `9`
`ARCHIVED-STALE` runs found and a human decision is owed.

**Launching workers — detached by default (Ryan, 2026-08-28).** Long archiver
passes are launched detached on a login node so they survive a dropped ssh
connection or a sleeping laptop; run in the foreground only for quick dry-runs.
One worker per checkout may run in parallel — disjoint `--root`s, never two
workers on one checkout (the lock enforces it, but launch disjoint on purpose):

```bash
ssh orc 'bash -lc "cd /home/rander39/projects/FLOWPanel.jl && \
  log=/home/rander39/archiver_worker_<tag>_\$(date +%Y%m%d).log; \
  nohup nice -n 19 env ZSTD_THREADS=2 ./scripts/run_archiver.sh \
      --root <checkout> --apply > \$log 2>&1 < /dev/null & \
  echo LAUNCHED pid=\$! log=\$log"'
```

Logs go on `/home` (shared, backed up), never `/tmp` (node-local). Poll the log
for the per-run labels and summary sentinels; the launcher MUST collect every
worker's final `STALE_COUNT=`, `VERIFY_FAIL_COUNT=` and exit status and refuse
to report success otherwise — a `VERIFY-FAIL` lost in one of N scrollbacks
means "do not trust this pass" went unread. `ZSTD_THREADS=2` with 3–4 workers
keeps total load inside login-node etiquette; back off if the arbiter
complains.

Its contract, in order:

1. **Three-way classification, not a binary** (Ryan, 2026-08-27). The two
   liveness signals are reported separately rather than OR-ed, because measured
   on the cluster they disagree in *both* directions:

   | class | condition | action |
   |---|---|---|
   | `LIVE` | a queue job name matches | never touched |
   | `RECENT` | no queue job, but touched < 24 h ago | **ask Ryan** — not archived |
   | `ARCHIVE` | no queue job, quiet ≥ 24 h | archived automatically |

   `RECENT` exists because "no Slurm job + recently written" is genuinely
   ambiguous: it is usually a run that finished an hour ago and is ready to
   archive, but it can also be a live job the name match missed. The archiver
   reports it with its quiet age and VTK size and moves on; Ryan approves by
   name with `--include-recent --only RUN[,RUN...]`. `--include-recent` without
   `--only` is refused outright — it would silently become "archive everything
   not in the queue", which is exactly the judgement the label escalates.

   **Why the queue is advisory and mtime is the real signal.** Measured
   2026-08-27: job `fp-018-csarc_n5_nt144_l2p4` normalises to a string that is a
   substring of **seven** different run directories, six of which were not being
   written at all — while job `fp-018gpu-n2_nt72_l3p0` normalises to
   `018gpu_n2_nt72_l3p0`, which is **not** a substring of the directory it was
   actively writing, `p018_csarc_n2_nt72_l3p0`. Do not "fix" this by trusting the
   name match harder.

   **The backstop for all of it**: the newest mtime under the run is snapshotted
   before the tar, re-checked after verification, and re-checked AGAIN
   immediately before the delete (`CHANGED-BEFORE-DELETE`), together with a
   fresh protect-list read (`PROTECTED-LATE`) and a fresh squeue poll
   (`LIVE-LATE`). If anything moved, the tarball is quarantined as
   `.changed.<host.pid.start>` and **nothing is deleted** — regardless of how
   the run was classified. That is what makes an approved `RECENT` run safe
   even if the approval was wrong.
2. **Tar the whole run directory** — VTK plus CSVs, TOMLs, `.pvd` and
   `monitors/` — streamed straight to the archive under a `.partial` name that
   is renamed only on success, so an interrupted transfer can never be mistaken
   for a complete archive.
3. **Verify by reading the tarball back**: the member path set and the total
   content byte count must both match the source. On failure the tarball is
   renamed `.corrupt.<host.pid.start>` (so the next pass re-archives rather than trusting
   it) and **nothing is deleted**.
4. **Only then** remove ParaView files from `/home` — but only those **outside
   the newest `--keep N` restartable steps (default 5)**. A step counts as
   restartable only if all four paths the warmstart loader reads exist, the same
   definition the sweeper uses, and *all* files at a kept step are kept (`.vtm`
   indices together with their pieces) so the retained state is both
   warm-startable and ParaView-openable. The CSV/TOML/pvd/monitors residue stays,
   and an `ARCHIVED.txt` breadcrumb records the tarball, the kept steps, and the
   restore command.

   A run with **no** complete restartable step has all its VTK removed — there is
   no restart to preserve, and unlike the sweeper's refusal to guess, here the
   verified tarball still holds every file. Runs with fewer steps than the window
   keep everything and free nothing; that is expected, not a failure.
   `--keep 0` strips all VTK and prints a warning: legal, but it means continuing
   the run requires a `--restore` first.

`--restore RUN` is the inverse and is tested to be byte-identical; it refuses to
overwrite existing files.

**`ARCHIVED-STALE` means stop and ask Ryan** (Ryan, 2026-08-27). It fires when a
tarball exists *and* ParaView files survive **outside the keep window** — i.e.
the archive landed but the delete did not complete. VTK inside the window is
expected and reads as `ARCHIVED-ALREADY`, which is the normal end state. Do not re-tar, do not delete, do not
"tidy it up": report it and wait. Once Ryan authorises it,
`--resume-delete RUN` finishes the job, but only after confirming every
remaining ParaView file is present in the tarball **byte for byte**; any
mismatch exits 8 and deletes nothing. A deliberate `--restore` produces the same
on-disk state, so restore leaves a `RESTORED.txt` marker that annotates the
alert — it never suppresses it.

Behaviour is pinned by `scripts/tests/run_archiver_test.sh` (synthetic fixture,
no cluster needed), which covers classification, the residue that must survive,
the restore round-trip, stale detection, corrupt-tarball quarantine, and every
exit code. **Run it after any edit to the archiver and before mirroring.**

**Do not hand-roll VTK deletions or archive moves.** A step is only
restartable if all four paths the warmstart loader reads exist, and `.vtm`
indices must be kept together with their pieces — a sweep that left index stubs
without pieces killed job 13036477 on 2026-08-04. The sweeper encodes that rule;
ad-hoc `rm` does not.

**Both scripts** must be mirrored to the cluster checkout whenever they change
(`/home/rander39/projects/FLOWPanel.jl/scripts/`); verify with `md5sum` on both
sides afterwards. `run_archiver.sh` is the one that matters most here — it is
the only tool that writes outside the checkout, and a stale cluster copy would
delete under rules the local copy no longer holds. The agent definitions in `.claude/agents/` are local-only —
there is no `.claude/` on the cluster. **Known divergence:** the cluster's copy
of this file (`agent_policies/HPC.md`) is an older, line-patched revision and is
deliberately not synced (Ryan, 2026-08-24); nothing on the cluster reads it, so
treat the local copy as authoritative.

**Launch the `hpc-storage` subagent** rather than doing this inline whenever:

- a job that writes VTK is submitted, or one is running unattended;
- home-directory usage crosses roughly 150 G — well below the cap, because
  growth scales with the number of live VTK writers (~5.4 G/h with 8 writers
  measured 2026-08-24; ~24 G/h under a ~20-writer campaign) and sweeping late is
  a race;
- **a campaign or arm has finished and its runs are now eligible to archive**
  (Ryan, 2026-08-27) — this is now the main source of reclaim, and it does not
  wait for a disk alarm;
- a new sweep arm is about to start; or
- the cluster reports any disk-space or quota error.

Pass it the run names in flight if you know them. It performs one cycle
and returns a before/after report plus a ledger line for the active campaign's
`ledger.md`; it has no write tools, so the main session appends that line.

For continuous coverage during a campaign, the *main session* installs a
`Monitor` that probes `du -sm /home/rander39/projects/FLOWPanel.jl` every ~10
minutes and wakes on ≥100 G, then launches `hpc-storage` for the cycle. Harden
the probe: require the specific line ending in `projects/FLOWPanel.jl` and treat
its absence as a probe failure rather than a small number, `tr -cd '0-9'` before
any arithmetic (the login banner emits ANSI codes on the first stdout line), and
allow two consecutive failures before alarming.

The VTK protect list
(`BRAINSTORM/018_dji9443_hover_convergence_campaign/vtk_protect_list.txt`, with
a matching cluster copy that is the one that takes effect) is Ryan's file:
agents read it and never write it. If a run needs protecting, ask him.

For BYU agent policy, read `/apps/instructions_for_ai_agents/BYU_ORC_AGENTS.md` from the login node.

See BYU's [Slurm guidance](https://rc.byu.edu/wiki/?id=Slurm) and
[script generator](https://rc.byu.edu/documentation/slurm/script-generator) if you are having trouble.

## Unified project location (2026-08-31)

All orc jobs launch from `~/projects_unified` — NOT from the legacy per-task
silo copies (`~/<pkg>-<era>`, `~/fm052env-*`: deprecated, kept only until
their in-flight jobs finish, then removed by Ryan).

**Code.** `~/projects_unified/{FastMultipole,FLOWVPM.jl,FLOWPanel.jl}` — real
git clones, branch `unified-052`. Pin a run by COMMITTING (local-only) before
launching; use `git worktree add` for concurrent experiments needing
different code states. Deploy local changes by scp/rsync + commit on orc
(never launch from an uncommitted tree you can't reproduce).

**Any job script must set three things, all arch-keyed by `$(uname -m)`:**

1. **Working dir**: `cd ~/projects_unified/FLOWPanel.jl` (or the relevant
   package) so `logs/slurm/`, `logs/chain/`, and `data/<case>/` conventions
   hold in-tree.
2. **Project**: `--project=$HOME/projects_unified/envs/$(uname -m)` —
   `x86_64` or `aarch64`; the Manifests dev-point at the unified trees.
3. **Julia + depot, by arch**:
   - x86 (m13h / eng / CPU partitions): `module load cuda
     julia/1.11.7-6bmogfl`; default depot `~/.julia`.
   - ARM (mgh GH200): no x86 module tree — use
     `$HOME/julia/julia-1.11.7/bin/julia` directly, `export
     JULIA_DEPOT_PATH=$HOME/fm052depot-gh200` (holds the ARM CUDA artifacts —
     do not delete or rename), and no `cuda` module: CUDA comes from CUDA.jl
     artifacts plus the node driver.

**Submission inputs by pool:**

| Pool | sbatch flags |
|---|---|
| H200 (default) | `-p m13h --gres=gpu:h200:1` |
| H200 (alternate) | `-p eng --qos=eng --gres=gpu:h200:1` |
| GH200 (ARM) | `-p mgh --qos=gpu --gres=gpu:gh200:1 -C arm` |

**Existing launchers** (e.g. `FLOWPanel.jl/examples/run_p018_screen_gpu052.slurm.sh
<arch> <case>`) already implement this pattern and take `*_REPO_OVERRIDE` /
`*_PROJECT_OVERRIDE` env vars for special cases — copy their arch-dispatch
header rather than reinventing it. CAVEAT: they still source GPU tuning knobs
from the legacy silo path `~/FLOWVPM-052-<arch>/scripts/fm052_common.sh`;
migrate that bundle into `~/projects_unified` before the silos are deleted.

**Validation**: envs + `FLOWVPMCUDAExt` verified on all three GPU pools and
x86 CPU (2026-09-01). After any env/depot change, re-run
`~/projects_unified/smoke_unified_{x86,gpu_m13h,arm}.slurm.sh` (adjust
partition/gres as needed).

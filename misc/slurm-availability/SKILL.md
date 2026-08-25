---
name: slurm-availability
description: This skill should be used when Slurm jobs are not starting, or when deciding what resources to request and where to submit. Triggers on "jobs stuck in the queue", "pending forever", "low priority", "which partition should I use", "how many cores should I ask for", "which nodes are free", "cluster availability", "idle nodes", "sinfo", "GPU availability", "are there free GPUs", "which QOS", "normal queue vs preemption", "standby", "when will my job start", and "Slurm capacity". Probes a cluster read-only over ssh, refreshes a per-partition CSV, and reports a candidate shortlist, GPU capacity by model, QOS reachability, and non-allocating start-time estimates.
allowed-tools: Bash(python3 *), Bash(ssh *), Read
---

# Slurm availability

Size the request, then find where it lands. Probe, report, and hand the decision
to whoever asked — do not pick the partition yourself.

`SCRIPT=~/.claude/skills/slurm-availability/slurm_availability.py`

## 1. Establish the job's real ask first

Before probing, read the actual requirements out of the sbatch script or the
user's description: threads (`--ntasks` / `THREADS`), `--mem`, `--gres`, and
`--time`. Establish one thing the CSV can never know: **does the job survive
requeue?** (checkpoint/restart, or `--requeue` plus row-level resume). That
single fact decides whether preemptible partitions are usable at all.

If the ask is unknown, say which thresholds you assumed.

## 2. Run the probe

```bash
python3 $SCRIPT --host orc --cpus <CORES> --mem-gb <GB> [--gpus <N>]
```

Needs a live ssh ControlMaster socket. If ssh hangs or prompts for 2FA, **stop**
— report that there is no live socket and ask the user to run `ssh <host> -fN`.
Never retry into a 2FA loop.

The script is local and self-contained: the probe is a string piped to
`ssh <host> bash -s`. Nothing is installed on the cluster. It runs only
`sinfo`, `scontrol show`, `sacctmgr show`, `squeue`, and — with `--eta` —
`sbatch --test-only`, which allocates nothing.

Flags: `--print` full table · `--top N` longer shortlist · `--gpus N` /
`--gpu-report` GPU breakdown · `--qos-map` QOS reachability · `--eta` start-time
estimates · `--time` walltime for `--eta` · `--host`, `--slurm-bin`, `--out`.

## 3. Read `access` — it decides the submission

Derived mechanically from the user's QOS associations intersected with each
partition's `AllowQos`:

- **`normal`** — a non-preemptible QOS. `qos_normal` lists the alternatives;
  pass **one** as `--qos=`. The job will not be preempted.
- **`preempt`** — reachable only through `standby`/`gstandby`, on someone
  else's hardware, `PreemptMode=REQUEUE`. The job **will** be requeued. Only
  send work here that step 1 established is requeue-tolerant.
- **`none`** — not submittable; in the CSV for inventory only.

Also check, per candidate: `maxtime` ≥ the job's walltime; `pend_jobs` (idle
nodes alongside a deep queue means a reservation or an unmet node feature, not
free capacity); and `features` if the job needs a specific microarchitecture.

## 4. Sizing the request: iterate the thresholds

**This is the step that turns the tool from a report into a sizing exercise.**
The ask and the availability are coupled — `nodes_fit` is only as meaningful as
the thresholds fed to it, and a threshold above a partition's per-node core
count zeroes that partition out *silently*.

Run the probe across the plausible range and report the tradeoff, rather than a
single verdict:

- **Cores.** `--cpus 32` reported `m9` as `nodes_fit=0`; `--cpus 28` reported
  **226**, because m9 nodes are 28-core. If the job can run on fewer threads,
  dropping the ask can buy the largest idle pool on the cluster.
- **Memory.** `max_free_mem_gb` shows which ceiling excludes whom. A 500 G ask
  restricts to big-memory nodes; 96 G opens most of the cluster.
- **GPU model.** Loosening "must be H200" to "any GPU" moved a test job from a
  six-day wait to minutes.

Phrase the finding as *"at 64 cores you wait; at 28 you start now"* and let the
user judge whether the smaller ask is scientifically acceptable.

## 5. GPU questions: always pass `--gpus N`

Adds **GPU CAPACITY BY MODEL**, aggregating each model across every reachable
partition with free/total and the QOS that reaches it — including fully-busy
partitions, which the shortlist filters out.

This matters because one model is often split across partitions reached by
*different* QOS. On BYU ORC the H200s are 32 in `m13h` (`--qos=gpu`) plus 8 in
`eng` (`--qos=eng`) — two pools, 40 cards. Never answer "how many `<model>` are
there" from the shortlist; use the by-model section or the `gpu_types` /
`gpu_free_types` CSV columns.

`gpus_total` is capacity including down/maint nodes; `gpus_idle` and
`gpu_free_types` count only schedulable nodes. Note that GPU partitions
typically **reject CPU-only jobs** ("without a need for GPUs").

## 6. QOS and partition mechanics — `--qos-map`

A QOS is an **entitlement, not a queue**. One QOS often reaches several
partitions, and `--qos-map` inverts the CSV to show which:

```
--qos=gpu       m13h,m13l,m9g,mgh
--qos=normal    m11-1,m11-2,m12,m14,m8,m8g,m9,m9g,mgh
--qos=eng       eng
```

- `-p a,b,c --qos=<one>` queues across all of them and takes the earliest free
  node. This is the cheapest latency win available.
- **Slurm silently drops partitions whose `AllowQos` excludes your QOS.** A job
  submitted `-p m13h,eng --qos=gpu` runs only on `m13h`; `eng` is ignored with
  no warning. Multi-partition works *within* one QOS family only.
- Pools reached by **disjoint** QOS cannot be queued together in one job — a
  job carries exactly one QOS. Submit two jobs and cancel the loser.
- Sites may block a QOS from direct request (ORC blocks `gstandby`: "request
  `standby` instead").

## 7. `--eta` — the definitive check before committing

```bash
python3 $SCRIPT --cpus 28 --mem-gb 96 --time 1-00:00:00 --eta
```

Runs `sbatch --test-only` for the top candidates individually and combined per
QOS family. It **creates no job**. Idle-node counts cannot see backfill,
fairshare, or preemption; this can — and it exposes gates the CSV has no way to
represent:

```
m9 (--qos=normal)                   2026-08-25T13:01  on m9
COMBINED m9,m12,m9g (--qos=normal)  2026-08-25T13:01  on m9
knlg (--qos=standby)   REJECTED: Cannot access partition knlg unless 'knl'
                                 node feature is requested with -C/--constraint
COMBINED ... (--qos=test)  REJECTED: timelimit 1440 min exceeds maximum of 60
```

`knlg` led the preempt shortlist with 11 idle nodes and an empty queue, and is
unusable without `-C knl`. **Treat a `REJECTED` row as disqualifying** and pass
the reason on — the QOS walltime caps, required constraints, and GPU-only
policies it reveals are exactly what makes a job pend forever.

Use `--eta` whenever the question is "when will this start" or "why is my job
not running", and before recommending any partition.

## 8. Report, do not decide

Relay the shortlist **and the CSV path** so the caller can read the full table.
Then hand the choice over, calling out: the preempt caveat on any
`access=preempt` candidate; any partition idle but deeply queued; any
walltime/`maxtime` conflict; and every `--eta` rejection.

Do not choose the partition. Required node features and whether the job
survives requeue live in the caller's context, not in the CSV.

## 9. Probe-fail discipline

`PROBE-FAIL` with a non-zero exit means the CSV was **left untouched** — a stale
file is never presented as current. Report the failure and its cause (no socket,
unresolved host, no Slurm at `--slurm-bin`). An empty or failed probe is never
"nothing is available".

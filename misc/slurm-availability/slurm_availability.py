#!/usr/bin/env python3
"""Refresh a per-partition Slurm availability CSV in place.

Answers "why won't my job start, and where should it go instead?" for any
Slurm cluster reachable over ssh.

Entirely local and self-contained: the only thing it depends on remotely is a
live ssh ControlMaster socket. The probe below is held as a string, piped to
the remote shell's stdin (`ssh <host> bash -s`), and never written to the
cluster -- nothing to deploy, nothing to keep in sync, nothing to go stale.

Writes ~/.claude/slurm/<host>_availability.csv with one row per Slurm
partition ("cluster"): node states, schedulable CPUs / memory / GPUs, queue
pressure, and -- the column that actually decides where to submit -- whether
this user reaches the partition through a normal QOS or only through
preemptible standby.

    python3 slurm_availability.py --host orc --cpus 16 --mem-gb 96

Set --cpus/--mem-gb/--gpus to the job's real ask, or the `nodes_fit` column is
meaningless: a 32-core threshold silently hides every 28-core partition.

Probe-fail discipline: if the probe comes back empty or truncated the CSV is
left untouched and the exit status is non-zero. An empty answer is never
reported as "nothing available".

Read-only against the cluster. Python 3 stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# Output lives outside any repo so the tool is usable from any session.
# Keyed by ssh host so several clusters can coexist.
OUT_DIR = Path.home() / ".claude" / "slurm"


def default_out(host: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", host)
    return OUT_DIR / f"{safe}_availability.csv"


# ORC's Slurm install. `ssh host bash -s` is a non-login shell, so the site's
# profile is not sourced and the binaries are usually off PATH.
DEFAULT_SLURM_BIN = "/apps/slurm/latest/bin"

# Read-only Slurm inventory probe, executed on the cluster via `ssh <host> bash
# -s` with this text on stdin. {slurm_bin} is substituted before sending.
# Output is sentinel-delimited (###BEGIN/###END) because login banners and MOTD
# land in the same stream. The parser reads only what sits between sentinels.
# Nothing here submits, cancels, or deletes.
#
# `field:|` in the sinfo -O format uses '|' as the field suffix (verified on
# slurm 24.05.8); Features is last because it is long and comma-laden.
PROBE_SH = r"""
set -uo pipefail

export PATH="{slurm_bin}:{slurm_bin_sbin}:$PATH"

# Distinguish "no Slurm here" from "probe produced nothing", which otherwise
# look identical downstream.
if ! command -v sinfo >/dev/null 2>&1; then
  echo "###BEGIN NOSLURM"
  echo "sinfo not found; PATH={slurm_bin} plus default"
  echo "###END NOSLURM"
  exit 0
fi

emit() {
  local name="$1"; shift
  echo "###BEGIN $name"
  "$@" 2>/dev/null || true
  echo "###END $name"
}

echo "###BEGIN WHOAMI"
echo "${USER:-$(id -un)}"
echo "###END WHOAMI"

emit NODES sinfo -N -h -O \
  'NodeList:|,Partition:|,StateLong:|,CPUsState:|,Memory:|,AllocMem:|,Gres:|,GresUsed:|,Features:|'

# One line per partition, key=value -- far easier to parse than the block form.
emit PARTITIONS scontrol show partition -o

# The QOS set this user may actually request (may be several assoc rows).
emit USERQOS sacctmgr -nP show assoc user="${USER:-$(id -un)}" format=Account,Partition,QOS

# Queue pressure. A pending job's Partition may be a comma-separated list.
emit QUEUE squeue -h -a -O 'Partition:|,StateCompact:|,UserName:|,NumNodes:|'
"""

# QOS names that only ever buy preemptible time.
PREEMPT_QOS = {"standby", "gstandby"}

# Node-state flag suffixes that mean "not schedulable right now".
UNAVAILABLE_FLAGS = set("*~#!%$@^")
PLANNED_FLAG = "-"          # reserved by the backfill scheduler for a queued job

FIELDS = [
    "cluster", "access", "qos_normal", "qos_preempt", "preempt_mode",
    "maxtime", "prio_tier",
    "nodes_total", "nodes_idle", "nodes_mixed", "nodes_alloc", "nodes_down",
    "nodes_planned",
    "cpus_idle", "cpus_total",
    "mem_free_gb", "max_free_mem_gb",
    "gpus_total", "gpus_alloc", "gpus_idle", "gpu_types", "gpu_free_types",
    "nodes_gpu", "nodes_gpu_free",
    "nodes_fit", "fit_cpus", "fit_mem_gb", "fit_gpus",
    "pend_jobs", "run_jobs", "my_pend", "my_run",
    "features", "access_note", "updated_utc",
]


# ---------------------------------------------------------------- probe ----

def run_probe(host: str, timeout: int, slurm_bin: str) -> str:
    """Run PROBE_SH on `host` over the existing ssh socket. No files involved."""
    # Plain .replace, never .format: the probe contains shell braces such as
    # ${USER:-$(id -un)} that str.format would try to interpret as fields.
    script = (PROBE_SH
              .replace("{slurm_bin_sbin}", slurm_bin.rstrip("/") + "/../sbin")
              .replace("{slurm_bin}", slurm_bin))
    try:
        proc = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=60",
             host, "bash -s"],
            input=script.encode("utf-8"),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=timeout,
        )
    except FileNotFoundError:
        sys.exit("PROBE-FAIL no ssh binary on PATH")
    except subprocess.TimeoutExpired:
        sys.exit(f"PROBE-FAIL ssh {host} timed out after {timeout}s "
                 f"(is the ControlMaster socket up? `ssh {host} -fN`)")
    if proc.returncode != 0:
        tail = proc.stderr.decode("utf-8", "replace").strip().splitlines()[-3:]
        sys.exit(f"PROBE-FAIL ssh {host} exited {proc.returncode}: "
                 + " / ".join(tail))
    return proc.stdout.decode("utf-8", "replace")


def split_blocks(text: str) -> dict[str, list[str]]:
    """Pull sentinel-delimited blocks out of the banner-polluted stream."""
    blocks: dict[str, list[str]] = {}
    current, buf = None, []
    for raw in text.splitlines():
        line = raw.replace("\x1b[0m", "").rstrip("\n")
        stripped = re.sub(r"\x1b\[[0-9;]*m", "", line).strip()
        if stripped.startswith("###BEGIN "):
            current, buf = stripped[len("###BEGIN "):].strip(), []
        elif stripped.startswith("###END "):
            if current is not None:
                blocks[current] = buf
            current, buf = None, []
        elif current is not None and stripped:
            buf.append(stripped)
    return blocks


# ---------------------------------------------------------------- parse ----

def parse_gres(spec: str) -> dict[str, int]:
    """`gpu:a100:8(IDX:0-2,4-7),gpu:h100:2` -> {'a100': 8, 'h100': 2}."""
    out: dict[str, int] = defaultdict(int)
    if not spec or spec in ("(null)", "N/A"):
        return out
    flat = re.sub(r"\([^)]*\)", "", spec)          # drop (IDX:..)/(S:..) blobs
    for tok in flat.split(","):
        tok = tok.strip()
        if not tok or not tok.startswith("gpu"):
            continue
        parts = tok.split(":")
        try:
            count = int(re.sub(r"[^0-9]", "", parts[-1]) or 0)
        except ValueError:
            continue
        model = parts[1] if len(parts) >= 3 else "gpu"
        out[model] += count
    return out


def split_state(state: str) -> tuple[str, str]:
    """`mixed-` -> ('mixed', '-');  `idle*` -> ('idle', '*')."""
    flags = ""
    while state and state[-1] in UNAVAILABLE_FLAGS | {PLANNED_FLAG}:
        flags = state[-1] + flags
        state = state[:-1]
    return state, flags


def parse_nodes(lines: list[str]) -> list[dict]:
    nodes = []
    for line in lines:
        f = line.split("|")
        if len(f) < 8:
            continue
        name, part, state, cpustate, mem, allocmem, gres, gresused = f[:8]
        feats = f[8] if len(f) > 8 else ""
        base, flags = split_state(state)
        try:
            a, i, o, t = (int(x) for x in cpustate.split("/"))
        except ValueError:
            a = i = o = t = 0
        mem_mb = int(mem) if mem.isdigit() else 0
        alloc_mb = int(allocmem) if allocmem.isdigit() else 0
        g_tot = parse_gres(gres)
        g_use = parse_gres(gresused)
        nodes.append({
            "name": name, "part": part.rstrip("*"),
            "state": base, "flags": flags,
            "down": bool(UNAVAILABLE_FLAGS & set(flags))
                    or base in ("down", "drained", "draining", "fail",
                                "failing", "maint", "reserved", "unknown",
                                "planned", "powering_down", "power_down"),
            "planned": PLANNED_FLAG in flags,
            "cpus_alloc": a, "cpus_idle": i, "cpus_total": t,
            "mem_mb": mem_mb, "free_mb": max(mem_mb - alloc_mb, 0),
            "gpus_total": g_tot, "gpus_used": g_use,
            "features": [x for x in feats.split(",") if x],
        })
    return nodes


def parse_partitions(lines: list[str]) -> dict[str, dict]:
    parts = {}
    for line in lines:
        kv = dict(re.findall(r"(\w+)=(\S*)", line))
        name = kv.get("PartitionName")
        if not name:
            continue
        parts[name] = kv
    return parts


def parse_userqos(lines: list[str]) -> set[str]:
    qos: set[str] = set()
    for line in lines:
        f = line.split("|")
        if f:
            qos |= {q.strip() for q in f[-1].split(",") if q.strip()}
    return qos


def parse_queue(lines: list[str], me: str) -> dict[str, dict[str, int]]:
    q: dict[str, dict[str, int]] = defaultdict(
        lambda: {"pend": 0, "run": 0, "my_pend": 0, "my_run": 0})
    for line in lines:
        f = line.split("|")
        if len(f) < 3:
            continue
        parts, state, user = f[0], f[1].strip(), f[2].strip()
        for p in (x.strip().rstrip("*") for x in parts.split(",")):
            if not p:
                continue
            mine = (user == me)
            if state == "PD":
                q[p]["pend"] += 1
                if mine:
                    q[p]["my_pend"] += 1
            elif state == "R":
                q[p]["run"] += 1
                if mine:
                    q[p]["my_run"] += 1
    return q


# ------------------------------------------------------------ aggregate ----

def classify(pinfo: dict, user_qos: set[str]) -> tuple[str, list[str], list[str], str]:
    allow = {q.strip() for q in pinfo.get("AllowQos", "ALL").split(",") if q.strip()}
    deny = {q.strip() for q in pinfo.get("DenyQos", "").split(",") if q.strip()}
    allowed = set(user_qos) if allow == {"ALL"} else (allow & user_qos)
    allowed -= deny
    non_preempt = sorted(allowed - PREEMPT_QOS)
    preempt = sorted(allowed & PREEMPT_QOS)
    if non_preempt:
        access = "normal"
    elif preempt:
        access = "preempt"
    else:
        access = "none"

    notes = []
    groups = pinfo.get("AllowGroups", "ALL")
    accounts = pinfo.get("AllowAccounts", "ALL")
    if groups != "ALL":
        notes.append(f"groups={groups}")
    if accounts != "ALL":
        notes.append(f"accounts={accounts}")
    if pinfo.get("State", "UP") != "UP":
        notes.append(f"state={pinfo['State']}")
    return access, non_preempt, preempt, ";".join(notes)


def build_rows(blocks, args, stamp) -> list[dict]:
    me = (blocks.get("WHOAMI") or [""])[0].strip()
    nodes = parse_nodes(blocks.get("NODES", []))
    parts = parse_partitions(blocks.get("PARTITIONS", []))
    user_qos = parse_userqos(blocks.get("USERQOS", []))
    queue = parse_queue(blocks.get("QUEUE", []), me)

    if "NOSLURM" in blocks:
        detail = (blocks["NOSLURM"] or [""])[0]
        sys.exit(f"PROBE-FAIL no Slurm on this host: {detail}. "
                 f"Point --slurm-bin at the site's Slurm bin directory.")
    if not nodes or not parts:
        sys.exit("PROBE-FAIL empty NODES or PARTITIONS block "
                 "(slurm binaries missing from PATH, or ssh returned only the banner)")
    if not user_qos:
        sys.exit(f"PROBE-FAIL no QOS associations found for user {me!r}")

    by_part: dict[str, list[dict]] = defaultdict(list)
    for n in nodes:
        by_part[n["part"]].append(n)

    fit_mb = args.mem_gb * 1024
    rows = []
    for pname, pinfo in parts.items():
        pn = by_part.get(pname, [])
        access, qn, qp, note = classify(pinfo, user_qos)

        gpus_tot: dict[str, int] = defaultdict(int)
        gpus_use: dict[str, int] = defaultdict(int)
        gpus_free_by_model: dict[str, int] = defaultdict(int)
        counts = {"idle": 0, "mixed": 0, "alloc": 0, "down": 0}
        planned = fit = gpu_free_nodes = gpu_nodes = gpus_free = 0
        cpus_idle = cpus_total = mem_free = max_free = 0
        feats: set[str] = set()

        for n in pn:
            for k, v in n["gpus_total"].items():
                gpus_tot[k] += v
            for k, v in n["gpus_used"].items():
                gpus_use[k] += v
            n_gpu_tot = sum(n["gpus_total"].values())
            n_gpu_free = max(n_gpu_tot - sum(n["gpus_used"].values()), 0)
            if n_gpu_tot:
                gpu_nodes += 1
            cpus_total += n["cpus_total"]
            feats.update(n["features"])
            if n["down"]:
                counts["down"] += 1
                continue
            if n["planned"]:
                planned += 1
            counts[{"idle": "idle", "mixed": "mixed",
                    "allocated": "alloc"}.get(n["state"], "down")] += 1
            cpus_idle += n["cpus_idle"]
            mem_free += n["free_mb"]
            max_free = max(max_free, n["free_mb"])
            gpus_free += n_gpu_free
            for k, v in n["gpus_total"].items():
                gpus_free_by_model[k] += max(v - n["gpus_used"].get(k, 0), 0)
            if n_gpu_free > 0:
                gpu_free_nodes += 1
            if (n["cpus_idle"] >= args.cpus and n["free_mb"] >= fit_mb
                    and n_gpu_free >= args.gpus):
                fit += 1

        q = queue.get(pname, {"pend": 0, "run": 0, "my_pend": 0, "my_run": 0})
        g_tot_all = sum(gpus_tot.values())
        g_use_all = sum(gpus_use.values())

        rows.append({
            "cluster": pname,
            "access": access,
            "qos_normal": ";".join(qn),
            "qos_preempt": ";".join(qp),
            "preempt_mode": pinfo.get("PreemptMode", ""),
            "maxtime": pinfo.get("MaxTime", ""),
            "prio_tier": pinfo.get("PriorityTier", ""),
            "nodes_total": len(pn),
            "nodes_idle": counts["idle"],
            "nodes_mixed": counts["mixed"],
            "nodes_alloc": counts["alloc"],
            "nodes_down": counts["down"],
            "nodes_planned": planned,
            "cpus_idle": cpus_idle,
            "cpus_total": cpus_total,
            "mem_free_gb": round(mem_free / 1024, 1),
            "max_free_mem_gb": round(max_free / 1024, 1),
            "gpus_total": g_tot_all,
            "gpus_alloc": g_use_all,
            # free GPUs on *schedulable* nodes only -- GPUs on down/maint
            # nodes are capacity (gpus_total), not availability.
            "gpus_idle": gpus_free,
            "gpu_types": ";".join(f"{k}:{v}" for k, v in sorted(gpus_tot.items())),
            # per-model free, schedulable nodes only -- lets a caller ask
            # "where are the H200s" across partitions with different QOS.
            "gpu_free_types": ";".join(f"{k}:{gpus_free_by_model.get(k, 0)}"
                                       for k in sorted(gpus_tot)),
            "nodes_gpu": gpu_nodes,
            "nodes_gpu_free": gpu_free_nodes,
            "nodes_fit": fit,
            "fit_cpus": args.cpus,
            "fit_mem_gb": args.mem_gb,
            "fit_gpus": args.gpus,
            "pend_jobs": q["pend"],
            "run_jobs": q["run"],
            "my_pend": q["my_pend"],
            "my_run": q["my_run"],
            "features": ";".join(sorted(feats)),
            "access_note": note,
            "updated_utc": stamp,
        })

    rank = {"normal": 0, "preempt": 1, "none": 2}
    rows.sort(key=lambda r: (rank[r["access"]], -r["nodes_fit"],
                             -r["cpus_idle"], r["cluster"]))
    return rows


# ----------------------------------------------------------------- write ---

def write_csv(rows: list[dict], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(out.parent), prefix=".{}.".format(out.name))
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, out)                 # atomic: never a half-written CSV
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def print_table(rows: list[dict], limit: int) -> None:
    cols = [("cluster", 12), ("access", 8), ("qos_normal", 22), ("nodes_fit", 9),
            ("nodes_idle", 10), ("nodes_total", 11), ("cpus_idle", 9),
            ("gpus_idle", 9), ("gpus_total", 10), ("gpu_types", 16),
            ("pend_jobs", 9),
            ("maxtime", 12)]
    print(" ".join(h[:w].ljust(w) for h, w in cols))
    print(" ".join("-" * w for _, w in cols))
    shown = [r for r in rows if r["access"] != "none"][:limit]
    for r in shown:
        print(" ".join(str(r[h])[:w].ljust(w) for h, w in cols))


def _kv(spec: str) -> dict[str, int]:
    """'h200:32;l40s:16' -> {'h200': 32, 'l40s': 16}"""
    out = {}
    for tok in (t for t in spec.split(";") if t):
        k, _, v = tok.rpartition(":")
        out[k or "gpu"] = int(v)
    return out


def print_gpu_capacity(rows: list[dict]) -> None:
    """Aggregate GPUs by model ACROSS partitions.

    The same model is often split over partitions reached by different QOS --
    on ORC the H200s are 32 in `m13h` (qos gpu) plus 8 in `eng` (qos eng). A
    per-partition shortlist filtered on nodes_fit>0 hides that split whenever
    the cards happen to be busy, which is exactly when the caller needs to know
    the other pool exists.
    """
    models: dict[str, list[tuple[str, int, int, dict]]] = defaultdict(list)
    for r in rows:
        if r["access"] == "none" or not r["gpus_total"]:
            continue
        free = _kv(r["gpu_free_types"])
        for model, tot in _kv(r["gpu_types"]).items():
            models[model].append((r["cluster"], free.get(model, 0), tot, r))
    if not models:
        return

    print("\nGPU CAPACITY BY MODEL (all reachable partitions, busy ones included):")
    order = sorted(models.items(),
                   key=lambda kv: (-sum(f for _, f, _, _ in kv[1]),
                                   -sum(t for _, _, t, _ in kv[1])))
    for model, entries in order:
        entries.sort(key=lambda e: (-e[1], -e[2]))
        tf, tt = sum(e[1] for e in entries), sum(e[2] for e in entries)
        print(f"  {model:<8} {tf:>3} free / {tt:>3} total")
        for cluster, free, tot, r in entries:
            qos = (r["qos_normal"] or r["qos_preempt"]).replace(";", "|")
            tag = "" if r["access"] == "normal" else "  [PREEMPT]"
            print(f"      {cluster:<10} {free:>3}/{tot:<3} free   "
                  f"qos: {qos}{tag}")


def qos_partition_map(rows: list[dict]) -> dict[str, dict[str, list[str]]]:
    """Invert the CSV: which partitions does each QOS actually reach?

    A QOS is an entitlement, not a queue. Knowing that one QOS spans several
    partitions is what makes the `-p a,b,c --qos=<one>` trick available.
    """
    out: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: {"normal": [], "preempt": []})
    for r in rows:
        # skip access=none and placeholder partitions with no nodes, e.g. "(auto)"
        if r["access"] == "none" or not int(r["nodes_total"]):
            continue
        for q in (x for x in r["qos_normal"].split(";") if x):
            out[q]["normal"].append(r["cluster"])
        for q in (x for x in r["qos_preempt"].split(";") if x):
            out[q]["preempt"].append(r["cluster"])
    return out


def print_qos_map(rows: list[dict]) -> None:
    qmap = qos_partition_map(rows)
    print("\nQOS -> PARTITIONS (one job carries ONE qos; -p a,b,c spans a row):")
    for q, d in sorted(qmap.items(),
                       key=lambda kv: (-len(kv[1]["normal"]), kv[0])):
        for kind in ("normal", "preempt"):
            if d[kind]:
                tag = "" if kind == "normal" else "  [PREEMPT]"
                print(f"  --qos={q:<10} {','.join(sorted(d[kind]))}{tag}")


# gstandby is blocked by site policy ("request standby instead"), and asking for
# a GPU QOS without --gres is rejected outright, so never auto-pick either.
def pick_qos(row: dict, wants_gpu: bool) -> str:
    opts = [q for q in (row["qos_normal"] or row["qos_preempt"]).split(";") if q]
    opts = [q for q in opts if q != "gstandby"] or opts
    if not wants_gpu:
        opts = [q for q in opts if q not in ("gpu",)] or opts
    return opts[0]


# Non-allocating start-time estimate. `sbatch --test-only` validates the request
# and reports when it *would* start; it creates no job. This sees backfill,
# fairshare, and preemption -- none of which idle-node counts can show.
ETA_SH_HEADER = """
set -uo pipefail
export PATH="{slurm_bin}:$PATH"
probe() {
  echo "###ETA $1"
  shift
  sbatch --test-only "$@" 2>&1 | head -3
  echo "###ETAEND"
}
"""


def build_eta_script(cands: list[tuple[str, str, str]], args,
                     slurm_bin: str) -> str:
    """cands: (label, partition-list, qos)."""
    gres = f" --gres=gpu:{args.gpus}" if args.gpus else ""
    body = ETA_SH_HEADER.replace("{slurm_bin}", slurm_bin)
    for label, parts, qos in cands:
        body += (f'probe "{label}" -p {parts} --qos={qos}{gres} '
                 f'-N1 -n{max(args.cpus, 1)} --mem={args.mem_gb}G '
                 f'-t {args.time} --wrap=hostname\n')
    return body


def parse_eta(text: str) -> dict[str, str]:
    out, label, buf = {}, None, []
    for line in text.splitlines():
        line = re.sub(r"\x1b\[[0-9;]*m", "", line).strip()
        if line.startswith("###ETA "):
            label, buf = line[len("###ETA "):], []
        elif line == "###ETAEND" and label is not None:
            blob = " ".join(buf)
            m = re.search(r"to start at (\S+).*?in partition (\S+)", blob)
            if m:
                pre = "  (preempts running work)" if "Preempts:" in blob else ""
                out[label] = f"{m.group(1)}  on {m.group(2)}{pre}"
            else:
                e = re.search(r"error: (.+?)(?:\s+sbatch:|$)", blob)
                reason = (e.group(1) if e else blob) or "no answer"
                reason = re.sub(r"\s+", " ", reason).strip()
                out[label] = "REJECTED: " + reason[:96]
            label = None
        elif label is not None:
            buf.append(line)
    return out


def print_eta(rows: list[dict], args, host: str, slurm_bin: str,
              timeout: int) -> None:
    """Estimate start times for the top candidates, individually and combined."""
    cands: list[tuple[str, str, str]] = []
    seen_group: set[str] = set()
    for access in ("normal", "preempt"):
        pool = [r for r in rows
                if r["access"] == access and r["nodes_fit"] > 0][:args.top]
        for r in pool:
            qos = pick_qos(r, args.gpus > 0)
            cands.append((f"{r['cluster']} (--qos={qos})", r["cluster"], qos))
        # One combined test per QOS family: the multi-partition latency win.
        by_qos: dict[str, list[str]] = defaultdict(list)
        for r in pool:
            for q in (r["qos_normal"] or r["qos_preempt"]).split(";"):
                if q and q != "gstandby" and (args.gpus > 0 or q != "gpu"):
                    by_qos[q].append(r["cluster"])
        for q, parts in by_qos.items():
            if len(parts) > 1 and q not in seen_group:
                seen_group.add(q)
                cands.append((f"COMBINED {','.join(parts)} (--qos={q})",
                              ",".join(parts), q))
    if not cands:
        print("\nESTIMATED START: no candidates fit; nothing to test.")
        return

    script = build_eta_script(cands, args, slurm_bin)
    try:
        proc = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=60",
             host, "bash -s"],
            input=script.encode("utf-8"), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=timeout)
    except subprocess.TimeoutExpired:
        print("\nESTIMATED START: probe timed out; skipped.")
        return
    etas = parse_eta(proc.stdout.decode("utf-8", "replace"))

    gres = f", {args.gpus} gpu" if args.gpus else ""
    print(f"\nESTIMATED START (sbatch --test-only; no job created) for "
          f"{args.cpus}c/{args.mem_gb}G{gres}, walltime {args.time}:")
    for label, _, _ in cands:
        print(f"  {label:<44} {etas.get(label, 'no answer')}")
    print("  NOTE: estimates move as the queue moves, and a COMBINED row is the"
          "\n  real win -- one qos across several partitions takes the earliest"
          "\n  free node. Slurm SILENTLY drops partitions your qos cannot reach.")


def print_shortlist(rows: list[dict], top: int, args) -> None:
    """Mechanical candidate list. No judgment -- the caller decides."""
    def group(access: str) -> list[dict]:
        # rows are already sorted access-rank, then -nodes_fit, then -cpus_idle
        return [r for r in rows if r["access"] == access and r["nodes_fit"] > 0][:top]

    def line(r: dict) -> str:
        gpu = (f"  {r['gpus_idle']} free {r['gpu_types']}"
               if r["gpus_idle"] else "")
        # '|' not ';' -- these are alternatives to pick ONE of, not a list to pass.
        qos = (r["qos_normal"] or r["qos_preempt"]).replace(";", "|")
        return (f"  {r['cluster']:<12} {r['nodes_fit']:>4} fit  "
                f"{r['nodes_idle']:>4} idle  {r['pend_jobs']:>5} pending  "
                f"maxtime {r['maxtime']:<12} qos: {qos}{gpu}")

    fit = (f"fits >={args.cpus}c/{args.mem_gb}G"
           + (f"/{args.gpus}gpu" if args.gpus else ""))

    normal, preempt = group("normal"), group("preempt")
    print(f"\nTOP CANDIDATES (access=normal, {fit}):")
    for r in normal:
        print(line(r))
    if not normal:
        print("  none -- no non-preemptible partition currently fits this job")

    print(f"\nPREEMPT FALLBACK (standby only -- job WILL be requeued):")
    for r in preempt:
        print(line(r))
    if not preempt:
        print("  none")

    print("\nDecision is yours: check walltime, node features, and whether the "
          "job\nsurvives requeue before choosing. Full table:")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", default="orc",
                    help="ssh alias of the cluster login node (default: orc)")
    ap.add_argument("--out", type=Path, default=None,
                    help="CSV path (default: ~/.claude/slurm/<host>_availability.csv)")
    ap.add_argument("--slurm-bin", default=DEFAULT_SLURM_BIN,
                    help=f"remote Slurm bin dir (default: {DEFAULT_SLURM_BIN})")
    # Thresholds should match the job's real ask -- the --mem it requests and
    # the thread count it will use. A threshold above a partition's per-node
    # core count silently reports that partition as unusable.
    ap.add_argument("--cpus", type=int, default=16,
                    help="idle cores a node needs to count in nodes_fit")
    ap.add_argument("--mem-gb", type=int, default=96,
                    help="free GB a node needs to count in nodes_fit")
    ap.add_argument("--gpus", type=int, default=0,
                    help="free GPUs a node needs to count in nodes_fit")
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--print", dest="do_print", action="store_true",
                    help="also print a ranked table to stdout")
    ap.add_argument("--limit", type=int, default=25,
                    help="rows to show with --print")
    ap.add_argument("--eta", action="store_true",
                    help="estimate start times via sbatch --test-only "
                         "(non-allocating; creates no job)")
    ap.add_argument("--time", default="1-00:00:00",
                    help="walltime used for --eta (default: 1-00:00:00)")
    ap.add_argument("--qos-map", action="store_true",
                    help="show which partitions each QOS reaches")
    ap.add_argument("--gpu-report", action="store_true",
                    help="show the by-model GPU breakdown even without --gpus")
    ap.add_argument("--top", type=int, default=5,
                    help="candidates per group in the shortlist (default: 5)")
    args = ap.parse_args()
    out = args.out or default_out(args.host)

    stamp = _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()
    blocks = split_blocks(run_probe(args.host, args.timeout, args.slurm_bin))
    rows = build_rows(blocks, args, stamp)
    write_csv(rows, out)

    reach = sum(1 for r in rows if r["access"] != "none")
    print(f"{args.host}: {len(rows)} partitions ({reach} reachable)  @ {stamp}")
    print_shortlist(rows, args.top, args)
    if args.gpus > 0 or args.gpu_report:
        print_gpu_capacity(rows)
    if args.qos_map:
        print_qos_map(rows)
    if args.eta:
        print_eta(rows, args, args.host, args.slurm_bin, args.timeout)
    if args.gpus > 0 or args.gpu_report or args.qos_map or args.eta:
        print("\nFull table:")
    print(f"  {out}")
    if args.do_print:
        print()
        print_table(rows, args.limit)
    return 0


if __name__ == "__main__":
    sys.exit(main())

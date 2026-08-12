#!/usr/bin/env python3
"""BRAINSTORM/019 single reproducible entry point: sigma-selection procedure.

Implements the pipeline required by
BRAINSTORM/019_sigma_selection_procedure.md ("P1 -- Formalized procedure
(PRE-REGISTERED 2026-08-06)" + "Reproducibility + polish contract"):

  constants     declare/derive all physical constants (with cross-checks)
  gamma_v       recompute measured peak bound circulation Gamma_v from data
  initializers  sigma_eq, sigma_stab, sigma_0 (R units)
  ingest        run registry + provenance + sigma_over_R + config checks
  margins       per-run outcome, time-to-ignition, stability margin M
  iterate       apply the pre-registered update rule sigma_{k+1}=sigma_k*sqrt(M_k/eps)
  gaps          target-grid coverage + machine-readable gap manifest
  figures       emit per-figure CSV data dirs (no .tex authored here)

Default (no subcommand) runs everything in the order above.

Pre-registered choices (fixed 2026-08-06, see item P1 section):
  eps = 0.2; M = max over steps & particles of dt*Z, read directly from a
  `max_dtZ` wake-health column when present, else the level proxy
  M~ = max over downward steps of -Delta ln(min_sigma_ratio), valid only
  while min_sigma_ratio > 0 and before ignition onset; p99.9 of the per-step
  -Delta ln values is reported alongside (commentary only, never in the
  update rule). Method is a labeled attribute of every reported M.

Hard-error contract: any ingested run whose provenance disagrees with the
invariant knobs (RPM=5400, mesh 45_185_ct4 for screens, NT per registry) or
whose two independent sigma derivations disagree by >2% aborts the pipeline
naming the run and key. Missing provenance aborts unless
--allow-missing-provenance is passed (harvest in progress).

Idempotent: byte-identical outputs given identical inputs (sorted keys,
fixed float formats, no timestamps).

Python 3 stdlib only (+ tomllib, py>=3.11).
"""

import argparse
import csv
import math
import os
import re
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = os.path.join(REPO, "data")
FIGDIR = os.path.join(REPO, "BRAINSTORM", "019_sigma_selection_procedure",
                      "figures")

# ---------------------------------------------------------------------------
# constants (single declaration block -- everything else is derived)
# ---------------------------------------------------------------------------

RPM = 5400.0                 # rotor speed, all registry runs (config-checked)
NT_BASE = 36                 # steps per revolution, base timestep
DT_BASE = 60.0 / RPM / NT_BASE   # = 3.0864e-4 s; cross-checked vs metadata
EPS = 0.2                    # pre-registered margin target (P1)
C_EQ = 2.0                   # viscous-fidelity multiplier on sigma_eq (P1)

# 018 anchor for the ambient-strain scale. Provenance: measured ambient
# strain from the B0 production run, sigma_blowup_mechanism.md section 2
# ("At the measured ambient Zbar = 3.2 s^-1 ... sigma* = 2.166 mm").
# ANCHOR ONLY: the item pre-registers that Zbar is recomputed from probe
# data in a later phase; until that probe exists this constant is the
# declared input, labeled as such in every output.
ZBAR_ANCHOR = 3.2            # s^-1

# nu discrepancy carried in the 018 record: prose used 1.5e-5 m^2/s, the
# case metadata records wake_nu = 1.4334e-5 m^2/s. The METADATA value is
# authoritative (it is what the solver ran with); both are printed.
NU_PROSE = 1.5e-5            # m^2/s (018 prose value, NOT used)

GAMMA_V_SANITY = 0.28        # m^2/s, expected order of peak bound circulation
GAMMA_V_RUN = "p018_L1_visc" # Gamma_v source: measured peak bound circulation
                             # from this run's monitor03 (see cmd gamma_v).
                             # A-priori BEM path is a later deliverable (P2).

MESH_SCREEN = "45_185_ct4"   # invariant mesh for all screen-class runs

# Ignition tripwire thresholds for outcome classification (matches the 018
# forensics convention: negative sigma OR runaway velocity).
IGNITE_MAX_U = 1000.0        # m/s
SIGMA_XCHECK_RTOL = 0.02     # tip-arc vs sgm0 cross-check tolerance
KNOWN_PAIR_ATOL = 1e-3       # tip-arc verification tolerance (absolute)
GAP_MATCH_RTOL = 0.07        # |dsigma/sigma| for a run to cover a grid cell

# ---------------------------------------------------------------------------
# run registry
# ---------------------------------------------------------------------------
# expected_* are the registry's own record (used to config-check provenance,
# and as clearly-labeled provisional values while provenance harvest is in
# progress). job_id "unknown" = pending Slurm-banner harvest.

REGISTRY = {
    # --- screens (8-rev, pulse-less startup, tripwire on) ---
    "scr_ufsig_0p05":  dict(cls="screen", viscous=False, nt=36, ov=2.0, pps=7,
                            role="sigma ladder 0.050R", job_id="unknown"),
    "scr_uf_n1_d3p4":  dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=14,
                            role="sigma 0.0299R, N1 d3.4", short="N1d3.4", job_id="unknown"),
    "scr_ufdt_nt36":   dict(cls="screen", viscous=False, nt=36, ov=2.0, pps=12,
                            role="dt study NT36 (sigma 0.0291R)",
                            job_id="unknown"),
    "scr_ufdt_nt72":   dict(cls="screen", viscous=False, nt=72, ov=2.0, pps=6,
                            role="dt study NT72 (sigma 0.0291R, half dt)",
                            job_id="unknown"),
    "scr_ufsig_0p02":  dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=21,
                            role="sigma ladder 0.020R", job_id="unknown"),
    "scr_uf_n2_d6p5":  dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=14,
                            role="sigma 0.0299R, N2 d6.5", short="N2d6.5", job_id="unknown"),
    "scr_uf_n2_d8p5":  dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=14,
                            role="sigma 0.0299R, N2 d8.5", short="N2d8.5", job_id="unknown"),
    "scr_uf_n1_d2p6":  dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=14,
                            role="sigma 0.0299R, N1 d2.6", short="N1d2.6", job_id="unknown"),
    "scr_uf_n1_d5p0":  dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=14,
                            role="sigma 0.0299R, N1 d5.0", short="N1d5.0", job_id="unknown"),
    "scr_b0":          dict(cls="screen", viscous=False, nt=36, ov=2.0, pps=4,
                            role="B0 carrier (sigma 0.0873R)",
                            job_id="unknown"),
    # --- production ---
    "p018_L1_visc":    dict(cls="production", viscous=True, nt=36, ov=2.4,
                            pps=11, role="L1 40-rev viscous certification "
                            "(no wake-health CSV; stability point only)",
                            job_id="unknown"),
    "p018_ufront_n1":  dict(cls="production", viscous=False, nt=36, ov=2.4,
                            pps=14, role="production inviscid 0.0299R",
                            job_id="unknown"),
    "p018_ufront_n1_visc": dict(cls="production", viscous=True, nt=36, ov=2.4,
                            pps=14, role="production viscous 0.0299R "
                            "(known FAILED)", job_id="13058534",
                            known_failed=True),
    "p018_ufront_s035_visc": dict(cls="production", viscous=True, nt=36,
                            ov=2.4, pps=12, role="production viscous 0.0349R hedge "
                            "(CANCELLED by user at step 997/1080; "
                            "monitors valid through ~996)",
                            job_id="13058988"),
    "p018_L2_visc":    dict(cls="production", viscous=True, nt=36, ov=None,
                            pps=None, expected_sigma_over_R=0.0193,
                            role="L2 viscous, forensics-only ingestion "
                            "(ignition point)", job_id="unknown",
                            known_failed=True, forensics=True),
    # --- Campaign E (019 P4 grid, launched 2026-08-06; all D=3.4 per the
    # C1-clearance ruling; screens carry the direct max_dtZ column) ---
    "scr_p019_s015v":  dict(cls="screen", viscous=True, nt=36, ov=2.4, pps=28,
                            role="P4 grid 0.015R viscous (ignited, OOM)",
                            job_id="13060963", known_failed=True),
    "scr_p019_s020v":  dict(cls="screen", viscous=True, nt=36, ov=2.4, pps=21,
                            role="P4 grid 0.020R viscous (ignited, OOM)",
                            job_id="13060964", known_failed=True),
    "scr_p019_s025v":  dict(cls="screen", viscous=True, nt=36, ov=2.0, pps=14,
                            role="P4 grid 0.025R viscous (ignited; Julia "
                            "crash post-ignition)", job_id="13060965",
                            known_failed=True),
    "scr_p019_s030v":  dict(cls="screen", viscous=True, nt=36, ov=2.4, pps=14,
                            role="P4 grid 0.030R viscous (COMPLETED 8 revs)",
                            job_id="13060966"),
    "scr_p019_s015":   dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=28,
                            role="P4 grid 0.015R inviscid (ignited, OOM)",
                            job_id="13061166", known_failed=True),
    "scr_p019_s025":   dict(cls="screen", viscous=False, nt=36, ov=2.0, pps=14,
                            role="P4 grid 0.025R inviscid (ignited, OOM)",
                            job_id="13061167", known_failed=True),
    "scr_p019_s038":   dict(cls="screen", viscous=False, nt=36, ov=2.4, pps=11,
                            role="P4 grid 0.038R inviscid (COMPLETED 8 revs)",
                            job_id="13061168"),
    "scr_p019_s038v":  dict(cls="screen", viscous=True, nt=36, ov=2.4, pps=11,
                            role="P4 grid 0.038R viscous (P4b reference; "
                            "COMPLETED 8 revs)", job_id="13061169"),
    "scr_p019_fid144": dict(cls="screen", viscous=True, nt=144, ov=2.4, pps=5,
                            role="P4b fidelity discriminator NT144 "
                            "(sigma 0.0209R = 1.18 sigma_eq; IGNITED "
                            "~step 910 = 6.4 revs, OOM)",
                            job_id="13061089", known_failed=True),
    "scr_p019_sstab":  dict(cls="screen", viscous=True, nt=36, ov=2.5, pps=14,
                            role="sigma_stab-alone initializer probe "
                            "(sigma 0.0312R; COMPLETED 8 revs)",
                            job_id="13064696"),
}

# tip-arc relation verification pairs: (overlap, p_per_step, NT) -> sigma/R
KNOWN_PAIRS = [
    (2.4, 14, 36, 0.02993),
    (2.0, 7, 36, 0.0499),
    (2.4, 21, 36, 0.01995),
    (2.0, 12, 36, 0.02909),
    (2.4, 28, 36, 0.01496),   # Campaign E s015(v), banner-verified
    (2.0, 14, 36, 0.02494),   # Campaign E s025(v)
    (2.4, 11, 36, 0.03808),   # Campaign E s038(v)
    (2.4, 5, 144, 0.02094),   # fid144 (NT=144)
    (2.5, 14, 36, 0.03117),   # sstab probe
]

# target grid for the P4 regime map (screen class), and the OV/PPS pair to
# submit for any uncovered cell (NT=36).
GRID_SIGMAS = [0.015, 0.020, 0.025, 0.030, 0.038]
GRID_JOBSPEC = {
    0.015: (2.4, 28),   # 0.01496R
    0.020: (2.4, 21),   # 0.01995R
    0.025: (2.0, 14),   # 0.02494R
    0.030: (2.4, 14),   # 0.02993R
    0.038: (2.4, 11),   # 0.03808R
}


def die(msg):
    print(f"HARD ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def tip_arc_sigma_over_R(ov, pps, nt):
    """sigma/R from the shedding geometry: overlap * (tip arc per step) /
    particles per step = ov * (2*pi/NT) / pps."""
    return ov * (2.0 * math.pi / nt) / pps


def fmt(x, spec=".6f"):
    if x is None:
        return "-"
    if isinstance(x, float) and math.isnan(x):
        return "nan"
    return format(x, spec)


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

def load_provenance(run):
    """Return (dict, source) from case_metadata.toml, else banner_config.toml,
    else (None, None)."""
    for tag in ("case_metadata", "banner_config"):
        p = os.path.join(DATA, run, f"{run}_{tag}.toml")
        if os.path.isfile(p):
            with open(p, "r") as fh:
                text = fh.read()
            # Julia writes bare `NaN`/`Inf`; TOML requires lowercase.
            text = re.sub(r"=\s*(-?)NaN\b", r"= \g<1>nan", text)
            text = re.sub(r"=\s*(-?)Inf\b", r"= \g<1>inf", text)
            doc = tomllib.loads(text)
            if tag == "banner_config":
                # banner schema is nested: [provenance] + [config];
                # flatten config to top level, carry job id/state along.
                meta = dict(doc.get("config", {}))
                prov = doc.get("provenance", {})
                if "job_id" in prov:
                    meta["_job_id"] = str(prov["job_id"])
                if "state" in prov:
                    meta["_state"] = str(prov["state"])
            else:
                meta = doc
            return meta, tag
    return None, None


def config_check(run, meta, spec):
    """Invariant-knob verification. Mismatch = hard error (item contract)."""
    if meta.get("RPM") is not None and float(meta["RPM"]) != RPM:
        die(f"config check failed: {run} key RPM = {meta['RPM']} != {RPM}")
    nt_exp = spec["nt"]
    if meta.get("NT") is not None and int(meta["NT"]) != nt_exp:
        die(f"config check failed: {run} key NT = {meta['NT']} != {nt_exp}")
    if spec["cls"] == "screen" and meta.get("mesh_key") is not None \
            and meta["mesh_key"] != MESH_SCREEN:
        die(f"config check failed: {run} key mesh_key = {meta['mesh_key']} "
            f"!= {MESH_SCREEN}")
    dt_exp = 60.0 / RPM / nt_exp
    if meta.get("dt") is not None and \
            abs(float(meta["dt"]) - dt_exp) > 1e-12:
        die(f"config check failed: {run} key dt = {meta['dt']} != {dt_exp}")


def verify_known_pairs():
    for ov, pps, nt, want in KNOWN_PAIRS:
        got = tip_arc_sigma_over_R(ov, pps, nt)
        if abs(got - want) > KNOWN_PAIR_ATOL:
            die(f"tip-arc relation verification failed: OV{ov}/PPS{pps}/"
                f"NT{nt} -> {got:.5f}, expected {want:.5f} (atol 1e-3)")


def build_registry(allow_missing, R_ref):
    """Ingest provenance for every registry run; derive sigma_over_R;
    config-check. Returns (runs dict, missing list)."""
    verify_known_pairs()
    runs, missing = {}, []
    for run in sorted(REGISTRY):
        spec = dict(REGISTRY[run])
        meta, source = load_provenance(run)
        rec = dict(run=run, spec=spec, meta=meta, prov_source=source)
        if meta is not None:
            config_check(run, meta, spec)
            r_ref = float(meta.get("R", R_ref))
            ov = meta.get("overlap")
            pps = meta.get("p_per_step")
            nt = int(meta.get("NT", spec["nt"]))
            sig_arc = (tip_arc_sigma_over_R(float(ov), float(pps), nt)
                       if ov is not None and pps is not None else None)
            # core size in meters: case_metadata's core_spreading_sgm0,
            # or the banner's sigma_m
            sgm0 = meta.get("core_spreading_sgm0", meta.get("sigma_m"))
            sig_sgm0 = (float(sgm0) / r_ref
                        if sgm0 is not None and r_ref else None)
            if sig_arc is not None and sig_sgm0 is not None:
                if abs(sig_arc - sig_sgm0) / sig_sgm0 > SIGMA_XCHECK_RTOL:
                    die(f"sigma cross-check failed for {run}: tip-arc "
                        f"{sig_arc:.5f} vs sgm0/R {sig_sgm0:.5f} (>2%)")
                rec["sigma_over_R"] = sig_sgm0
                rec["sigma_method"] = "sgm0/R (x-checked vs tip-arc)"
            elif sig_sgm0 is not None:
                rec["sigma_over_R"] = sig_sgm0
                rec["sigma_method"] = "sgm0/R"
            elif sig_arc is not None:
                rec["sigma_over_R"] = sig_arc
                rec["sigma_method"] = "tip-arc"
            else:
                die(f"{run}: provenance present but no sigma derivable")
            # banner also carries its own derived sigma_over_R: cross-check
            if meta.get("sigma_over_R") is not None:
                sig_banner = float(meta["sigma_over_R"])
                if abs(rec["sigma_over_R"] - sig_banner) / sig_banner \
                        > SIGMA_XCHECK_RTOL:
                    die(f"sigma cross-check failed for {run}: pipeline "
                        f"{rec['sigma_over_R']:.5f} vs banner sigma_over_R "
                        f"{sig_banner:.5f} (>2%)")
            if "core_spreading_active" in meta and \
                    bool(meta["core_spreading_active"]) != spec["viscous"]:
                die(f"config check failed: {run} key core_spreading_active "
                    f"= {meta['core_spreading_active']} != registry viscous "
                    f"{spec['viscous']}")
            rec["viscous"] = bool(meta.get("core_spreading_active",
                                           spec["viscous"]))
        else:
            missing.append(run)
            # provisional sigma from the registry's own record, labeled
            if spec.get("ov") is not None:
                rec["sigma_over_R"] = tip_arc_sigma_over_R(
                    spec["ov"], spec["pps"], spec["nt"])
                rec["sigma_method"] = "registry OV/PPS (PROVISIONAL)"
            else:
                rec["sigma_over_R"] = spec["expected_sigma_over_R"]
                rec["sigma_method"] = "registry constant (PROVISIONAL)"
            rec["viscous"] = spec["viscous"]
        rec["nt"] = spec["nt"]
        rec["dt"] = 60.0 / RPM / spec["nt"]
        rec["cls"] = spec["cls"]
        rec["role"] = spec["role"]
        rec["job_id"] = (meta or {}).get("_job_id",
                                         spec.get("job_id", "unknown"))
        if spec.get("job_id", "unknown") != "unknown" and \
                rec["job_id"] != spec["job_id"]:
            die(f"config check failed: {run} key job_id = {rec['job_id']} "
                f"!= registry {spec['job_id']}")
        runs[run] = rec
    if missing and not allow_missing:
        print("MISSING PROVENANCE (no case_metadata / banner_config):")
        for r in missing:
            print(f"  {r}")
        die("missing provenance; rerun with --allow-missing-provenance "
            "while harvest is in progress")
    return runs, missing


# ---------------------------------------------------------------------------
# shared physical constants derived from data
# ---------------------------------------------------------------------------

def load_reference_meta():
    meta, src = load_provenance(GAMMA_V_RUN)
    if meta is None:
        die(f"reference run {GAMMA_V_RUN} has no provenance file")
    return meta


def derived_constants():
    meta = load_reference_meta()
    R = float(meta["R"])
    rho = float(meta["rho"])
    nu_eff = float(meta["wake_nu"])
    dt_meta = float(meta["dt"])
    if abs(dt_meta - DT_BASE) > 1e-12:
        die(f"dt cross-check failed: 60/RPM/NT = {DT_BASE!r} vs metadata "
            f"dt = {dt_meta!r}")
    return dict(R=R, rho=rho, nu_eff=nu_eff, dt=DT_BASE, meta=meta)


# ---------------------------------------------------------------------------
# gamma_v
# ---------------------------------------------------------------------------

def compute_gamma_v():
    """Measured peak bound circulation Gamma_v [m^2/s] from
    data/<run>/monitors/<run>_monitor03_bound_circulation_system1.csv.

    Estimator (exact, pre-registered here):
      1. settled window = the last 25% of the distinct steps present in the
         CSV (the file itself may already start mid-run after a restart);
      2. per station (blade, section): median over the window of
         |circulation_te|;
      3. per blade: max of the station medians over sections;
      4. Gamma_v = max over blades.
    Streams with the csv module (files are large); two passes: one for the
    step range, one accumulating windowed values.
    """
    path = os.path.join(DATA, GAMMA_V_RUN, "monitors",
                        f"{GAMMA_V_RUN}_monitor03_bound_circulation_system1.csv")
    if not os.path.isfile(path):
        die(f"gamma_v source missing: {path}")
    steps = set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                steps.add(int(row["step"]))
            except (KeyError, ValueError):
                continue
    if not steps:
        die(f"gamma_v source empty: {path}")
    ordered = sorted(steps)
    cut = ordered[int(math.floor(0.75 * len(ordered)))]
    acc = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                step = int(row["step"])
                if step < cut:
                    continue
                key = (int(row["blade"]), int(row["section"]))
                val = abs(float(row["circulation_te"]))
            except (KeyError, ValueError):
                continue
            if math.isfinite(val):
                acc.setdefault(key, []).append(val)
    per_blade = {}
    for (blade, _sec), vals in acc.items():
        vals.sort()
        n = len(vals)
        med = vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])
        per_blade[blade] = max(per_blade.get(blade, 0.0), med)
    gamma_v = max(per_blade.values())
    return dict(gamma_v=gamma_v, per_blade=per_blade,
                window=(cut, ordered[-1]), n_steps=len(ordered))


# ---------------------------------------------------------------------------
# initializers
# ---------------------------------------------------------------------------

def compute_initializers(consts, gamma_v):
    nu = consts["nu_eff"]
    R = consts["R"]
    sigma_eq = math.sqrt(nu / ZBAR_ANCHOR)
    sigma_stab = math.sqrt(gamma_v * consts["dt"] / (2.0 * math.pi))
    sigma_0 = max(C_EQ * sigma_eq, sigma_stab)
    return dict(sigma_eq=sigma_eq, sigma_stab=sigma_stab, sigma_0=sigma_0,
                sigma_eq_R=sigma_eq / R, sigma_stab_R=sigma_stab / R,
                sigma_0_R=sigma_0 / R)


# ---------------------------------------------------------------------------
# margins
# ---------------------------------------------------------------------------

def wake_health_path(run):
    p = os.path.join(DATA, run, "monitors",
                     f"{run}_monitor04_wake_health_system1.csv")
    return p if os.path.isfile(p) else None


def read_health_series(run, rec):
    """Return list of dicts(step, time, msr, max_u, gos2) or None.

    For the forensics-only run the death_trajectory.csv columns are remapped
    (min_sig_rel -> min_sigma_ratio, max_Gos2 -> max_gamma_over_sigma2; time
    reconstructed as step*dt)."""
    if rec["spec"].get("forensics"):
        p = os.path.join(DATA, f"{run}_forensics", "death_trajectory.csv")
        if not os.path.isfile(p):
            return None
        out = []
        with open(p, newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    step = int(float(row["step"]))
                    out.append(dict(step=step, time=step * rec["dt"],
                                    msr=float(row["min_sig_rel"]),
                                    max_u=float(row["max_u"]),
                                    gos2=float(row["max_Gos2"])))
                except (KeyError, ValueError):
                    continue
        return out or None
    p = wake_health_path(run)
    if p is None:
        return None
    out = []
    with open(p, newline="") as fh:
        rd = csv.DictReader(fh)
        has_dtz = "max_dtZ" in (rd.fieldnames or [])
        for row in rd:
            try:
                d = dict(step=int(row["step"]), time=float(row["time"]),
                         msr=float(row["min_sigma_ratio"]),
                         max_u=float(row["max_u"]),
                         gos2=float(row["max_gamma_over_sigma2"]))
                if has_dtz:
                    d["max_dtZ"] = float(row["max_dtZ"])
                out.append(d)
            except (KeyError, ValueError):
                continue
    return out or None


def p_quantile(sorted_vals, q):
    if not sorted_vals:
        return None
    i = min(len(sorted_vals) - 1,
            max(0, int(math.ceil(q * len(sorted_vals))) - 1))
    return sorted_vals[i]


def score_margin(run, rec, series):
    """Outcome classification + M per the pre-registered estimator.

    Both estimators are always computed when their inputs exist (direct
    `max_dtZ` column; level proxy from min_sigma_ratio) and stored as
    M_direct / M_proxy for the required proxy-vs-direct cross-validation;
    the primary M is direct when available, else proxy (P1)."""
    res = dict(run=run, outcome="stable", t_ignite_step=None,
               t_ignite_time=None, M=None, M_p999=None, method="none",
               M_direct=None, M_proxy=None, M_direct_p999=None,
               M_proxy_p999=None)
    if series is None:
        if rec["spec"].get("known_failed"):
            res["outcome"] = "ignited (provenance)"
        res["method"] = "none (no wake-health CSV)"
        return res
    onset = None
    for d in series:
        if (math.isfinite(d["msr"]) and d["msr"] < 0.0) or \
                d["max_u"] > IGNITE_MAX_U:
            onset = d
            break
    if onset is not None:
        res["outcome"] = "ignited"
        res["t_ignite_step"] = onset["step"]
        res["t_ignite_time"] = onset["time"]
    elif rec["spec"].get("known_failed"):
        res["outcome"] = "ignited (provenance; onset beyond CSV)"
    elif rec["spec"].get("in_flight"):
        res["outcome"] = "survivor (in flight, partial)"
    else:
        res["outcome"] = "survivor"
    # direct estimator when the max_dtZ column exists (added for 019)
    pre = [d["max_dtZ"] for d in series if "max_dtZ" in d
           and math.isfinite(d["max_dtZ"])
           and (onset is None or d["step"] < onset["step"])]
    if pre:
        res["M_direct"] = max(pre)
        res["M_direct_p999"] = p_quantile(sorted(pre), 0.999)
    # level proxy from min_sigma_ratio: valid rows are msr > 0 strictly
    # before ignition onset; M~ = max over downward steps of -dln(msr).
    valid = [d for d in series if math.isfinite(d["msr"]) and d["msr"] > 0.0
             and (onset is None or d["step"] < onset["step"])]
    dlns = []
    for a, b in zip(valid, valid[1:]):
        if b["step"] == a["step"] + 1:
            dlns.append(-(math.log(b["msr"]) - math.log(a["msr"])))
    down = [x for x in dlns if x > 0.0]
    if down:
        res["M_proxy"] = max(down)
        res["M_proxy_p999"] = p_quantile(sorted(dlns), 0.999)
    if res["M_direct"] is not None:
        res["M"] = res["M_direct"]
        res["M_p999"] = res["M_direct_p999"]
        res["method"] = "direct max_dtZ"
    elif res["M_proxy"] is not None:
        res["M"] = res["M_proxy"]
        res["M_p999"] = res["M_proxy_p999"]
        res["method"] = "proxy -dln(min_sr)"
    else:
        res["method"] = "proxy -dln(min_sr) (no valid downward steps)"
    return res


def compute_margins(runs):
    out = {}
    for run in sorted(runs):
        rec = runs[run]
        series = read_health_series(run, rec)
        rec["series"] = series
        out[run] = score_margin(run, rec, series)
    return out


# ---------------------------------------------------------------------------
# proxy-vs-direct cross-validation (required P1 deliverable)
# ---------------------------------------------------------------------------

def compute_xval(runs, margins):
    """Rows for every run carrying BOTH estimators (pre-onset windows)."""
    rows = []
    for run in sorted(runs):
        m = margins[run]
        if m["M_direct"] is None or m["M_proxy"] is None:
            continue
        rows.append(dict(run=run, sigma_R=runs[run]["sigma_over_R"],
                         outcome=m["outcome"],
                         M_direct=m["M_direct"], M_proxy=m["M_proxy"],
                         ratio=m["M_proxy"] / m["M_direct"]))
    return rows


# ---------------------------------------------------------------------------
# p-exponent fit (P1 addendum): M(sigma) ~ sigma^-p, Gamma_implied
# ---------------------------------------------------------------------------

def compute_pfit(runs, margins, consts):
    """Least-squares fit of ln M = a - p ln(sigma/R) per viscosity row,
    over NT=36 screen-class runs with a DIRECT M (proxy points are biased
    low and are excluded from the fit; they still appear on the margin
    figure with their method label).  Only SURVIVOR points are fitted:
    an ignited run's pre-onset M is a divergence transient (measured
    12-28 >> the ignition threshold 2), not a stationary margin, so it
    would corrupt the scaling law; ignited points are reported alongside.
    Gamma_implied = 2*pi*sigma^2*M/dt per point (sigma in meters,
    per-run dt)."""
    R = consts["R"]
    fits = []
    for visc in (False, True):
        pts = []
        for run in sorted(runs):
            rec, m = runs[run], margins[run]
            if rec["cls"] != "screen" or rec["viscous"] != visc:
                continue
            if m["M_direct"] is None:
                continue
            if rec["nt"] != NT_BASE:
                continue    # fid144 has its own dt; excluded from the fit
            sig_m = rec["sigma_over_R"] * R
            gimp = 2.0 * math.pi * sig_m ** 2 * m["M_direct"] / rec["dt"]
            pts.append(dict(run=run, sigma_R=rec["sigma_over_R"],
                            M=m["M_direct"], outcome=m["outcome"],
                            gamma_implied=gimp))
        fitted = [p_ for p_ in pts if p_["outcome"].startswith("survivor")]
        fit = dict(viscous=visc, points=pts, fitted=fitted, p=None, lnA=None)
        if len(fitted) >= 2:
            xs = [math.log(p_["sigma_R"]) for p_ in fitted]
            ys = [math.log(p_["M"]) for p_ in fitted]
            n = len(xs)
            mx, my = sum(xs) / n, sum(ys) / n
            sxx = sum((x - mx) ** 2 for x in xs)
            sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
            slope = sxy / sxx if sxx > 0 else float("nan")
            fit["p"] = -slope
            fit["lnA"] = my - slope * mx
        fits.append(fit)
    # fid144 point, reported alongside (never fitted)
    extra = []
    for run in sorted(runs):
        rec, m = runs[run], margins[run]
        if rec["cls"] == "screen" and rec["nt"] != NT_BASE \
                and m["M_direct"] is not None:
            sig_m = rec["sigma_over_R"] * R
            extra.append(dict(run=run, sigma_R=rec["sigma_over_R"],
                              M=m["M_direct"], outcome=m["outcome"],
                              gamma_implied=2.0 * math.pi * sig_m ** 2
                              * m["M_direct"] / rec["dt"]))
    return dict(fits=fits, off_dt=extra)


# ---------------------------------------------------------------------------
# iterate
# ---------------------------------------------------------------------------

def compute_iterates(runs, margins, init):
    """Procedure table: start at sigma_0, at each iterate read M from the
    nearest probed sigma (any registry run with a finite M), pass/fail vs
    eps, would-be update sigma_{k+1} = sigma_k*sqrt(M_k/eps)."""
    probed = [(runs[r]["sigma_over_R"], r) for r in sorted(runs)
              if margins[r]["M"] is not None]
    rows = []
    sig = init["sigma_0_R"]
    for k in range(6):
        if not probed:
            break
        dist, s_probe, r_probe = min(
            (abs(math.log(s / sig)), s, r) for s, r in probed)
        M = margins[r_probe]["M"]
        ok = M <= EPS
        nxt = sig * math.sqrt(M / EPS)
        rows.append(dict(k=k, sigma_R=sig, probe_run=r_probe,
                         probe_sigma_R=s_probe, M=M,
                         method=margins[r_probe]["method"],
                         verdict="PASS (sigma* = sigma_k)" if ok else "FAIL",
                         sigma_next_R=nxt))
        if ok:
            break
        sig = nxt
    return rows


# ---------------------------------------------------------------------------
# gaps
# ---------------------------------------------------------------------------

def compute_gaps(runs):
    cells = []
    for tgt in GRID_SIGMAS:
        for visc in (False, True):
            match = None
            for r in sorted(runs):
                rec = runs[r]
                if rec["cls"] != "screen" or rec["viscous"] != visc:
                    continue
                if abs(rec["sigma_over_R"] - tgt) / tgt <= GAP_MATCH_RTOL:
                    if match is None or \
                            abs(rec["sigma_over_R"] - tgt) < \
                            abs(runs[match]["sigma_over_R"] - tgt):
                        match = r
            tagbase = f"s{int(round(tgt*1000)):03d}"
            cell = dict(cell=f"{tagbase}_{'visc' if visc else 'inv'}",
                        target=tgt, viscous=visc, covered_by=match)
            if match is None:
                ov, pps = GRID_JOBSPEC[tgt]
                cell.update(ov=ov, pps=pps,
                            actual=tip_arc_sigma_over_R(ov, pps, NT_BASE),
                            case_tag=f"scr_p019_{tagbase}"
                                     + ("v" if visc else ""))
            else:
                cell["actual"] = runs[match]["sigma_over_R"]
            cells.append(cell)
    return cells


def write_gap_manifest(cells):
    path = os.path.join(HERE, "p019_gap_manifest.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cell", "sigma_over_R_target", "sigma_over_R_actual",
                    "OV", "PPS", "viscous", "case_tag"])
        for c in cells:
            if c["covered_by"] is not None:
                continue
            w.writerow([c["cell"], f"{c['target']:.3f}",
                        f"{c['actual']:.5f}", f"{c['ov']:.1f}",
                        str(c["pps"]), str(c["viscous"]).lower(),
                        c["case_tag"]])
    return path


# ---------------------------------------------------------------------------
# figures (CSV data dirs only; .tex authored later)
# ---------------------------------------------------------------------------

def display_label(rec):
    lab = (f"{'screen' if rec['cls']=='screen' else 'prod'}"
           f"-{rec['sigma_over_R']:.4f}"
           f"-{'visc' if rec['viscous'] else 'inv'}"
           f"-NT{rec['nt']}")
    if rec["spec"].get("short"):
        lab += f"-{rec['spec']['short']}"
    return lab


def write_figures(runs, margins, iterates, cells, init, xval, pfit):
    def ensure(d):
        os.makedirs(d, exist_ok=True)
        return d

    # fig_margin_curve
    d = ensure(os.path.join(FIGDIR, "fig_margin_curve"))
    with open(os.path.join(d, "margin_curve.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sigma_over_R", "M", "M_p999", "M_direct", "M_proxy",
                    "method", "class", "viscous", "outcome", "label"])
        for run in sorted(runs):
            rec, m = runs[run], margins[run]
            w.writerow([f"{rec['sigma_over_R']:.5f}",
                        fmt(m["M"], ".6e"), fmt(m["M_p999"], ".6e"),
                        fmt(m["M_direct"], ".6e"), fmt(m["M_proxy"], ".6e"),
                        m["method"], rec["cls"],
                        str(rec["viscous"]).lower(), m["outcome"],
                        display_label(rec)])
    with open(os.path.join(d, "xval_proxy_vs_direct.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "sigma_over_R", "outcome", "M_direct", "M_proxy",
                    "ratio_proxy_over_direct"])
        for r in xval:
            w.writerow([r["run"], f"{r['sigma_R']:.5f}", r["outcome"],
                        f"{r['M_direct']:.6e}", f"{r['M_proxy']:.6e}",
                        f"{r['ratio']:.4f}"])
    with open(os.path.join(d, "pfit.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["viscous", "n_fitted", "p_exponent", "lnA",
                    "runs_fitted"])
        for f in pfit["fits"]:
            w.writerow([str(f["viscous"]).lower(), len(f["fitted"]),
                        fmt(f["p"], ".4f"), fmt(f["lnA"], ".4f"),
                        ";".join(p_["run"] for p_ in f["fitted"])])
    with open(os.path.join(d, "gamma_implied.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "sigma_over_R", "viscous", "nt", "M_direct",
                    "gamma_implied_m2s", "outcome"])
        for f in pfit["fits"]:
            for p_ in f["points"]:
                w.writerow([p_["run"], f"{p_['sigma_R']:.5f}",
                            str(f["viscous"]).lower(), str(NT_BASE),
                            f"{p_['M']:.6e}",
                            f"{p_['gamma_implied']:.6f}", p_["outcome"]])
        for p_ in pfit["off_dt"]:
            w.writerow([p_["run"], f"{p_['sigma_R']:.5f}", "true",
                        str(runs[p_["run"]]["nt"]), f"{p_['M']:.6e}",
                        f"{p_['gamma_implied']:.6f}", p_["outcome"]])
    # plot-ready numeric table for pgfplots (no quoted strings): flags are
    # 0/1 so .tex filters never depend on free-text columns.
    with open(os.path.join(d, "margin_points.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sigma_over_R", "M", "viscous01", "ignited01",
                    "direct01", "screen01", "nt"])
        for run in sorted(runs):
            rec, m = runs[run], margins[run]
            if m["M"] is None:
                continue
            w.writerow([f"{rec['sigma_over_R']:.5f}", f"{m['M']:.6e}",
                        1 if rec["viscous"] else 0,
                        1 if m["outcome"].startswith("ignited") else 0,
                        1 if m["method"].startswith("direct") else 0,
                        1 if rec["cls"] == "screen" else 0,
                        rec["nt"]])
    with open(os.path.join(d, "scales.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "value"])
        w.writerow(["eps", f"{EPS:.3f}"])
        w.writerow(["sigma_eq_over_R", f"{init['sigma_eq_R']:.5f}"])
        w.writerow(["sigma_stab_over_R", f"{init['sigma_stab_R']:.5f}"])
        w.writerow(["sigma_0_over_R", f"{init['sigma_0_R']:.5f}"])

    # fig_regime_map
    d = ensure(os.path.join(FIGDIR, "fig_regime_map"))
    with open(os.path.join(d, "regime_cells.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sigma_over_R", "viscous", "class", "outcome",
                    "t_ignite_step", "t_ignite_time_s", "label"])
        for run in sorted(runs):
            rec, m = runs[run], margins[run]
            w.writerow([f"{rec['sigma_over_R']:.5f}",
                        str(rec["viscous"]).lower(), rec["cls"],
                        m["outcome"], fmt(m["t_ignite_step"], "d"),
                        fmt(m["t_ignite_time"], ".6f"),
                        display_label(rec)])
        for c in cells:
            if c["covered_by"] is None:
                w.writerow([f"{c['actual']:.5f}",
                            str(c["viscous"]).lower(), "screen",
                            "planned (gap)", "-", "-", c["case_tag"]])
    with open(os.path.join(d, "regime_points.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sigma_over_R", "viscous01", "ignited01", "pending01",
                    "t_ign_rev", "screen01", "nt"])
        for run in sorted(runs):
            rec, m = runs[run], margins[run]
            t_rev = (m["t_ignite_step"] / rec["nt"]
                     if m["t_ignite_step"] is not None else float("nan"))
            pending = 1 if (rec["spec"].get("in_flight")
                            and rec.get("series") is None) else 0
            w.writerow([f"{rec['sigma_over_R']:.5f}",
                        1 if rec["viscous"] else 0,
                        1 if m["outcome"].startswith("ignited") else 0,
                        pending,
                        f"{t_rev:.3f}" if math.isfinite(t_rev) else "nan",
                        1 if rec["cls"] == "screen" else 0, rec["nt"]])
    for run in sorted(runs):
        rec, m = runs[run], margins[run]
        if rec.get("series") and m["outcome"].startswith("ignited"):
            p = os.path.join(d, f"traj_{run}.csv")
            with open(p, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["step", "time", "min_sigma_ratio",
                            "max_gamma_over_sigma2", "max_dtZ"])
                for row in rec["series"]:
                    w.writerow([row["step"], f"{row['time']:.8f}",
                                f"{row['msr']:.6e}", f"{row['gos2']:.6e}",
                                fmt(row.get("max_dtZ"), ".6e")])

    # fig_iterates
    d = ensure(os.path.join(FIGDIR, "fig_iterates"))
    with open(os.path.join(d, "iterates.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["k", "sigma_over_R", "probe_run", "probe_sigma_over_R",
                    "M", "method", "verdict", "sigma_next_over_R"])
        for r in iterates:
            w.writerow([r["k"], f"{r['sigma_R']:.5f}", r["probe_run"],
                        f"{r['probe_sigma_R']:.5f}", f"{r['M']:.6e}",
                        r["method"], r["verdict"],
                        f"{r['sigma_next_R']:.5f}"])

    # provenance appendix
    with open(os.path.join(FIGDIR, "provenance_appendix.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "run", "job_id", "provenance_source",
                    "sigma_method", "role"])
        for run in sorted(runs):
            rec = runs[run]
            w.writerow([display_label(rec), run, rec["job_id"],
                        rec["prov_source"] or "MISSING",
                        rec["sigma_method"], rec["role"]])


# ---------------------------------------------------------------------------
# printing
# ---------------------------------------------------------------------------

def print_constants(consts):
    print("=== constants ===")
    print(f"  R          = {consts['R']:.4f} m   "
          f"(from data/{GAMMA_V_RUN}/{GAMMA_V_RUN}_case_metadata.toml)")
    print(f"  RPM        = {RPM:.0f}")
    print(f"  NT         = {NT_BASE}  (base; per-run NT in registry)")
    print(f"  dt         = {DT_BASE:.9e} s  (= 60/RPM/NT; metadata "
          f"cross-check PASSED)")
    print(f"  rho        = {consts['rho']:.4f} kg/m^3  (metadata)")
    print(f"  nu_eff     = {consts['nu_eff']:.6e} m^2/s  (metadata wake_nu, "
          f"AUTHORITATIVE)")
    print(f"               [018 prose carried {NU_PROSE:.4e} m^2/s -- "
          f"discrepancy noted, prose value NOT used]")
    print(f"  Zbar       = {ZBAR_ANCHOR:.2f} s^-1  (018 anchor, "
          f"sigma_blowup_mechanism.md sec.2; recompute pending)")
    print(f"  eps        = {EPS:.2f}  (pre-registered)")
    print(f"  Gamma_v    : measured peak bound circulation from "
          f"{GAMMA_V_RUN} monitor03 (see gamma_v)")


def print_gamma_v(gv):
    print("=== gamma_v ===")
    lo, hi = gv["window"]
    print(f"  source: {GAMMA_V_RUN} monitor03_bound_circulation "
          f"({gv['n_steps']} steps in CSV; settled window steps "
          f"{lo}..{hi} = last 25%)")
    for b in sorted(gv["per_blade"]):
        print(f"  blade {b}: peak-over-stations median |circulation_te| "
              f"= {gv['per_blade'][b]:.6f} m^2/s")
    print(f"  Gamma_v = {gv['gamma_v']:.6f} m^2/s  "
          f"(sanity target ~{GAMMA_V_SANITY} m^2/s)")


def print_initializers(init, consts):
    R = consts["R"]
    print("=== initializers ===")
    print(f"  {'scale':<12} {'meters':>12} {'sigma/R':>10}  definition")
    print(f"  {'sigma_eq':<12} {init['sigma_eq']:>12.6f} "
          f"{init['sigma_eq_R']:>10.5f}  sqrt(nu_eff/Zbar)")
    print(f"  {'sigma_stab':<12} {init['sigma_stab']:>12.6f} "
          f"{init['sigma_stab_R']:>10.5f}  sqrt(Gamma_v*dt/2pi)")
    print(f"  {'sigma_0':<12} {init['sigma_0']:>12.6f} "
          f"{init['sigma_0_R']:>10.5f}  max({C_EQ:.0f}*sigma_eq, sigma_stab)")


def print_ingest(runs, missing):
    print("=== ingest (run registry) ===")
    print(f"  {'run':<24} {'class':<10} {'visc':<5} {'NT':<3} "
          f"{'sigma/R':>8} {'sigma method':<30} {'prov':<14} job")
    for run in sorted(runs):
        rec = runs[run]
        print(f"  {run:<24} {rec['cls']:<10} "
              f"{str(rec['viscous']).lower():<5} {rec['nt']:<3} "
              f"{rec['sigma_over_R']:>8.5f} {rec['sigma_method']:<30} "
              f"{rec['prov_source'] or 'MISSING':<14} {rec['job_id']}")
    if missing:
        print("  MISSING PROVENANCE list (harvest in progress):")
        for r in missing:
            print(f"    {r}")


def print_margins(runs, margins):
    print("=== margins ===")
    print(f"  {'run':<24} {'sigma/R':>8} {'class':<10} {'visc':<5} "
          f"{'outcome':<38} {'t_ign(step)':>11} {'M':>10} {'M_p99.9':>10} "
          f"method")
    for run in sorted(runs):
        rec, m = runs[run], margins[run]
        print(f"  {run:<24} {rec['sigma_over_R']:>8.5f} {rec['cls']:<10} "
              f"{str(rec['viscous']).lower():<5} {m['outcome']:<38} "
              f"{fmt(m['t_ignite_step'], 'd'):>11} "
              f"{fmt(m['M'], '.4e'):>10} {fmt(m['M_p999'], '.4e'):>10} "
              f"{m['method']}")


def print_xval(xval):
    print("=== proxy-vs-direct cross-validation (P1 required) ===")
    if not xval:
        print("  no run carries both estimators yet")
        return
    print(f"  {'run':<24} {'sigma/R':>8} {'outcome':<20} {'M_direct':>10} "
          f"{'M_proxy':>10} {'proxy/direct':>12}")
    for r in xval:
        print(f"  {r['run']:<24} {r['sigma_R']:>8.5f} {r['outcome']:<20} "
              f"{r['M_direct']:>10.4e} {r['M_proxy']:>10.4e} "
              f"{r['ratio']:>12.4f}")
    print("  note: pre-registered expectation is proxy biased LOW "
          "(particle-switch bias); ratio << 1 confirms the direct column "
          "is load-bearing.")


def print_pfit(pfit):
    print("=== p-exponent fit M(sigma) ~ sigma^-p (P1 addendum) ===")
    for f in pfit["fits"]:
        tag = "viscous" if f["viscous"] else "inviscid"
        if f["p"] is None:
            print(f"  {tag}: <2 direct-M NT36 screen SURVIVOR points, fit "
                  f"skipped ({len(f['fitted'])} fittable of "
                  f"{len(f['points'])} with direct M)")
        else:
            print(f"  {tag}: p = {f['p']:.3f} over {len(f['fitted'])} "
                  f"survivor points "
                  f"({', '.join(p_['run'] for p_ in f['fitted'])})")
        for p_ in f["points"]:
            infit = "fitted" if p_ in f["fitted"] else "NOT fitted"
            print(f"    {p_['run']:<24} sigma/R {p_['sigma_R']:.5f}  "
                  f"M {p_['M']:.4e}  Gamma_implied {p_['gamma_implied']:.4f} "
                  f"m^2/s  [{p_['outcome']}; {infit}]")
    for p_ in pfit["off_dt"]:
        print(f"  off-dt point (not fitted): {p_['run']} sigma/R "
              f"{p_['sigma_R']:.5f}  M {p_['M']:.4e}  Gamma_implied "
              f"{p_['gamma_implied']:.4f} m^2/s  [{p_['outcome']}]")
    print("  interpretation (pre-registered): p ~ 2 -> compact sub-sigma "
          "strainer (sigma_stab regime, Gamma_implied ~ Gamma_v); "
          "p ~ 0 -> ambient-dominated (resolved).")


def print_iterates(iterates):
    print("=== iterate ===")
    print(f"  {'k':<2} {'sigma_k/R':>9} {'probe run':<24} "
          f"{'probe s/R':>9} {'M':>10} {'verdict':<24} {'sigma_k+1/R':>11}")
    for r in iterates:
        print(f"  {r['k']:<2} {r['sigma_R']:>9.5f} {r['probe_run']:<24} "
              f"{r['probe_sigma_R']:>9.5f} {r['M']:>10.4e} "
              f"{r['verdict']:<24} {r['sigma_next_R']:>11.5f}")
    print("  note: sigma_k+1 is the would-be update sigma_k*sqrt(M_k/eps); "
          "on PASS it is informational only.")


def print_gaps(cells, manifest_path):
    print("=== gaps (target grid, screen class) ===")
    print(f"  {'cell':<10} {'target':>7} {'visc':<5} {'status':<40}")
    for c in cells:
        if c["covered_by"] is not None:
            st = (f"covered by {c['covered_by']} "
                  f"(sigma/R {c['actual']:.5f})")
        else:
            st = (f"GAP -> OV{c['ov']:.1f}/PPS{c['pps']} "
                  f"(sigma/R {c['actual']:.5f}) tag {c['case_tag']}")
        print(f"  {c['cell']:<10} {c['target']:>7.3f} "
              f"{str(c['viscous']).lower():<5} {st}")
    print(f"  gap manifest written: {manifest_path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("subcommand", nargs="?", default="all",
                    choices=["all", "constants", "initializers", "gamma_v",
                             "ingest", "margins", "iterate", "gaps",
                             "figures"])
    ap.add_argument("--allow-missing-provenance", action="store_true",
                    help="do not exit nonzero on runs lacking provenance "
                         "(harvest in progress)")
    args = ap.parse_args()
    cmd = args.subcommand

    consts = derived_constants()
    if cmd in ("all", "constants"):
        print_constants(consts)
    if cmd == "constants":
        return

    need_registry = cmd in ("all", "ingest", "margins", "iterate", "gaps",
                            "figures")
    need_gv = cmd in ("all", "gamma_v", "initializers", "iterate", "figures")

    gv = compute_gamma_v() if need_gv else None
    if cmd in ("all", "gamma_v"):
        print_gamma_v(gv)
    if cmd == "gamma_v":
        return

    init = compute_initializers(consts, gv["gamma_v"]) if gv else None
    if cmd in ("all", "initializers"):
        print_initializers(init, consts)
    if cmd == "initializers":
        return

    runs, missing = build_registry(args.allow_missing_provenance,
                                   consts["R"]) \
        if need_registry else ({}, [])
    if cmd in ("all", "ingest"):
        print_ingest(runs, missing)
    if cmd == "ingest":
        return

    margins = compute_margins(runs)
    xval = compute_xval(runs, margins)
    pfit = compute_pfit(runs, margins, consts)
    if cmd in ("all", "margins"):
        print_margins(runs, margins)
        print_xval(xval)
        print_pfit(pfit)
    if cmd == "margins":
        return

    if cmd in ("all", "iterate", "figures"):
        iterates = compute_iterates(runs, margins, init)
    if cmd in ("all", "iterate"):
        print_iterates(iterates)
    if cmd == "iterate":
        return

    cells = compute_gaps(runs)
    if cmd in ("all", "gaps"):
        manifest = write_gap_manifest(cells)
        print_gaps(cells, manifest)
    if cmd == "gaps":
        return

    write_figures(runs, margins, iterates, cells, init, xval, pfit)
    print(f"=== figures ===")
    print(f"  CSV data dirs written under {FIGDIR}")
    if missing:
        print("REMINDER -- MISSING PROVENANCE:")
        for r in missing:
            print(f"  {r}")


if __name__ == "__main__":
    main()

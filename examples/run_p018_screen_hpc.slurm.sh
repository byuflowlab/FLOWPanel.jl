#!/usr/bin/env bash
#SBATCH --job-name=fp-018-screen
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# BRAINSTORM/018 phase_14 SCREENING runner (Ryan 2026-08-05: short time-marched
# runs for parameter SEARCH; slot limits waived for the screen ladders).
#
# One 8-rev case per job on the PRODUCTION mesh, freestream pulse ZEROED
# (1-rev spinup only). The per-step monitors give every stage window (revs
# 1-2 / 2-4 / 4-8) from a single job — no warm chaining. PRE-REGISTERED
# EVIDENCE CLASS: screen results rank parameters and prune ladders; they can
# NEVER satisfy M1/M2 (those need >=15 settled revs on the standard staged
# startup). Screen CT is not a CT prediction (no freestream pulse); only
# rung-to-rung deltas at matched windows are meaningful.
#
# Submit ONE case per job, from the top level of the cluster checkout:
#   sbatch --job-name=fp-018-scr-<tag> examples/run_p018_screen_hpc.slurm.sh scr_<tag>
# Case arms are UNCONDITIONAL exports (ops rule: no ${VAR:-} defaults on
# defining knobs). Verify every job's banner after submission (mandatory).

set -euo pipefail
THREADS=64
# 026 Phase 1b: P018_REPO_OVERRIDE submits from a git worktree carrying code
# the live checkout must not receive (pair with P018_PROJECT_OVERRIDE below,
# and give the worktree a `data` symlink to the shared data root).
EXPECTED_REPO="${P018_REPO_OVERRIDE:-/home/rander39/projects/FLOWPanel.jl}"
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

CASE="${1:-}"
[[ -n "$CASE" ]] || { echo "ERROR: pass a case tag, e.g. scr_b0" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg

# BRAINSTORM/025 (2026-08-20): the codebase default filament regularization
# changed to Gaussian. Pin this campaign to the legacy Vatistas family so
# results stay comparable; override at submission via --export if intended.
export FLOWPANEL_FILAMENT_REG="${FLOWPANEL_FILAMENT_REG:-vatistas}"
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

# ---- Fixed screen configuration (identical for every case) ------------------
export RHPC_MESH=45_185_ct4          # mesh-ladder cases override below
export RPM=5400
export RHPC_FORMULATION=velocity
export BERNOULLI_ONLY=true
export RUN_KJ=false
export SAVE_VTK=true                 # false suppresses ALL monitors/CSVs
export TRUNCATION_DEPTH_R=4
export PARTICLE_SHEDDING=sigma_overlap
export DAS_KINEMATIC_ARC=false
# Screen schedule: 1-rev spinup, NO freestream pulse, 8 revs total.
export SPINUP_REVS=1
export SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0
export MAGVINF_PEAK=0.0
export MAGVINF_END=0.0
export FREESTREAM_RAMP_REVS=0
export FREESTREAM_HOLD_REVS=0
export FREESTREAM_WITHDRAW_REVS=0
export SETTLE_REVS=0
export NREVS=8
export CONVERGENCE_REVS=2
export CONVERGENCE_MEAN_TOL=0.005
export CONVERGENCE_PTP_TOL=0.02
# Defaults per B0 carrier; case arms override unconditionally.
export NT=36
export RELAX_RLXF=0.3

case "$CASE" in
  # ---- Calibration ladder: known HPC verdicts re-measured at screen length --
  # (same production mesh, so calibration is purely about run length)
  scr_b0)     export OVERLAP=2.0; export P_PER_STEP=4; export MERGE_R_FACTOR=0.0120; export NWAKEROWS=4; export DAS_ETA_KINEMATIC=1.0 ;;
  scr_ov1p0)  export OVERLAP=1.0; export P_PER_STEP=4; export MERGE_R_FACTOR=0.0060; export NWAKEROWS=4; export DAS_ETA_KINEMATIC=1.0 ;;
  scr_nrows1) export OVERLAP=2.0; export P_PER_STEP=4; export MERGE_R_FACTOR=0.0120; export NWAKEROWS=1; export DAS_ETA_KINEMATIC=1.0 ;;
  scr_nrows2) export OVERLAP=2.0; export P_PER_STEP=4; export MERGE_R_FACTOR=0.0120; export NWAKEROWS=2; export DAS_ETA_KINEMATIC=1.0 ;;
  scr_nt72)   export NT=72; export OVERLAP=2.0; export P_PER_STEP=2; export MERGE_R_FACTOR=0.0120; export NWAKEROWS=4; export DAS_ETA_KINEMATIC=2.0; export RELAX_RLXF=0.15 ;;

  # ---- Campaign A: uniform-d_front map at sigma=0.03R (OV 2.4/PPS 14) -------
  # N=2 feasibility floor at this sigma is D>=6.45 (phase_13 s4b) — no N=2 arm
  # below it.
  scr_uf_n1_d2p6) export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=2.6 ;;
  scr_uf_n1_d3p4) export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4 ;;
  scr_uf_n1_d5p0) export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=5.0 ;;
  scr_uf_n2_d6p5) export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=2; export DAS_UNIFORM_DSIGMA=6.5 ;;
  scr_uf_n2_d8p5) export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=2; export DAS_UNIFORM_DSIGMA=8.5 ;;

  # ---- Campaign B: dt ladder under the ufront law (first CLEAN dt ladder —
  # d_front is dt-independent by construction). sigma pinned at 0.0291R via
  # OVERLAP 2.0 with PPS 12/6/3; rlxf halves with dt.
  scr_ufdt_nt36)  export NT=36;  export OVERLAP=2.0; export P_PER_STEP=12; export MERGE_R_FACTOR=0.0040; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export RELAX_RLXF=0.3 ;;
  scr_ufdt_nt72)  export NT=72;  export OVERLAP=2.0; export P_PER_STEP=6;  export MERGE_R_FACTOR=0.0040; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export RELAX_RLXF=0.15 ;;
  scr_ufdt_nt144) export NT=144; export OVERLAP=2.0; export P_PER_STEP=3;  export MERGE_R_FACTOR=0.0040; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export RELAX_RLXF=0.075 ;;

  # ---- Campaign C: sigma ladder at ufront N=1 D=3.4 (0.03R rung = scr_uf_n1_d3p4).
  # 0.02R is the sigma_eq rung that blew up at long horizons; the screen reads
  # the WakeHealthMonitor contraction RATE (ignition ~step 1250 is beyond the
  # 8-rev horizon — pre-registered limitation).
  scr_ufsig_0p05) export OVERLAP=2.0; export P_PER_STEP=7;  export MERGE_R_FACTOR=0.0069;  export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4 ;;
  scr_ufsig_0p02) export OVERLAP=2.4; export P_PER_STEP=21; export MERGE_R_FACTOR=0.00275; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4 ;;

  # ---- Campaign D: mesh ladder at the ufront carrier (anchor 45_185_ct4 =
  # scr_uf_n1_d3p4). Chordwise anchor'd by the 2c/2d HPC convergence (to 249);
  # spanwise untested (Phase 10). ct4 cap recipe ONLY.
  scr_mesh_c145) export RHPC_MESH=45_145_ct4; export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4 ;;
  scr_mesh_c249) export RHPC_MESH_FILE=dji9443_20260725_45_249_capped_captess4.msh; export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4 ;;
  scr_mesh_s60)  export RHPC_MESH_FILE=dji9443_20260731_60_185_capped_captess4.msh; export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4 ;;
  scr_mesh_s80)  export RHPC_MESH_FILE=dji9443_20260731_80_185_capped_captess4.msh; export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4 ;;

  # ---- Campaign E (BRAINSTORM 019): sigma regime map, {sigma} x {inviscid,
  # viscous}, ufront N=1 carrier. D=3.4 (C1 bound) at EVERY rung, incl.
  # 0.038R (Ryan 2026-08-06: C1 clearance outranks the tip chord-band cap;
  # the s038 pair therefore exceeds 1.5c at the tip — labeled attribute).
  # MERGE_R_FACTOR = 0.1375*sigma/R (Campaign A/C convention). Viscous arms:
  # CoreSpreading, beta 1e9, WAKE_NU default (matches p018_*_visc). All arms
  # carry WAKE_HEALTH_DTZ=true (max_dtZ tripwire column, default-off in src).
  # Inviscid 0.020R reuses scr_ufsig_0p02; inviscid 0.030R reuses
  # scr_uf_n1_d3p4 / scr_ufdt_nt36 — no new arms for those cells.
  scr_p019_s015)  export OVERLAP=2.4; export P_PER_STEP=28; export MERGE_R_FACTOR=0.00206; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true ;;
  scr_p019_s015v) export OVERLAP=2.4; export P_PER_STEP=28; export MERGE_R_FACTOR=0.00206; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;
  scr_p019_s020v) export OVERLAP=2.4; export P_PER_STEP=21; export MERGE_R_FACTOR=0.00275; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;
  scr_p019_s025)  export OVERLAP=2.0; export P_PER_STEP=14; export MERGE_R_FACTOR=0.00343; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true ;;
  scr_p019_s025v) export OVERLAP=2.0; export P_PER_STEP=14; export MERGE_R_FACTOR=0.00343; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;
  scr_p019_s030v) export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041;  export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;
  scr_p019_s038)  export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true ;;
  scr_p019_s038v) export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;
  scr_p019_s038v_cpu40) export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=40 ;;

  # ---- Campaign E addendum (019 fidelity discriminator, Ryan 2026-08-06):
  # NT=144 viscous run INSIDE the fidelity band (sigma = 0.0209R = 1.18
  # sigma_eq = 0.59*2sigma_eq) but ABOVE the Dt/4 stability scale
  # (1.35*sigma_stab(NT144) = 0.0155R). Pre-registered: stability arm says it
  # survives with M <= eps; fidelity assertion says its core-size field
  # migrates to the local Burgers equilibrium instead of staying at shed
  # sigma. Submit with --time=36:00:00 (1152 steps).
  scr_p019_fid144) export NT=144; export OVERLAP=2.4; export P_PER_STEP=5; export MERGE_R_FACTOR=0.00288; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export RELAX_RLXF=0.075; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;

  # ---- Campaign E addendum 2 (Ryan 2026-08-06): probe at exactly sigma_stab
  # (0.03105R predicted; OV2.5/PPS14 -> sigma/R 0.03117, +0.4%) to answer
  # "would sigma_stab alone work as the initializer" by measurement.
  scr_p019_sstab) export OVERLAP=2.5; export P_PER_STEP=14; export MERGE_R_FACTOR=0.00429; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;

  # ---- BRAINSTORM 020 Phase 2 (2026-08-06, Ryan-authorized this session):
  # Stage B — exponential (exact-in-Z) integrator prototype at sigma/R=0.02
  # viscous (clone of scr_p019_s020v + WAKE_EXPINT=true; closure does not
  # exist yet). Pre-registered readouts in phase_01_theory.md sec.7.2.
  scr_p020_exp_s020v) export OVERLAP=2.4; export P_PER_STEP=21; export MERGE_R_FACTOR=0.00275; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export WAKE_EXPINT=true ;;
  # Phase-2R corrected frozen-gradient geometric rerun. Kept under a new tag
  # so the superseded 13065481 artifacts remain available for comparison.
  scr_p020r_geom_s020v) export OVERLAP=2.4; export P_PER_STEP=21; export MERGE_R_FACTOR=0.00275; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export WAKE_EXPINT=true ;;
  # Leg 3 — SFS viscous-null twin: clone of scr_p019_s025v + SFS_OFF=true
  # (no SFS-off run exists anywhere in data/; sec.7.3 fallback pair, with
  # scr_p019_s025v itself as the SFS-on member).
  scr_p020_s025v_nosfs) export OVERLAP=2.0; export P_PER_STEP=14; export MERGE_R_FACTOR=0.00343; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export SFS_OFF=true ;;

  # ---- Illustration render reruns (Ryan 2026-08-27): exact clones of four
  # earlier arms whose particle VTPs were swept, resubmitted under _rr tags so
  # the original run dirs (and their 019/020 monitor CSVs, which figures cite)
  # are never touched. Purpose is ParaView animation frames only — the cluster
  # dep stack has moved since the originals (julia 1.12.6 env sync 08-24), so
  # ignition steps may shift; do not score these against 019/020 verdicts.
  scr_p019_s025_rr)       export OVERLAP=2.0; export P_PER_STEP=14; export MERGE_R_FACTOR=0.00343; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true ;;
  scr_p019_s025v_rr)      export OVERLAP=2.0; export P_PER_STEP=14; export MERGE_R_FACTOR=0.00343; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;
  scr_p019_s030v_rr)      export OVERLAP=2.4; export P_PER_STEP=14; export MERGE_R_FACTOR=0.0041;  export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9 ;;
  scr_p020r_geom_s020v_rr) export OVERLAP=2.4; export P_PER_STEP=21; export MERGE_R_FACTOR=0.00275; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export WAKE_EXPINT=true ;;

  # ---- BRAINSTORM/026 Phase 1 no-split discriminator arms (2026-09-03):
  # warm-start continuations of the scr_p019_s038v_gpu40 stretching ignition
  # (patient zero idx 102340, growth 850-998, ignition 985-1010) from step 950.
  # Clones of scr_p019_s038v_gpu40 minus its device pfield: CPU arms (driver
  # default VPM_ARRAYTYPE=array) because euler_exp has no GPU path; all four
  # share the backend so arms differ by exactly one knob. NREVS=29.5833 ends
  # the continuation at total step 1100 (>= design gate step 1070). Submit with
  # RESTART_STEP=950 RESTART_NAME=scr_p019_s038v_gpu40
  # RESTART_PATH=data/p026_restart_gpu40_s950 (extracted from the archive
  # tarball; see BRAINSTORM/026 design doc sec.9-10).
  scr_p026ph1_ctrl_gpu40)    export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=29.5833 ;;
  scr_p026ph1_exp_gpu40)     export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=29.5833; export WAKE_EXPINT=true ;;
  scr_p026ph1_proj_gpu40)    export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=29.5833; export SFS_NO_BACKSCATTER_PROJECT=true ;;
  scr_p026ph1_expproj_gpu40) export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=29.5833; export WAKE_EXPINT=true; export SFS_NO_BACKSCATTER_PROJECT=true ;;

  # ---- 026 Phase 1 LineGauss-event arms (2026-09-03): same discriminators
  # replayed on the independent LineGauss ignition (patient zero idx 160932,
  # root region, ignition ~490-516) from step 450. Source run = the LG twin
  # of scr_p019_s038v_gpu40 (SAME run name, 052-h200 silo, archived as
  # __052-h200; LineGauss filament reg, rlxf 0.3, otherwise identical env).
  # NREVS=15.6944 ends at total step 600. Submit with RESTART_STEP=450
  # RESTART_NAME=scr_p019_s038v_gpu40 RESTART_PATH=data/p026_restart_lg_s450.
  scr_p026ph1_ctrl_lg)    export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=15.6944; export FLOWPANEL_FILAMENT_REG=linegauss ;;
  scr_p026ph1_exp_lg)     export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=15.6944; export FLOWPANEL_FILAMENT_REG=linegauss; export WAKE_EXPINT=true ;;
  scr_p026ph1_proj_lg)    export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=15.6944; export FLOWPANEL_FILAMENT_REG=linegauss; export SFS_NO_BACKSCATTER_PROJECT=true ;;
  scr_p026ph1_expproj_lg) export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=15.6944; export FLOWPANEL_FILAMENT_REG=linegauss; export WAKE_EXPINT=true; export SFS_NO_BACKSCATTER_PROJECT=true ;;

  # ---- 026 Phase 1b Task 2 RK3 discriminator arms (2026-09-04): same two
  # ignition continuations as the Phase 1 ctrl arms, stepping the wake with
  # low-storage RK3 (WAKE_INTEGRATOR=rk3; panel solve frozen within the step,
  # particle UJ+SFS re-evaluated per stage, ~3x UJ cost). Prediction on
  # record (design doc sec.14): RK3's dtZ ceiling ~0.84 vs Euler's 2/3 ->
  # expect DELAY, NOT ARREST. Submit exactly as the Phase 1 arms
  # (RESTART_STEP/RESTART_NAME/RESTART_PATH), from the Phase-1b worktree via
  # P018_REPO_OVERRIDE/P018_PROJECT_OVERRIDE.
  scr_p026ph1_rk3_gpu40) export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=29.5833; export WAKE_INTEGRATOR=rk3 ;;
  scr_p026ph1_rk3_lg)    export OVERLAP=2.4; export P_PER_STEP=11; export MERGE_R_FACTOR=0.00524; export NWAKEROWS=1; export DAS_UNIFORM_DSIGMA=3.4; export WAKE_HEALTH_DTZ=true; export WAKE_HEALTH_ATTRIBUTION=true; export CORE_SPREADING_ACTIVE=true; export WAKE_CORE_BETA=1e9; export NREVS=15.6944; export FLOWPANEL_FILAMENT_REG=linegauss; export WAKE_INTEGRATOR=rk3 ;;

  *) echo "ERROR: unknown screen case '$CASE'" >&2; exit 2 ;;
esac

export RUN_NAME="$CASE"
# Restart chaining (Ryan 2026-08-29, after 13508968 hit its 12 h wall at step
# 423/1475). With RESTART_STEP >= 0 the driver resumes from the on-disk state at
# that step and APPENDS to the same VTK series, so the existing run dir MUST
# survive — never wipe or rename it in this mode. RESTART_NAME/RESTART_PATH
# default to this run's own dir; override only to chain from a different tag.
if [[ "${RESTART_STEP:--1}" -ge 0 ]]; then
  export RESTART_STEP
  export RESTART_NAME="${RESTART_NAME:-$RUN_NAME}"
  export RESTART_PATH="${RESTART_PATH:-data/$RESTART_NAME}"
  [[ -d "$RESTART_PATH" ]] || { echo "ERROR: RESTART_STEP=$RESTART_STEP but $RESTART_PATH does not exist" >&2; exit 2; }
  probe="$RESTART_PATH/${RESTART_NAME}_body1/${RESTART_NAME}_body1.${RESTART_STEP}.vtu"
  [[ -f "$probe" ]] || { echo "ERROR: restart state missing: $probe" >&2; exit 2; }
  echo "RESTART MODE: resuming $RUN_NAME from step $RESTART_STEP of $RESTART_PATH (run dir PRESERVED)"
elif [[ -d "data/$RUN_NAME" ]]; then
  if [[ "${RHPC_KEEP_PREV:-false}" == "true" ]]; then
    rm -rf "data/$RUN_NAME.prev"; mv "data/$RUN_NAME" "data/$RUN_NAME.prev"
  else
    rm -rf "data/$RUN_NAME"
  fi
fi

echo "BRAINSTORM/018 phase_14 SCREEN — case $CASE (8 revs, no freestream pulse)"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
echo "  mesh:${RHPC_MESH_FILE:-$RHPC_MESH} formulation:$RHPC_FORMULATION RPM:$RPM NT:$NT depth:${TRUNCATION_DEPTH_R}R rlxf:$RELAX_RLXF"
echo "  das_eta:${DAS_ETA_KINEMATIC:-nan} overlap:${OVERLAP} pps:${P_PER_STEP} merge_r:${MERGE_R_FACTOR} nrevs:${NREVS} spinup:${SPINUP_REVS}"
echo "  nwakerows:${NWAKEROWS} das_chord:${DAS_CHORD_FRACTION:-nan} das_uniform:${DAS_UNIFORM_DSIGMA:-nan} visc:${CORE_SPREADING_ACTIVE:-false} expint:${WAKE_EXPINT:-false} integrator:${WAKE_INTEGRATOR:-euler-family} sfs_off:${SFS_OFF:-false}"
echo "  started $(date '+%F %T')"

julia --project="${P018_PROJECT_OVERRIDE:-.}" -t "$THREADS" examples/rotor_hover_pressure_comparison.jl

echo "Screen case $CASE finished ($(date '+%F %T'))."
echo "Artifacts: data/$RUN_NAME/${RUN_NAME}_CT_vs_rev.csv + monitors/"

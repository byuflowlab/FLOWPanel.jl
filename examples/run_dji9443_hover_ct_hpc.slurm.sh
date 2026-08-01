#!/usr/bin/env bash
#SBATCH --job-name=fp-p2e-hoverct
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm/slurm-%x-%j.out
#SBATCH --error=logs/slurm/slurm-%x-%j.err

# DJI 9443 Phase 2e — unsteady hover CT convergence on the Phase 2d production
# mesh (flat root cap + round CapUMinTess=4 tip cap, 45 span, n_airfoil=185),
# PanelParticleWake with nwakerows=1, velocity vs Green-reconstruction solve
# formulations, converging timestep (NT) and wake truncation depth.
#
# Submit ONE case per job, from the top level of the cluster checkout:
#
#   sbatch --job-name=fp-p2e-vel-nt36-d4   examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_vel_nt36_d4
#   sbatch --job-name=fp-p2e-green-nt36-d4 examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_green_nt36_d4
#   sbatch --job-name=fp-p2e-vel-nt72-d4 --time=48:00:00 \
#          examples/run_dji9443_hover_ct_hpc.slurm.sh p2e_vel_nt72_d4
#
# Study directive: at most SIX active study jobs at a time (raised from three by
# Ryan on 2026-07-25; check `squeue` before submitting).
# The NT=72 case is 2x the steps of NT=36; give it a longer --time on submission.
#
# MEMORY: the production mesh has 36,752 triangular panels, so the dense
# operators are large (measured by examples/ preflight, not guessed):
#   Backslash G (lu! in place) ........................ 10.1 GB
#   GreenReconstruction :area_mean build peak (B + K) .. 20.1 GB (K retained: 10.1 GB)
#   => velocity case steady ~10 GB; green case peak ~30 GB
# --mem=64G covers both with headroom for the particle wake and FMM trees. A
# velocity-only case can be submitted leaner with `sbatch --mem=24G ...`.
#
# `logs/slurm/` must already exist before submission (Slurm opens the log paths
# before the script runs).

set -euo pipefail
THREADS=64
EXPECTED_REPO=/home/rander39/projects/FLOWPanel.jl
[[ "$PWD" == "$EXPECTED_REPO" ]] || { echo "ERROR: submit from $EXPECTED_REPO; current dir is $PWD" >&2; exit 2; }

CASE="${1:-}"
[[ -n "$CASE" ]] || { echo "ERROR: pass a case tag, e.g. p2e_vel_nt36_d4" >&2; exit 2; }

export JULIA_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
       BLAS_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" MPLBACKEND=Agg
# Non-interactive submission shells don't put julia on PATH (Manifest pins 1.11.7);
# fall back to the site spack julia-1.11.7 binary (the shared /apps/juliaup launcher is broken).
command -v julia >/dev/null 2>&1 || \
  export PATH="/apps/spack/root/opt/spack/linux-rhel9-haswell/gcc-13.2.0/julia-1.11.7-6bmogflhr2w6mi2zerinukr2gpnpr2rs/juliaup/julia-1.11.7+0.x64.linux.gnu/bin:$PATH"

# ---- Fixed Phase 2e configuration (identical for every case) ----------------
# Mesh: Phase 2d production recipe. NEVER substitute a round-ct3 "capped" mesh.
export RHPC_MESH=45_185_ct4
export RPM=5400
export BERNOULLI_ONLY=true          # PressureBernoulli(unsteady=false) -> ForceMonitor
export RUN_KJ=false
export SAVE_VTK=true                # HPC policy: leave VTK on
# Validated staged-startup recipe (Item 005 stable-wake baseline).
export SPINUP_REVS=1.5
export SPINUP_START_FRACTION=0.4
export MAGVINF_START=0.0
export MAGVINF_PEAK=5.0
export MAGVINF_END=0.0
export FREESTREAM_RAMP_REVS=1.0
export FREESTREAM_HOLD_REVS=1.5
export FREESTREAM_WITHDRAW_REVS=4
export SETTLE_REVS=12               # => ~20 revs total; required_revs drives length
export NREVS=10                     # below the schedule total; schedule wins
# Phase 2e convergence metric (final 5 revs; see the phase plan).
export CONVERGENCE_REVS=5
export CONVERGENCE_MEAN_TOL=0.005
export CONVERGENCE_PTP_TOL=0.02

# ---- Das construction -------------------------------------------------------
# The driver defaults DAS_KINEMATIC_ARC=true (offset follows the arc the TE
# sweeps). Every Phase 2e case through 2026-07-28 -- the whole eta ladder and
# both filter cases -- ran the LEGACY TANGENT construction, so it is pinned false
# here to keep that ladder reproducible and internally comparable. Cases that
# deliberately test the arc construction set it back to true individually.
export DAS_KINEMATIC_ARC=false

# ---- Per-case knobs ---------------------------------------------------------
# Stock corrected-Pedrizzetti rlxf is 0.3 at NT=36. At NT=72 the timestep halves,
# so the per-step factor is halved to preserve the physical relaxation rate.
case "$CASE" in
  p2e_vel_nt36_d4)   export RHPC_FORMULATION=velocity; export NT=36; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3  ;;
  p2e_green_nt36_d4) export RHPC_FORMULATION=green;    export NT=36; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3  ;;
  p2e_vel_nt72_d4)   export RHPC_FORMULATION=velocity; export NT=72; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.15 ;;
  p2e_green_nt72_d4) export RHPC_FORMULATION=green;    export NT=72; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.15 ;;
  p2e_vel_nt36_d6)   export RHPC_FORMULATION=velocity; export NT=36; export TRUNCATION_DEPTH_R=6; export RELAX_RLXF=0.3  ;;
  p2e_green_nt36_d6) export RHPC_FORMULATION=green;    export NT=36; export TRUNCATION_DEPTH_R=6; export RELAX_RLXF=0.3  ;;
  p2e_vel_nt36_d8)   export RHPC_FORMULATION=velocity; export NT=36; export TRUNCATION_DEPTH_R=8; export RELAX_RLXF=0.3  ;;

  # ---- Batch 2: near-wake relaxation filter, 2x2 factorial in (sigma, Das) ----
  # All: green, NT=36, truncation 4R, relaxation DISABLED inside 2.0R of the
  # rotor plane (BRAINSTORM/006: the depth curve saturates at 2.0R, giving the
  # best settled CT 0.0686 vs 0.0640 unfiltered).
  #
  # The filtered wake settles into a ~9-rev LIMIT CYCLE, not a plateau, so the
  # run is extended to 20 settle revs (~26.5 total) and judged on a 10-rev
  # cycle-mean. The 5-rev criterion cannot fire on a limit cycle by construction.
  #
  # sigma: measured on job 12894164 the stock config gives particle
  # sigma/R = 0.022 (min) / 0.145 (median) / 0.379 (max) -- the tip vortex is
  # smeared over a third of the radius. sigma = dist*OVERLAP/P_PER_STEP with
  # dist ~ one step's TE travel, so sigC = stock (3/2), sigF = 2/8 -> ~6x finer
  # (sigma_med ~0.024R, ~230k particles, ~4x cost). MERGE_R_FACTOR must shrink
  # with sigma or merging undoes the refinement.
  # MEMORY: the relaxation filter makes the particle count NEVER SATURATE. The
  # unfiltered baseline levels at ~44k particles by step 719 (the truncation
  # cylinder sweeps them out); with relaxation off inside 2R the particles are
  # not realigned, scatter, and linger inside the cylinder -- 56k by step 420
  # (sigC) and 118k by step 219 (sigF), both still climbing when they OOM'd.
  # All four Batch-2 cases died OUT_OF_MEMORY. CORRECTED 2026-07-27: this is
  # DIVERGENCE, not a memory ceiling. The OOM is thrown from merge_particles!
  # once the blown-up cloud wrecks the spatial hash; ParaView confirms the wake
  # exploded, and resubmitting at 384 GB failed identically (12915036/37), so
  # memory is the symptom. An earlier reading of "sane max|CT|" was invalid --
  # buffered stdout is discarded on SIGKILL, so a killed run's log stops well
  # before the failure.
  #
  # OVERLAP=2 WAS TRIED AND DIVERGED (job 12894481: CT ran away to ~90, then
  # OOM'd inside merge_particles! when the exploded cloud wrecked the spatial
  # hash). Ryan's stability fallback OVERLAP=4 is therefore the default here:
  # sigma_med ~0.048R, still 3x finer than stock, particle count unchanged
  # (p_per_step sets count; overlap only scales sigma). MERGE_R_FACTOR is
  # scaled to hold r_merge/sigma at the stock ratio (0.02/0.145 = 0.138).
  p2e_green_f2p0_sigC_das0p2)
    export RHPC_FORMULATION=green; export NT=36; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export RELAX_FILTER_DOWNSTREAM_R=2.0; export SETTLE_REVS="${SETTLE_REVS:-20}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=0.2 ;;
  p2e_green_f2p0_sigF_das0p2)
    export RHPC_FORMULATION=green; export NT=36; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export RELAX_FILTER_DOWNSTREAM_R=2.0; export SETTLE_REVS="${SETTLE_REVS:-20}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=0.2
    export OVERLAP="${OVERLAP:-4}"; export P_PER_STEP=8; export MERGE_R_FACTOR="${MERGE_R_FACTOR:-0.0066}" ;;
  p2e_green_f2p0_sigC_das4p0)
    export RHPC_FORMULATION=green; export NT=36; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export RELAX_FILTER_DOWNSTREAM_R=2.0; export SETTLE_REVS="${SETTLE_REVS:-20}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=4.0 ;;
  p2e_green_f2p0_sigF_das4p0)
    export RHPC_FORMULATION=green; export NT=36; export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export RELAX_FILTER_DOWNSTREAM_R=2.0; export SETTLE_REVS="${SETTLE_REVS:-20}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=4.0
    export OVERLAP="${OVERLAP:-4}"; export P_PER_STEP=8; export MERGE_R_FACTOR="${MERGE_R_FACTOR:-0.0066}" ;;

  # ---- Staged studies (BRAINSTORM/006 "Staged studies (2026-07-27)") --------
  # Carrier for ALL of these is the UNFILTERED config (the only one that
  # completes at this mesh size) at NT=36, truncation 4R, SETTLE_REVS=12 ->
  # 719 steps, matching the completed baselines for like-for-like comparison.
  # Formulation defaults to `velocity`: Green sits +0.9% (inside cycle scatter)
  # and costs a second ~10 GB dense operator, so it is not worth carrying
  # through a discretization sweep -- re-apply it to the winning config.
  # Override any of RHPC_FORMULATION / SETTLE_REVS / NT via the environment.

  # (A) Das length convergence. Clean single-parameter sweep: sigma is set by
  # per-step TE travel and is INDEPENDENT of Das. Goal is the eta at which CT
  # stops changing -- NOT the eta nearest 0.072.
  p2e_das0p2|p2e_das0p5|p2e_das1p0|p2e_das2p0|p2e_das4p0)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    case "$CASE" in
      p2e_das0p2) export DAS_ETA_KINEMATIC=0.2 ;;
      p2e_das0p5) export DAS_ETA_KINEMATIC=0.5 ;;
      p2e_das1p0) export DAS_ETA_KINEMATIC=1.0 ;;
      p2e_das2p0) export DAS_ETA_KINEMATIC=2.0 ;;
      p2e_das4p0) export DAS_ETA_KINEMATIC=4.0 ;;
    esac ;;

  # (A2) Das arc construction A/B, 2026-07-28. Identical to p2e_das1p0 except
  # the offset follows the arc the TE actually sweeps instead of the tangent to
  # it (at eta=1.0: radius r and azimuth theta, vs r*sqrt(1+theta^2) and
  # atan(theta) -- a 0.015R radial / 1.1deg correction). Confirms the winning
  # eta=1.0 result (CT = 0.0713) under the correct construction.
  # SEPARATE CASE TAG ON PURPOSE: the launcher does `rm -rf data/$RUN_NAME`, so
  # re-running the p2e_das1p0 tag would destroy that result.
  p2e_das1p0_arc)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=1.0
    export DAS_KINEMATIC_ARC=true ;;

  # (D) Temporal convergence at the CONVERGED eta, 2026-07-28. The critical test
  # of the eta=1.0 result: hold eta=1.0 (semantics "shed one timestep back")
  # and halve dt, which halves Das in PHYSICAL units. At NT=36 one step of TE
  # travel at 0.75R is 1.02 chords, so the near wake is barely resolved.
  #   CT stays ~0.0713  => eta=1 is genuinely converged, 0.0713 is real.
  #   CT falls toward 0.0615 => the eta plateau was "near wake pushed out of
  #                             influence" and the +16% is a coarse-dt artifact.
  # rlxf halved to hold the physical relaxation rate as dt halves (same
  # convention as the p2e_*_nt72_d4 cases). ~2x the steps => needs --time=48:00:00.
  p2e_nt72_das1p0)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-72}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF="${RELAX_RLXF:-0.15}"
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=1.0 ;;

  # (A+B) Combination, 2026-07-28. The two levers that each worked on the
  # production mesh: eta=1.0 supplies the magnitude (0.0713, in band) and the
  # 0.5R filter supplies the lowest ripple measured here (1.64% / 4.34%). First
  # candidate to satisfy BOTH halves of BRAINSTORM/006's acceptance.
  # Tangent construction on purpose, so it is directly comparable to both
  # parents (p2e_das1p0 and p2e_filt0p5), which were tangent.
  p2e_das1p0_filt0p5)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=1.0
    export RELAX_FILTER_DOWNSTREAM_R=0.5 ;;

  # (B) Filter-depth survivability. Shallower filter = smaller unrelaxed
  # near-disk band; on 40_40 the depth curve was nearly flat (0.0672 at 0.5R vs
  # 0.0691 at 2.0R), so a shallow filter may survive here and still raise CT.
  p2e_filt0p5|p2e_filt1p0)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    case "$CASE" in
      p2e_filt0p5) export RELAX_FILTER_DOWNSTREAM_R=0.5 ;;
      p2e_filt1p0) export RELAX_FILTER_DOWNSTREAM_R=1.0 ;;
    esac ;;

  # (C) Unfiltered refined sigma -- separates "filter destabilises" from
  # "under-resolved wake cannot tolerate reduced damping". NO relaxation filter.
  p2e_sigF_nofilt)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export OVERLAP="${OVERLAP:-4}"; export P_PER_STEP=8; export MERGE_R_FACTOR="${MERGE_R_FACTOR:-0.0066}" ;;

  # (E) shed_with_induced_velocity null-test. Driver ships `false`, so after the
  # freestream is withdrawn the newly shed row gets essentially NO convection
  # before being shifted back a row; `true` lets it respond to self-induction.
  p2e_shedind)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export SHED_WITH_INDUCED_VELOCITY=true ;;

  # (F) Coupling ablations -- instability diagnostics, NOT production configs.
  # NOTE body->particle and wakerow->particle velocity GRADIENTS are already off
  # by default (BODY_HESSIAN_TO_PARTICLES / PANEL_WAKE_HESSIAN_TO_PARTICLES),
  # so the only live gradient is particle-self, i.e. vortex stretching.
  #   nostretch  : disable particle self-Hessian -> no vortex stretching
  #   nobodyonwake: remove body-induced VELOCITY on the wake
  p2e_abl_nostretch|p2e_abl_nobodyonwake)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export RELAX_FILTER_DOWNSTREAM_R="${RELAX_FILTER_DOWNSTREAM_R:-2.0}"
    case "$CASE" in
      p2e_abl_nostretch)    export PARTICLE_HESSIAN_SELF=false ;;
      p2e_abl_nobodyonwake) export BODY_ON_WAKE=false ;;
    esac ;;

  # ---- BRAINSTORM 014 (2026-07-28): first-wake-row offset attribution -------
  # All on the unfiltered carrier (velocity, NT=36, 4R, SETTLE_REVS=12, 719
  # steps, tangent Das) so each is like-for-like against the completed ladder:
  # eta=0.2 -> 0.06148, eta=1.0 -> 0.07133, eta=4.0 -> 0.07190.

  # (R1, HIGHEST PRIORITY) Floor test: the hard-coded min_displacement=0.01R
  # clamps Das over ~68% of the shedding span at eta=0.2, so the eta ladder was
  # not a clean single-parameter sweep. One change vs p2e_vel_nt36_d4: floor off.
  # Pre-registered interpretation (BRAINSTORM/014): CT -> 0.0713 => floor
  # artifact; CT stays ~0.0615 => eta is the real lever; between => both.
  p2e_das0p2_nofloor)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=0.2
    export DAS_MIN_DISPLACEMENT_R=1e-6 ;;

  # (R2a) Tiny rigid Das + one FREE convecting row: particle handoff ~ 1.0 step
  # of TE travel (between the completed eta=1.0 -> 0.4 and eta=4.0 -> 1.6
  # points, both ~0.0713-0.0719), but the rigid attached sheet is ~10x shorter
  # than eta=1.0's. CT ~ 0.0713 => only the total handoff distance matters and
  # the rigid sheet should be as short as the Kutta enforcement permits;
  # CT -> ~0.0615 => the rigid Das row itself carries the effect. Floor off so
  # the tiny Das is real and not replaced by 0.01R.
  p2e_nrows2_dassmall)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=0.1
    export DAS_MIN_DISPLACEMENT_R=1e-6
    export NWAKEROWS=2 ;;

  # (R2b) Fixed Das (eta=1.0), MORE free rows: handoff ~ 3.4 steps. Tests
  # insensitivity to panel-wake extent with the rigid sheet held fixed.
  p2e_nrows4_das1p0)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=1.0
    export NWAKEROWS=4 ;;

  # (R4, 2026-07-29) Sheet-extent asymptote hunt: nwakerows ladder UP at the
  # fixed eta=1.0 rigid offset. 2026-07-29 harvest: nrows 1 -> 4 raised CT
  # 0.0713 -> 0.0743 (+4.2%), so panel-sheet extent is its own lever; the
  # methodology is convergent only if CT(N) plateaus at some N*. 36 rows = one
  # full revolution of sheet. Unfiltered, stock relaxation (2026-07-28 ruling).
  p2e_nrows8_das1p0|p2e_nrows16_das1p0|p2e_nrows36_das1p0|p2e_nrows72_das1p0)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=1.0
    case "$CASE" in
      p2e_nrows8_das1p0)  export NWAKEROWS=8  ;;
      p2e_nrows16_das1p0) export NWAKEROWS=16 ;;
      p2e_nrows36_das1p0) export NWAKEROWS=36 ;;
      p2e_nrows72_das1p0) export NWAKEROWS=72 ;;
    esac ;;

  # (D-matched, 2026-07-29) Temporal convergence done PROPERLY (006 study D):
  # NT=72 with Das and sigma pinned at their NT=36 physical values — eta
  # scaled 1.0 -> 2.0 (frozen Das = 0.41 local chords again) and OVERLAP
  # scaled 3 -> 6 (sigma = dist*OVERLAP/PPS with dist halved). Isolates dt.
  # Context: naive p2e_nt72_das1p0 gave 0.07337 (+2.9% vs NT=36's 0.07133),
  # OUTSIDE the log-law pre-registration (0.0708-0.0712), but it halved Das,
  # sigma/chord AND rlxf at once, so it attributes nothing by itself.
  # Needs --time=48:00:00 on submission.
  p2e_nt72_das2p0_ov6)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-72}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF="${RELAX_RLXF:-0.15}"
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=2.0
    export OVERLAP=6 ;;

  # (R3) Un-frozen Das: re-derived from the CURRENT kinematics every step, so
  # the 0.4x spin-up factor disappears (eta_eff = eta = 1.0, i.e. Das = one full
  # step of TE travel at operating RPM -- 2.5x the frozen eta=1.0 length).
  # Tests whether 0.0713 survives making eta physically meaningful.
  p2e_das1p0_refresh)
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"; export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R=4; export RELAX_RLXF=0.3
    export SETTLE_REVS="${SETTLE_REVS:-12}"; export CONVERGENCE_REVS="${CONVERGENCE_REVS:-10}"
    export DAS_ETA_KINEMATIC=1.0
    export DAS_REFRESH=true ;;

  # ---- BRAINSTORM 018: publishable hover convergence campaign ---------------
  # See BRAINSTORM/018_dji9443_hover_convergence_campaign.md (+ subdir) for the
  # phase plan, decision rules, and the restart-chaining recipe. Rulings baked
  # in here (Ryan 2026-07-30): sigma_overlap shedding; stock relaxation with NO
  # near-disk filter anywhere; DAS_REFRESH never true; NWAKEROWS=4 by the
  # handoff-distance criterion (d >= 4*sigma); flat 64G.
  #
  # sigma = 2*pi*R/NT * OVERLAP/PPS (driver line ~336). At NT=36:
  # sigma/c(0.75R) = 1.363*OVERLAP/PPS. To pin sigma when NT doubles, HALVE
  # P_PER_STEP (the SigmaOverlap policy sheds the correct per-step count on its
  # own). MERGE_R_FACTOR = 0.138*(sigma/R) holds r_merge/sigma at stock ratio.
  #
  # Every knob below honors an environment override (sbatch --export=ALL,...),
  # which is how the final-settings cases (trunc6/green/nomerge/final) receive
  # Das*/NT*/sigma* once earlier phases pick them, and how restart segments
  # extend a run (P018_RUN_NAME / P018_SETTLE_REVS / RESTART_*).
  p018_*)
    unset RELAX_FILTER_DOWNSTREAM_R   # ruling: filter OFF in every 018 run
    export RHPC_FORMULATION="${RHPC_FORMULATION:-velocity}"
    export NT="${NT:-36}"
    export TRUNCATION_DEPTH_R="${TRUNCATION_DEPTH_R:-4}"
    export RELAX_RLXF="${RELAX_RLXF:-0.3}"
    export SETTLE_REVS="${P018_SETTLE_REVS:-12}"
    export CONVERGENCE_REVS=10
    export PARTICLE_SHEDDING=sigma_overlap
    export DAS_ETA_KINEMATIC="${DAS_ETA_KINEMATIC:-1.0}"
    export NWAKEROWS="${NWAKEROWS:-4}"
    export OVERLAP="${OVERLAP:-2.0}"
    export P_PER_STEP="${P_PER_STEP:-4}"
    export MERGE_R_FACTOR="${MERGE_R_FACTOR:-0.0120}"
    case "$CASE" in
      # Phase 1 baseline: eta=1.0, NT=36, sigma/c=0.68 (OVERLAP 2.0 / PPS 4).
      p018_b0) : ;;
      # Phase 2: physical-Das ladder at the B0 carrier (B0 supplies eta=1.0).
      p018_das0p5) export DAS_ETA_KINEMATIC=0.5 ;;
      p018_das2p0) export DAS_ETA_KINEMATIC=2.0 ;;
      p018_das4p0) export DAS_ETA_KINEMATIC=4.0 ;;
      # Floor-contamination bound for the smallest ladder point (optional).
      p018_das0p5_floorhalf) export DAS_ETA_KINEMATIC=0.5; export DAS_MIN_DISPLACEMENT_R=0.005 ;;
      # Phase 3: NT ladder with Das AND sigma pinned in PHYSICAL units
      # (eta prop. to NT, PPS inversely prop. to NT, rlxf halved per doubling).
      # Defaults assume Das* = 0.41c (eta*=1.0); override DAS_ETA_KINEMATIC at
      # submission if Phase 2 chose differently. nt72 needs --time=48:00:00;
      # nt144 runs as restart-chained 48 h segments.
      p018_nt72)  export NT=72;  export DAS_ETA_KINEMATIC="${DAS_ETA_KINEMATIC:-2.0}"; export P_PER_STEP=2; export RELAX_RLXF=0.15  ;;
      p018_nt144) export NT=144; export DAS_ETA_KINEMATIC="${DAS_ETA_KINEMATIC:-4.0}"; export P_PER_STEP=1; export RELAX_RLXF=0.075 ;;
      # Phase 4: sigma ladder with h refined FASTER than sigma (q ~ 1.2:
      # OVERLAP rises down the ladder). L1_warm is L1 warm-started from B0
      # (pass RESTART_* at submission). L1 ~45 h, L2 ~90 h: restart-chain.
      p018_L1|p018_L1_warm) export OVERLAP=2.4;  export P_PER_STEP=11; export MERGE_R_FACTOR=0.0053 ;;
      p018_L2)              export OVERLAP=2.88; export P_PER_STEP=26; export MERGE_R_FACTOR=0.0027 ;;
      # Phase 5: single N=8 spot-check of the handoff criterion (no ladder --
      # the legacy-sigma nwakerows ladder was non-monotone and rejected).
      p018_nrows8) export NWAKEROWS=8 ;;
      # Phases 6-9: final-settings carriers -- pass Das*/NT*/OVERLAP/PPS/
      # MERGE_R_FACTOR via --export at submission.
      p018_trunc6)       export TRUNCATION_DEPTH_R=6 ;;
      p018_green)        export RHPC_FORMULATION=green ;;
      p018_green_coarse) export RHPC_FORMULATION=green ;;
      p018_nomerge)      export MERGE_PARTICLES=false ;;
      # Phases 10/11: mesh refinement at final settings (pass Das*/NT*/sigma*
      # via --export). Spanwise rungs 60/80 must be generated locally (Mac
      # OpenVSP, 2d cap recipe: flat root + round tip captess4) and scp'd to
      # examples/data/ BEFORE submission. All cold starts (no cross-mesh warm
      # start). chord249/span60 dense G ~20 GB; span80 ~34 GB (consider
      # --mem=128G on submission).
      p018_chord145) export RHPC_MESH=45_145_ct4 ;;
      p018_chord249) export RHPC_MESH_FILE=dji9443_20260725_45_249_capped_captess4.msh ;;
      p018_span60)   export RHPC_MESH_FILE=dji9443_20260731_60_185_capped_captess4.msh ;;
      p018_span80)   export RHPC_MESH_FILE=dji9443_20260731_80_185_capped_captess4.msh ;;
      p018_final) : ;;
      *) echo "ERROR: unknown p018 case tag '$CASE'" >&2; exit 2 ;;
    esac ;;

  *) echo "ERROR: unknown case tag '$CASE'" >&2; exit 2 ;;
esac

export RUN_NAME="$CASE"
# 018 restart segments write to their own directory (P018_RUN_NAME=<case>_s2
# etc.) so the wipe below cannot destroy the restart source. See
# BRAINSTORM/018_.../ops_reference.md for the chaining recipe.
[[ -n "${P018_RUN_NAME:-}" ]] && export RUN_NAME="$P018_RUN_NAME"

# Clear the previous attempt's run directory so a rerun cannot leave stale steps
# past the current run length (HPC policy). Destructive by default: Ryan
# (2026-07-27) "in general, you can wipe old data. We don't have space to keep
# all paraview files."
#
# Opt in with RHPC_KEEP_PREV=true to move the old run to data/<name>.prev
# instead. Worth doing when the previous attempt FAILED and has not been
# inspected yet -- a plain wipe destroyed three OOM'd Batch-2 histories on
# 2026-07-27 minutes before they were wanted for ParaView. Harvest the CT series
# from logs/slurm/*.out first and the VTK is usually expendable.
if [[ -d "data/$RUN_NAME" ]]; then
  if [[ "${RHPC_KEEP_PREV:-false}" == "true" ]]; then
    rm -rf "data/$RUN_NAME.prev"
    mv "data/$RUN_NAME" "data/$RUN_NAME.prev"
    echo "  previous run preserved at data/$RUN_NAME.prev"
  else
    rm -rf "data/$RUN_NAME"
  fi
fi

echo "DJI 9443 Phase 2e unsteady hover CT — case $CASE"
echo "  repo:$PWD threads:$THREADS host:$(hostname) job:${SLURM_JOB_ID:-none}"
echo "  mesh:$RHPC_MESH formulation:$RHPC_FORMULATION RPM:$RPM NT:$NT depth:${TRUNCATION_DEPTH_R}R rlxf:$RELAX_RLXF"
echo "  relax_filter:${RELAX_FILTER_DOWNSTREAM_R:-off}R das_eta:${DAS_ETA_KINEMATIC:-0.2} overlap:${OVERLAP:-3} pps:${P_PER_STEP:-2} merge_r:${MERGE_R_FACTOR:-0.02} settle:${SETTLE_REVS} das_arc:${DAS_KINEMATIC_ARC}"
echo "  das_min_R:${DAS_MIN_DISPLACEMENT_R:-0.01} das_refresh:${DAS_REFRESH:-false} nwakerows:${NWAKEROWS:-1}"
echo "  started $(date '+%F %T')"

julia --project=. -t "$THREADS" examples/rotor_hover_pressure_comparison.jl

echo "Case $CASE finished ($(date '+%F %T'))."
echo "Artifacts: data/$RUN_NAME/${RUN_NAME}_CT_vs_rev.csv"
echo "           data/$RUN_NAME/${RUN_NAME}_CT_per_rev.csv"
echo "           data/$RUN_NAME/${RUN_NAME}_case_metadata.toml"

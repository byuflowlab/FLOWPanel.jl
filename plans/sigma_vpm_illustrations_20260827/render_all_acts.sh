set -u
# Session-5 (2026-08-29) act structure: regime-based three-act story.
# Acts I-II = the s025/s025v pair (same sigma_shed, inviscid vs viscous).
# HUD rev convention: rev = (step - spinup)/nt, i.e. spinup rev excluded.
cd /Users/ryan/Dropbox/research/projects/FLOWPanel.jl
S=/private/tmp/claude-502/-Users-ryan-Dropbox-research-projects-FLOWPanel-jl/2c6862dd-d8a5-430e-a379-84e453dda661/scratchpad
C=plans/sigma_vpm_illustrations_20260827/wake_health_csvs
D=$HOME/Dropbox/research/notebooks/img/20260827_sigma_vpm_illustrations
R=/Users/ryan/Dropbox/research/projects/FLOWPanel.jl
COMMON="--color-by abssigma --yaw 0 --sigma-log --sigma-clim 0.1:10 --sigma-clim-units shed --arrow-scale 1.0 --arrow-cap 0 --rotor-radius 0.11938 --font-size 18 --nt 36 --spinup-steps 36"
# Ryan 2026-08-29: Act-I gifs use FULLY OPAQUE arrows and a matched 266:288 window
ACTI_EXTRA="--arrow-opacity-min 1.0"

run () { # name glob steps sigma_shed hud title [extra...]
  n=$1; g=$2; s=$3; sh=$4; h=$5; ti=$6; shift 6
  echo "=== $n"
  python3 scripts/illus_render_particle_series.py --glob "$g" --steps "$s" \
    --sigma-shed "$sh" --hud-csv "$h" $COMMON "$@" --out "$D/$n" --title "$ti" 2>&1 | tail -1
}

# actI_blowup: REVERTED to the approved 266:305 faded-arrow version (Ryan 2026-08-29)
run actI_blowup "$S/vtp_s025/*.vtp" 266:305 2.967e-3 \
  "$C/scr_p019_s025_rr_monitor04_wake_health_system1.csv" \
  'Act I blow-up: inviscid sigma/R=0.0249 - strain thinning unopposed, sigma through zero'

# opaque-arrow 266:288 twin of the blow-up, used ONLY as the composite's left pane
run actI_blowup_opq "$S/vtp_s025/*.vtp" 266:288 2.967e-3 \
  "$C/scr_p019_s025_rr_monitor04_wake_health_system1.csv" \
  'Act I blow-up: inviscid sigma/R=0.0249 - strain thinning unopposed, sigma through zero' \
  $ACTI_EXTRA

run actI_fixed "$S/vtp_s025v/*.vtp" 266:288 2.967e-3 \
  "$C/scr_p019_s025v_rr_monitor04_wake_health_system1.csv" \
  'Act I fixed: same sigma/R=0.0249, viscosity ON - sigma arrested, wake healthy' \
  $ACTI_EXTRA
# milestone overlay REMOVED (Ryan 2026-08-29): the side-by-side composite is the
# Act-I "fixed" presentation; single actI_fixed.gif is a building block.

run actII_runaway "$S/vtp_s025v/*.vtp" 285:323 2.967e-3 \
  "$C/scr_p019_s025v_rr_monitor04_wake_health_system1.csv" \
  'Act II: viscosity helped - but Gamma runaway at the equilibrium (Dt Z -> O(1))'

run actIII_drift_death "$S/vtp_ufront_s035/*.vtp" 987:996 4.236e-3 \
  "$R/data/p018_ufront_s035_visc/monitors/p018_ufront_s035_visc_monitor04_wake_health_system1.csv" \
  'Act III drift death: viscous sigma/R=0.0356 - slow sigma drift to the viscous floor'

# actIII_fixed_matched / actIII_fixed_steady: from the task-4 long sigma* run
# (scr_p019_s038v_gpu40) once it lands -- add here with its sigma_shed and CSV.

# Composite Act-I side-by-side (inviscid | viscous), matched frames:
mkdir -p "$D/actI_sbs"; rm -f "$D/actI_sbs"/frame_*.png
for s in $(seq 266 288); do magick "$D/actI_blowup_opq/frame_00$s.png" "$D/actI_fixed/frame_00$s.png" +append "$D/actI_sbs/frame_00$s.png"; done
magick -delay 25 "$D/actI_sbs"/frame_*.png -loop 0 "$D/actI_sidebyside.gif"

# RETIRED from the act set 2026-08-29 (files may remain on disk):
#   actI_fixed (OLD s030v version -- s030v ignites rev 8.47 on the rerun stack),
#   actII_blowup (p018_L2 forensics corpse), actII_fixed_geom (expint; F7 coda
#   TikZ only), actIII_fast_ignition (renamed actII_runaway).
echo "=== ALL DONE"

#!/usr/bin/env julia

# Data-only validation notes for the rotor-hover pressure comparison runs.
#
# Source used for extraction:
#   data/rotor_hover_pressure_comparison_nt72_sfs
# via:
#   env RUN_NAME=rotor_hover_pressure_comparison_nt72_sfs RUN_STAGE1_LOG=true \
#     STAGE1_LOG_PATH=/private/tmp/rotor_nt72_stage1.csv \
#     julia --project=. examples/rotor_hover_pressure_comparison_replay.jl
#
# The saved nt72_sfs run contains 665 frame samples and ends at
# t = 0.09222222222222222 s (step 664). The comparable nt36 run contains
# 720 frame samples and ends at t = 0.09986111111111111 s (step 719).
#
# Measured on the last three saved samples of nt72_sfs:
# - Kutta-Joukowski CT: -0.050964309, -0.050475825, -0.049400557
# - KJ mean CT: -0.05028023
# - Bernoulli CT: -0.13667324, 0.00029200244, 0.00031150275
# - Laplace MD CT: -5.0069342e-5, -5.0069342e-5, -5.0069342e-5
# - Laplace lamb CT: -0.0029944656, -0.0029947082, -0.0029934218
#
# Stage-1 pressure diagnostics from the same replayed samples:
# - Bernoulli cancellation ratio: 1.8169 -> 735.637 -> 689.298
# - Laplace MD cancellation ratio: 3788.996 -> 3788.996 -> 3788.996
# - Laplace lamb cancellation ratio: 4.1212 -> 4.0982 -> 4.1040
# - Pressure to panel-force correlation stays weak:
#   Bernoulli ~ 0, Laplace MD ~ -0.0029, Laplace lamb ~ 0.061
#
# Suggestions for the discrepancy versus nt36:
# 1. The nt72_sfs dataset stops about 0.00764 s earlier than nt36, so the
#    comparison is not at the same late-time azimuth. That alone can change
#    pressure recovery and wake-induced pressure balance.
# 2. The KJ estimate is stable near -0.05, but the pressure-derived CTs are
#    either near zero (Laplace MD), much smaller than KJ (Laplace lamb), or
#    strongly step-sensitive (Bernoulli). That points to pressure reconstruction
#    / cancellation sensitivity rather than an unstable KJ estimate.
# 3. The pressure-panel correlations are weak while the cancellation ratios are
#    large, which is consistent with a pressure recovery issue, a Kutta-correction
#    mismatch, or an incompletely converged wake state.
# 4. The nt72_sfs name suggests an SFS variant, but the saved replay metadata
#    does not encode the viscous/SFS provenance. If the generating script changed
#    any wake model detail between nt36 and nt72_sfs, that is a plausible source
#    of the difference.

const ROTOR_HOVER_VALIDATE_NOTE = raw"""
Rotor hover validation note
===========================

nt72_sfs replay summary
-----------------------
Saved steps: 665
Final saved step: 664
Final saved time: 0.09222222222222222 s

Kutta-Joukowski CT on the last three replayed samples:
  step 662: -0.050964309
  step 663: -0.050475825
  step 664: -0.049400557
  mean:     -0.05028023

Pressure-recovered CT on the same samples:
  Bernoulli:   -0.13667324, 0.00029200244, 0.00031150275
  Laplace MD:  -5.0069342e-5, -5.0069342e-5, -5.0069342e-5
  Laplace lamb:-0.0029944656, -0.0029947082, -0.0029934218

Why nt36 and nt72_sfs may disagree
----------------------------------
The saved nt72_sfs run ends earlier than the saved nt36 run, so the two
datasets are not compared at the same physical time or azimuth.

Within nt72_sfs itself, the KJ force is stable near -0.05, while the
pressure-derived monitors do not agree with one another or with KJ. The
pressure balance is highly cancellation-sensitive, which makes it vulnerable
to small differences in wake state, pressure reconstruction, and Kutta handling.
"""

if abspath(PROGRAM_FILE) == @__FILE__
    println(ROTOR_HOVER_VALIDATE_NOTE)
end

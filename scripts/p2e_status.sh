#!/usr/bin/env bash
#
# Phase 2e status probe, run ON the cluster. Exists because nested ssh quoting
# kept corrupting remote readings three separate ways: the login MOTD banner
# polluting `tr -dc '0-9'` scalar extraction, `tail -N` clipping the squeue
# header and job rows, and $(...) expansion collapsing inside a shell function.
# Keep all parsing here, and have the caller read sentinel-prefixed lines only.
#
# Usage (from the local machine):
#   ssh orc /home/rander39/projects/FLOWPanel.jl/scripts/p2e_status.sh
#   ssh orc /home/rander39/projects/FLOWPanel.jl/scripts/p2e_status.sh brief
#
# Sentinel lines (safe to grep): P2E_NJOBS=, P2E_NBAD=, P2E_JOB=, P2E_STEP=, P2E_CT=

set -uo pipefail
export PATH=/apps/slurm/latest/bin:$PATH
cd /home/rander39/projects/FLOWPanel.jl || exit 2

MODE="${1:-full}"

njobs=$(squeue -u rander39 --noheader 2>/dev/null | grep -c fp-p2e)
echo "P2E_NJOBS=${njobs}"

# Divergence detection. NaN/Inf alone is NOT sufficient: job 12894481 diverged to
# CT ~ 90 in large-but-finite arithmetic and slipped straight past a NaN/Inf
# check, only surfacing as an OutOfMemoryError inside merge_particles! once the
# blown-up particle cloud wrecked the spatial hash. Flag on magnitude too:
# hover CT is O(0.06), so |CT| > 1 is unambiguously diverged.
nbad=0
for f in logs/slurm/*p2e*.out; do
    [ -f "$f" ] || continue
    if grep -qE 'NaN|Inf' "$f" 2>/dev/null; then
        nbad=$((nbad+1)); echo "P2E_DIVERGED=$(basename "$f") reason=nan_inf"; continue
    fi
    mx=$(grep -o 'CF = (-\?[0-9.eE+-]*' "$f" 2>/dev/null | sed 's/CF = (//' \
         | awk 'BEGIN{m=0}{v=$1<0?-$1:$1; if(v>m)m=v}END{printf "%.4g", m}')
    case "$mx" in ''|0) continue;; esac
    if awk "BEGIN{exit !($mx > 1.0)}"; then
        nbad=$((nbad+1)); echo "P2E_DIVERGED=$(basename "$f") reason=magnitude max_abs_CT=${mx}"
    fi
done
echo "P2E_NBAD=${nbad}"

squeue -u rander39 --noheader -o "%i %j %T %M" 2>/dev/null | grep fp-p2e | while read -r id name state elapsed; do
    echo "P2E_JOB=${id} ${name} ${state} ${elapsed}"
done

[ "$MODE" = brief ] && exit 0

for d in data/p2e_green_f2p0_* data/p2e_vel_nt36_d4 data/p2e_green_nt36_d4; do
    [ -d "$d" ] || continue
    b=$(basename "$d")
    hi=$(ls "$d/${b}_body1"/*.vtu 2>/dev/null | sed 's/.*\.\([0-9]*\)\.vtu/\1/' | sort -n | tail -1)
    np=$(ls "$d/${b}_wake1_particles"/*.vtp 2>/dev/null | wc -l | tr -d ' ')
    echo "P2E_STEP=${b} step=${hi:-none} nvtp=${np}"
done

# Per-revolution CT (thrust convention) from the ForceMonitor stream in each log.
for f in logs/slurm/*p2e*.out; do
    [ -f "$f" ] || continue
    case "$f" in *nt72*) nt=72;; *) nt=36;; esac
    tag=$(basename "$f" | sed 's/slurm-fp-p2e-//; s/-[0-9]*\.out//')
    awk -v nt="$nt" -v tag="$tag" '
        /ForceMonitor\[i_system=1, step=/ {
            l=$0; sub(/.*step=/,"",l); sub(/\].*/,"",l); s=l+0; want=1; next
        }
        want && /CF = / {
            t=$0; sub(/.*CF = \(/,"",t); sub(/,.*/,"",t); ct=-(t+0)
            r=int((s-1)/nt); sum[r]+=ct; n[r]++; want=0
        }
        END {
            m=0; for (r=0; r<400; r++) if (n[r]==nt) a[m++]=sprintf("r%d=%.5f", r, sum[r]/n[r])
            out=""; for (i=(m>6?m-6:0); i<m; i++) out=out" "a[i]
            printf "P2E_CT=%s nrev=%d%s\n", tag, m, out
        }
    ' "$f"
done

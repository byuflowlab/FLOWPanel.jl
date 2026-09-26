# Warm-start R4 harvest: fgs-wsr4-13890195

Windows (transients INCLUDED, Ryan 2026-09-24): A = steps 1-36 (from the very first step), B = the fourth revolution (steps 109-144, restarted from the family checkpoint at step 108 where restart legs exist).

**Reporting note (Ryan 2026-09-24):** solver warm-start histories are not serialized in the restart checkpoint, so the first (order+1) steps of each restarted WARM leg are effectively cold — that history-fill transient is INSIDE Window B's transient-included statistics by design (never excluded). Interpret the first few Window-B steps of warm arms accordingly; the per-step traces make the refill visible.

## Window A

| arm | t_solve [s] | niter_first | t_project [s] | setup t_setup / t_prime [s] | unconverged | nsolves≠1 | bcerr>tol | uncertified |
|---|---|---|---|---|---|---|---|---|
| fgs_cold | 14.62 ± 6.7 (med 14.44) | 28.25 ± 2 (med 28) | 0 ± 0 (med 0) | 315.351271 / 0 | 0 | 0 | 0 | 0 |
| fgs_prev | 12.58 ± 9.3 (med 12.4) | 19.89 ± 14 (med 20) | 0.001012 ± 0.0012 (med 0.001032) | 317.513279 / 0 | 0 | 0 | 0 | 0 |
| fgs_proj1 | 11.7 ± 10 (med 11.35) | 16.25 ± 19 (med 15) | 0.001052 ± 0.0013 (med 0.001114) | 317.137698 / 0 | 0 | 0 | 0 | 0 |
| fgs_proj2 | 11.15 ± 11 (med 10.68) | 13.94 ± 23 (med 12) | 0.001111 ± 0.0014 (med 0.001177) | 319.512115 / 0 | 0 | 0 | 0 | 0 |
| ilu_nfcache_cold | 11.1 ± 2.6 (med 11.06) | 16.25 ± 1 (med 16) | 0 ± 0 (med 0) | 81.6707541 / 28.2211751 | 0 | 0 | 0 | 0 |
| ilu_nfcache_prev | 10.22 ± 2.8 (med 10.19) | 11.78 ± 6 (med 12) | 0 ± 0 (med 0) | 79.4091126 / 28.2381778 | 0 | 0 | 0 | 0 |
| ilu_nfcache_proj1 | 9.918 ± 3.6 (med 9.895) | 10.31 ± 7 (med 10) | 0.0003606 ± 0.00056 (med 0.0003318) | 81.8585148 / 28.1256612 | 0 | 0 | 0 | 0 |
| ilu_nfcache_proj2 | 9.121 ± 3 (med 9.114) | 9.639 ± 8 (med 9) | 0.0003472 ± 0.00057 (med 0.0003438) | 82.238598 / 27.6698095 | 0 | 0 | 0 | 0 |

## Window B

| arm | t_solve [s] | niter_first | t_project [s] | setup t_setup / t_prime [s] | unconverged | nsolves≠1 | bcerr>tol | uncertified |
|---|---|---|---|---|---|---|---|---|
| fgs_cold | 12.97 ± 5.1 (med 12.88) | 21.89 ± 1 (med 22) | 0 ± 0 (med 0) | 315.351271 / 0 | 0 | 0 | 1 | 0 |
| fgs_prev | 10.08 ± 8.7 (med 9.861) | 9.5 ± 14 (med 9) | 0.001001 ± 0.0012 (med 0.001008) | 317.513279 / 0 | 0 | 0 | 0 | 0 |
| fgs_proj1 | 8.854 ± 9.2 (med 8.503) | 4.194 ± 19 (med 3) | 0.001135 ± 0.0015 (med 0.001173) | 317.137698 / 0 | 0 | 0 | 0 | 0 |
| fgs_proj2 | 8.633 ± 9.7 (med 8.245) | 3.222 ± 20 (med 2) | 0.001163 ± 0.0016 (med 0.001238) | 319.512115 / 0 | 0 | 0 | 0 | 0 |
| ilu_nfcache_cold | 10.64 ± 1.3 (med 10.65) | 15 ± 0 (med 15) | 0 ± 0 (med 0) | 81.6707541 / 28.2211751 | 0 | 0 | 0 | 0 |
| ilu_nfcache_prev | 9.268 ± 3.3 (med 9.244) | 7.083 ± 9 (med 7) | 0 ± 0 (med 0) | 79.4091126 / 28.2381778 | 0 | 0 | 0 | 0 |
| ilu_nfcache_proj1 | 8.732 ± 3.4 (med 8.68) | 5.111 ± 11 (med 5) | 0.0003381 ± 0.00051 (med 0.0003102) | 81.8585148 / 28.1256612 | 0 | 0 | 0 | 0 |
| ilu_nfcache_proj2 | 8.25 ± 3.1 (med 8.215) | 4.417 ± 12 (med 4) | 0.0003192 ± 0.00052 (med 0.0003161) | 82.238598 / 27.6698095 | 0 | 0 | 0 | 0 |

## Cross-arm solution agreement (rel-L2 vs fgs_cold, per step)

| arm | winA mean | winA max | winB mean | winB max |
|---|---|---|---|---|
| fgs_prev | 3.590e-05 | 1.979e-04 | 1.023e-05 | 1.501e-05 |
| fgs_proj1 | 8.446e-05 | 2.391e-04 | 4.379e-05 | 6.280e-05 |
| fgs_proj2 | 4.552e-05 | 1.438e-04 | 3.054e-05 | 6.874e-05 |
| ilu_nfcache_cold | 6.855e-05 | 2.413e-04 | 5.956e-04 | 1.132e-03 |
| ilu_nfcache_prev | 7.374e-05 | 3.960e-04 | 5.985e-04 | 1.140e-03 |
| ilu_nfcache_proj1 | 6.399e-05 | 1.932e-04 | 5.977e-04 | 1.137e-03 |
| ilu_nfcache_proj2 | 6.719e-05 | 3.831e-04 | 5.981e-04 | 1.140e-03 |

Known context: FGS and Krylov converge to slightly different wake-on fixed points (~2e-3, rigid_motion_tree_reuse_item.md §5) — reported, not chased. Window-B deltas within a solver family share the family checkpoint, so they isolate the initial guess; ilu-vs-fgs window-B deltas ALSO carry the two checkpoints' divergence.

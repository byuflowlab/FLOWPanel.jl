# v15 j1-b1 bounded harvest (job 13694724)

Source: `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v15-13694724`.
Harvested root regular files and completed `j1-b1/`; fixture directories and
`j4-b1/` were intentionally excluded. The copied set is 43 regular files,
totaling 1,643,851 bytes for `j1-b1/` plus 128,214 bytes of root evidence.

| Trial group | n | median solve seconds | accepted | evaluator | BC rel-L2 |
|---|---:|---:|---:|---|---:|
| uninstrumented, batches 1+3 | 20 | 37.5594 | 20/20 | certified_fmm | 4.77952513e-7 |
| instrumented, batches 2+4 | 20 | 37.4689 | 20/20 | certified_fmm | 4.77952513e-7 |
| batch 1, uninstrumented | 10 | 37.5386 | 10/10 | certified_fmm | 4.77952513e-7 |
| batch 2, instrumented | 10 | 37.5275 | 10/10 | certified_fmm | 4.77952513e-7 |
| batch 3, uninstrumented | 10 | 37.5687 | 10/10 | certified_fmm | 4.77952513e-7 |
| batch 4, instrumented | 10 | 37.4598 | 10/10 | certified_fmm | 4.77952513e-7 |

All rows report 27 iterations, 81 estimated inner sweeps, 28 estimated FMM
passes, zero relative solution delta, and direct/FMM delta as NaN because the
authoritative evaluator was certified FMM. Instrumented stage medians are
diagnostic only; do not use them for speed claims.

The capability probe on m12-2-26 succeeded with `/usr/bin/perf`:
`perf_probe_exit=0` and `perf_userspace_probe_exit=0`. It successfully
read cycles, instructions, cache references, and cache misses on the probe
command. This is only PMU capability evidence; no workload-scoped bandwidth
or cache measurement has been performed.

Exact per-file remote SHA256 values are in `remote-sha256.txt`; all 43
local hashes matched, as recorded in `sha256-verification.txt`.

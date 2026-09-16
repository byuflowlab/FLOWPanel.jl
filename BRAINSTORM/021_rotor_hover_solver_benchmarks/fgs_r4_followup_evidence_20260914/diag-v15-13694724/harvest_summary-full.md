# v15 harvest — job 13694724

Harvested 2026-09-16 UTC from `/home/rander39/projects/FLOWPanel.jl/data/p021-cold-20260910/diag-v15-13694724/` over the existing SSH ControlMaster. No remote files were modified; job was not resubmitted.

| Arm | Result files | Rows | Uninstrumented median (s) | Instrumented median (s) | Audit |
|---|---:|---:|---:|---:|---|
| j1/b1 | 18 | 40 | 37.5594 | 37.4689 | PASS |
| j4/b1 | 18 | 40 | 16.9358 | 17.0866 | PASS |
| j8/b1 | 18 | 40 | 13.4642 | 13.3843 | PASS |
| j16/b1 | 18 | 40 | 12.1994 | 11.9373 | PASS |
| j32/b1 | 18 | 40 | 11.3673 | 11.4110 | PASS |
| j64/b1 | 18 | 40 | 10.9626 | 10.9483 | PASS |

Root `COMPLETED` and all six arm `results/status.toml` markers are present. The scheduler `.out` and `.err` files were harvested and are zero bytes.

All six arm provenance files report Julia 1.11.7, BLAS threads 1, stage `verify`, and the same clean campaign pins: FLOWPanel `campaign/p021-cold-source-20260915-v15` (`39ec4e3630bc6c04f0865a1d3feecce7130904fe`), FastMultipole `campaign/p021-r4-diag-source-20260914-v10` (`87cbc8460b51f24ddf34cc5f41a1d1b6682bf04a`), and FLOWVPM `campaign/p021-cold-exec-20260910-v1` (`05c658f7804ec5f9b68d4cb9826a9f97cfecb373`). Arm Julia threads are 1/4/8/16/32/64 respectively.

`remote-sha256-full.txt` contains 131 remote digests: 23 root files, 16 `results/` files per arm (96 total), and two arm-level files per arm (`process.log` and `numactl_show.txt`, 12 total). `sha256-verification-full.txt` reports 131/131 `OK` and no failures. The two separately harvested scheduler logs also match remote SHA256 (`e3b0c442...2b855`, empty files). `analysis/audit-full.txt` is the prepared numerical/completeness audit and exits 0. The audit is data-only; provenance is independently checked above.

Cleanup provenance is preserved in `silo-provenance/`: remote-to-local SHA256 matches for `pins.toml`, `env/Project.toml`, `env/Manifest.toml`, `p021-r4-diag-flowpanel-v15.sha256`, and `p021-r4-diag-fastmultipole-v10.sha256`. The remote silo inventory at harvest consisted of `FLOWPanel.jl/`, `FastMultipole/`, `counters-v16/`, `env/`, the v15 FLOWPanel/FastMultipole content manifests, and `pins.toml`; no deletion was performed.

# v17 Julia FIFO protocol verification

Date: 2026-09-15 local MDT.

The v16 `counter_command` issued a FIFO command and then used an unbounded
`readline` for its perf acknowledgement. The Python smoke had a 15-second
`select` bound, but no Julia command had that bound. An initial asynchronous
Julia-read prototype was rejected because a timed-out read task remained
blocked.

The v17 helper uses POSIX `poll(2)` followed by raw one-byte `read(2)` calls.
Each byte of the required `ack\\n` shares one monotonic 15-second deadline.
This prevents both missing acknowledgements and partial acknowledgements from
blocking the Julia process. It bypasses Julia stream buffering after readiness
polls, so buffered prefetch cannot cause a later poll to wait for data Julia
already owns.

The narrow local driver ran with two Julia threads and passed. Its real-FIFO
checks cover valid `ack\\n`, missing ACK, malformed `nack\\n`, and a partial
`a`; the last two fail and both timeout paths complete in under one second.
`bash -n benchmark/run_r4_counters.slurm.sh` and `git diff --check` also pass.

The launcher now labels paths and Slurm artifacts v17. Before fixtures it runs
the extracted Julia helper under installed `perf stat --control`, with one
Julia thread, BLAS thread limits, and `--startup-file=no`; this remains the
required Linux/cluster gate. It has not run locally because macOS has neither
the Linux perf interface nor `/proc` activity support. No full FLOWPanel solver
control or cluster operation was repeated in this bounded verification.

Pin status: v17 worktree is modified and uncommitted; no tag, deployment, or
submission has occurred.

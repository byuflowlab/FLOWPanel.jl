# What a different triangular solver could change about FGS scaling

Date: 2026-09-26.

Follow-up to [the cooperative leaf scaling plan](fgs_cooperative_leaf_scaling_plan_20260926.md).

Ryan's question:

> "Preserve the GS recurrence through a different triangular solver"; I worry that for larger problems, the number of leaves increases, and we don't change the amount of work that can be done in parallel. E.g., larger problems introduce more work, but our rate of performing that work remains the same. What do you anticipate coming out of a different triangular solver?

**Your concern is valid. A different scheduler for the same leaf dependencies would not solve it.** If increasing problem size lengthens the dependency chain without widening it, more leaves simply mean more work at roughly the same limited throughput.

What I meant by a different triangular solver is an **algebraic transformation that changes the work and its dependencies**, while producing the same GS update. That possibility is real, but my recommendation should have been more explicitly conditional.

A GS sweep solves

$$
(D+L)x^{k+1}=b-Ux^k.
$$

The current implementation evaluates this by forward substitution through the leaf DAG. Preserving the solution of that triangular system does **not** require preserving that particular sequence of operations.

For example, consider a simple chain:

$$
x_i=a_i x_{i-1}+c_i.
$$

Ordinary substitution has a dependency chain of length $N$. But each segment can independently construct a relation between its input and output. Those relations can then be composed in a parallel tree, giving logarithmic dependency depth in this scalar example. The result is mathematically the same. For general block systems, however, those relations become matrices and can become expensive to construct and store.

The proposed subdomain/interface approach applies a related idea:

1. Precompute how each subdomain responds to its incoming interface values.
2. For each new RHS, compute independent subdomain particular solutions.
3. Solve a reduced system for the interface values.
4. Complete the subdomain solutions independently.

The SPIKE-based method I cited does this for triangular systems. It explicitly introduces fill and retains a reduced triangular solve as a potential sequential bottleneck. [Torun et al.](https://user.ceng.metu.edu.tr/~manguoglu/PDFs/dmpGS_SISC.pdf)

**What I would hope to gain is a change from one long global chain to many shorter local chains plus a substantially smaller interface problem.** With 64 subdomains, an illustrative timing model is

$$
T_{64}(N)\approx T_{\mathrm{local}}(N/64)
+T_{\mathrm{interface}}(N)
+T_{\mathrm{assembly/recovery}}(N,64).
$$

This helps only if the interface solve stays inexpensive and the additional work is manageable. Your objection applies directly if the interface grows into another long, expensive chain. A recursively treated interface could help, but adds further complexity and storage.

For our FGS implementation, **we have not established that this trade is favorable**. The near-field graph has dense leaf-to-leaf blocks; eliminating interiors could produce large interface maps. Reusing those maps across sweeps would amortize setup, but would not eliminate their memory traffic. I would therefore treat this as a feasibility investigation, not an anticipated performance win.

Before implementing it, I would measure three things across increasing mesh sizes:

- **Current work versus critical-path growth:** does the available parallelism actually remain constant?
- **Interface size and predicted fill** for spatial partitions, particularly at 16/32/64 subdomains.
- **Estimated local-plus-interface execution cost**, including additional storage and traffic, against current dagteam.

There is also a useful distinction for your original proposal: **multiple threads per leaf can raise throughput without fixing its asymptotic scaling.** If leaf sizes remain constant while the chain grows, a four-worker leaf product may shorten each chain step by a bounded factor, but it does not make the chain wider. It is a practical improvement for the 64-thread target, not by itself a scalable algorithmic cure.

My expectation is therefore: cooperative leaf products offer the more direct near-term experiment; a transformed triangular solve offers a possible structural improvement **only if the interface analysis supports it**. If that analysis fails, parallel local GS with a coarse correction is the more promising direction to investigate, accepting a changed stationary iteration and requiring new convergence evidence.

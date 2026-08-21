# Theory — Sparse near-field ILU preconditioning for panel-method GMRES

**Status:** Phase 0 W6 design basis, revised 2026-08-12 for the implemented
FMM-direct-list ILU(0) preconditioner.

## 1. Setting and idea

The body solve is a dense boundary-element system

$$
G\,x = b, \qquad G \in \mathbb{R}^{N\times N},
$$

where $G_{ij}$ is the influence of panel $j$ on control point $i$ (normal velocity for
Neumann bodies, doublet potential for Dirichlet bodies — the same operator `apply_G!`
applies matrix-free). The matrix-free Krylov path never forms $G$; the FMM applies it in
$O(N)$. A preconditioner must therefore be built from something cheaper than $G$ itself.

The physical opening: panel influence kernels decay with distance (the doublet potential
kernel like $1/r^2$, its velocity like $1/r^3$), so after a spatial tree permutation the
strong local couplings tend to concentrate near the diagonal — the **near field**. This
locality need not be visible in the original mesh ordering. The FMM already
formalizes exactly this split: its tree partitions every interaction into a direct
(near-field) list and a multipole (far-field) list. Write

$$
G = S + F,
$$

where $S$ contains the near-field entries (tree-neighbor panel pairs, including the
self-influence diagonal and the strong same-leaf couplings) and $F$ the far field. The far
field is smooth and often numerically low-rank collectively; its individual entries are
not assumed to be small. $S$ is **sparse**: with mean near-field neighborhood
size $m \ll N$, it has $O(Nm)$ nonzeros.

The preconditioner is an incomplete LU factorization of $S$:

$$
S \approx \tilde{L}\tilde{U}, \qquad M^{-1} v = \tilde{U}^{-1}\big(\tilde{L}^{-1} v\big),
$$

applied by two sparse triangular solves ($O(\mathrm{nnz}(S))$ per apply — comparable to one
FMM near-field pass, far cheaper than an FMM iteration).

## 2. Why it should work here

Right-preconditioned GMRES solves $G M^{-1} y = b$, $x = M^{-1} y$. With $M = \tilde{L}\tilde{U} \approx S$,

$$
G M^{-1} = (S + F)\,M^{-1} \approx I + E_{\mathrm{ilu}} + F S^{-1},
$$

where $E_{\mathrm{ilu}}$ is the incomplete-factorization error and $F S^{-1}$ the
preconditioned far field. Favorable clustering and mesh-independent iteration counts are
benchmark hypotheses, not consequences of this split: neither term is guaranteed small
for this nonsymmetric panel operator. Two properties of our specific operator make $S$
carry most of the conditioning:

1. **Panel-size variation.** The rotor mesh's panel areas vary strongly root→tip; the
   self/neighbor influence magnitudes scale with panel size, and that variation lives
   entirely in $S$'s diagonal blocks. (This is also why plain diagonal equilibration is
   worth testing first as a cheap control — see Phase 2b ladder.)
2. **Trailing-edge/wake coupling.** For `RigidWakeBody`, the attached wake ties
   shedding-edge panels together with $O(1)$ entries (the linear Kutta structure). These
   are geometrically local to the TE strip and can be included in $S$'s pattern (§4).

This is standard practice in FMM-accelerated BEM: near-field (a.k.a. "mesh-neighbor" or
"pattern-sparsified") ILU and sparse-approximate-inverse preconditioners are the classic
remedies for GMRES on dense BEM systems (Carpentieri–Duff–Giraud's SPAI/Frobenius-norm
work for electromagnetic BEM; Benzi's preconditioning survey, §BEM; block-diagonal /
near-field ILU variants throughout the fast-BEM literature). Block Jacobi is the crudest
member of this family — it keeps only the block diagonal of $S$ and re-solves it exactly;
ILU keeps the whole near-field pattern approximately, which is why it typically wins when
inter-cell near-field coupling matters (our smoke result — Jacobi *hurting* — suggests
exactly that).

## 3. Construction

**Pattern.** Build dedicated target/source `FastMultipole.Tree`s and call the exported
`build_interaction_lists` with `Barba()`, fixed `leaf_size` and MAC, and
`nearfield=true`, `farfield=true`, `self_induced=true`. Expand each directed
`(target_branch, source_branch)` in `direct_list` into its Cartesian panel pairs. Map the
branch ranges through each tree's `sort_index_list`, add every missing diagonal, do not
symmetrize, and factor in one common spatial tree ordering. A linear allocation guard
rejects more than `max_pattern_entries=512N` requested pairs before triplet allocation or
kernel evaluation. There is no radius scan or $N\times N$ fallback.

**Entries.** $S_{ij}$ is computed by the *same* kernel path as `_G!`/`apply_G!` (the
side-aware `induced` self limit included), restricted to pattern pairs. Assembly cost is
$O(Nm)$ kernel evaluations — the same near-field work one dense assembly row-block does,
and embarrassingly parallel.

**Factorization.** Phase 0 uses **ILU(0)** from `ILUZero.jl`: zero fill preserves the
sparse pattern and gives predictable factor memory. ILUT is deferred and requires an
explicit fill cap before adoption. `ilu0` supports the 3-argument
`LinearAlgebra.ldiv!(y, F, v)` — i.e. they drop into the existing preconditioner contract
(`JacobiPreconditioner`/`FGSPreconditioner` precedent) with **no Krylov-side changes**:
wrap and pass as `N=`, `ldiv=true`.

**Routing.** Right preconditioning (`N=`), per the 2026-08-12 honesty ruling: the monitored
residual stays the true $\|Gx - b\|$, and stopping means what it says. (The apply is a fixed
linear map, so plain GMRES is fine — FGMRES not required.)

## 4. Wake and formulation subtleties

- **Attached-wake terms**: for shedding-edge panels the operator includes wake-strip
  contributions. Include them in $S$ for pattern pairs touching the TE strip (same kernel
  path — they come along for free if $S_{ij}$ is evaluated through `induced` on the
  constructed body). If omitted, the preconditioner is merely weaker at the TE, not wrong —
  $M$ never has to be exact.
- **Dirichlet bodies**: assemble $S$ from the doublet-potential operator with sources
  zeroed, mirroring `_apply_dirichlet_G!`; the pattern is unchanged.
- **Geometry changes**: $S$ and its ILU are geometry-frozen, like the FGS tree and leaf
  matrices. Rigid-body motion (the whole campaign) never invalidates them; a re-build hook
  is out of scope.

## 5. Cost model and expected accounting (rulings 7–8)

| component | cost | CSV column |
| --- | --- | --- |
| pattern + $S$ assembly | $O(Nm)$ kernel evals | `t_assembly` (of the preconditioner; noted) |
| ILU factorization | $O(Nm^2)$ in the usual sparse-row implementation; effectively linear only for bounded $m$ | `t_precond` |
| memory | $O(Nm)$ | `mem_state_bytes` (preserves the $O(N)$-vs-$O(N^2)$ story iff $m \ll N$) |
| apply (per GMRES iteration) | two sparse triangular solves, $O(\mathrm{nnz})$ | inside `t_solve_min` |

The apply is *serial-ish* (triangular solves) — expect weaker multi-thread scaling than the
FMM apply; this is a real benchmark axis, not a footnote.

## 6. Risks / failure modes

- **ILU instability** on nonsymmetric indefinite $S$ (zero/small pivots): reject
  missing/nonfinite/small pivots with a recommendation to use row equilibration or a
  diagonal shift. A small $\|SM^{-1}v-v\|$ is not required for ordinary ILU(0), because
  dropped fill can make that quantity substantial without making the preconditioner
  invalid.
- **Pattern too thin** (small `leaf_size` → tiny $m$): preconditioner degenerates toward
  Jacobi and inherits its failure. Knob: widen to two rings of neighbor leaves.
- **TE rows**: if wake terms are omitted and TE conditioning dominates, iteration counts
  stay high near the TE-coupled modes — diagnosable by comparing with/without wake terms
  in $S$.

## 7. Validation plan (hooks into the existing harness)

1. Unit: $S$ restricted-pattern entries ≡ dense `_G!` entries on the octa/diamond fixtures
   (both BC types); `ldiv!` wrapper linearity + side-effect-freeness (existing test
   patterns).
2. Solve: GMRES+ILU ≡ Backslash strengths at Phase-1 thresholds on a smoke mesh.
3. Benchmark (2b protocol): iterations & time vs {none, Jacobi-right, FGS, ILU} at matched
   true-residual targets across the ladder; `t_precond`, memory, and multi-thread scaling
   recorded per the schema.

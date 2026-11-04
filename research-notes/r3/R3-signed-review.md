# Adversarial review: R3 signed circuits

Reviewed November 4, 2026, America/Chicago. Target: `/home/chris/code/eth-research/notes/approaches/R3-signed.md`. The published-paper sources were read, not edited. No commits were made.

## Verdict

**The principal mathematics survives. The claims need tighter model boundaries.** In particular, I could not break the paired-cycle rectangles: their disjointness includes their zero cells, which is exactly what the approximate-rank proof needs. The normalization and factor of two in Theorem C are correct, even for signed, unnormalized approximants.

| Claim | Verdict | Qualification |
|---|---|---|
| A, A′ | **HOLDS** | For the stated smooth, genuinely vtree-respecting circuit model and the published geometric hypotheses. |
| C, C1–C4 | **HOLDS** | State the nonzero-denominator domain. No positivity or normalization of Q is needed. |
| S | **HOLDS-WITH-FIXES** | The displayed proof establishes the shared-amplitude, structured model and the stated sum over separate vtrees. It does not justify unrestricted use of the words “squared” or “PSD.” |
| B, B1, B3 | **HOLDS** | Minor-table construction is sound; make the omitted root and smoothing conventions explicit. This is a reformulation of known determinant DP, not a new algorithm. |
| K | **HOLDS-WITH-FIXES** | The citation exists and supports the tropical/polynomial claim. The extension to arbitrary nondecomposable Boolean-valued nonnegative PCs is not established. |
| M | **HOLDS-WITH-FIXES** | Correct with unary/indicator leaves, or an explicitly charged polynomial-size multilinear representation of every input function. Not with free arbitrary multivariate leaves. |
| X | **HOLDS-WITH-FIXES** | Correct for λ = sum f > 0, 0 ≤ ε < 1, and the stated input-arity cap. The definition of λ is essential. |

**Score: 8/10.** The strongest contribution here is C. S's core proof is also good. Most repairs are definitions and overclaims, not failed algebra.

**Severity summary:** no FATAL defect found in the restricted core theorems; three MAJOR scope/transfer defects; five MINOR statement/accounting defects below. An unrestricted interpretation of “squaring gains nothing” is false, with an explicit construction given below.

## Findings and concrete fixes

### MAJOR 1 — Define the input model once; the scope-only version of L1 is insufficient

Locations: §0–1 and Proposition M. Theorem X already states the needed input-arity cap.

L1 defines structure by saying every gate's scope is a vtree-node scope, then assumes that product children are the two corresponding vtree children. It also implicitly assumes a computation cannot terminate in an opaque multivariate input strictly above the chosen cut.

The latter implication is false under the scope-only definition. For example, let a vtree contain the subtree `{x1,x2}` and let a circuit input at that subtree compute the two-bit equality function. Multiply it by a unary function of x3 at the root. Every gate has a vtree-node scope, but at the `{x1}` cut there are no gates of that scope, while the output matrix can have rank two. The claimed expansion over `G_v` is then false. Unbounded-fanin products can likewise skip the intermediate scopes used in the proof.

The bigger issue is M if an arbitrary joint input function has unit cost: a single input can store UST. M's conversion of a leaf into a linear literal expression is then unavailable. This is not a subtle cancellation phenomenon. X explicitly excludes the full-scope input counterexample through its input-arity cap; its balanced-product argument does not require expanding each smaller input table.

**Fix:** stipulate unary/indicator inputs and binary, vtree-aligned products. Alternatively, define the vtree leaves as disjoint variable blocks, require circuit inputs to occur only at those leaves, charge their representation cost, and retain the arity cap where a balanced cut is needed. M needs the additional multilinear-representation assumption for a block input. State whether size counts gates, wires, or bounded-fanin gates.

The conventional unary-input interpretation makes the proofs work. Also, with the stated input-arity cap, the first balanced node on the heavy path cannot lie strictly below a larger opaque input; that can rescue the particular balanced-cut application without rescuing L1's assertion for **every** node.

### MAJOR 2 — S is a structured-amplitude theorem, not a theorem about every model called squared/SOS/PSD

Locations: verdict item 2, §4, §7's “any PSD factorisation” sentence.

The actual proof is stronger and cleaner than a square-expansion proof, but its hypothesis is specific: across the selected split, the output amplitudes lie in a shared left function space of dimension at most the charged circuit size. Neither an arbitrary arithmetic circuit followed by a square nor an arbitrary nested matrix-valued PSD architecture has been shown to satisfy that hypothesis with this cost.

The following coverage **is** justified:

* A normalized square of a smooth structured-decomposable real circuit, with small inputs.
* An explicit sum of such squares on separate vtrees, charging `S = sum_k s_k`, not silently counting their DAG union once.
* A PSD form `d(x)^T K d(x)/Z`, `K ⪰ 0`, whose full-scope outputs d come from one shared smooth structured circuit. Linear mixing preserves the cut space; factoring K does not require charging a fictitious expanded circuit.
* The compatible-square class `Σ²_cmp`, and the top-level `psd` class defined in Loconte–Mengel–Vergari's paper. Their Proposition 2 is precisely about this top-level compatible-output PSD model.

The paper also discusses PSD units **inside** monotone circuits and `μSOCS`, whose implicit number of squares can be exponential. A bound in terms of an explicitly expanded sum is not automatically a bound in terms of the original compact circuit. Positive-unital/quantum circuit architectures likewise require their own interface-to-size argument. The note correctly leaves the latter open in O4; the earlier blanket language should not override that qualification.

There is a concrete counterexample to dropping decomposability of the inner circuit. On Boolean edge indicators let

`H(x) = [z^(N−1)] product_e (1 − x_e + x_e z)`.

This is the cardinality indicator `[|x| = N−1]` and has a polynomial-size ordinary arithmetic circuit. Let

`d(x) = det L_r(x) · H(x) / sqrt(τ)`.

The determinant counts spanning trees contained in x. On an assignment with exactly N−1 edges it is one for a tree and zero otherwise. Consequently `d(x)^2 = UST(x)` exactly. The determinant and H have polynomial-size ordinary arithmetic circuits. Their product is not decomposable. Thus unrestricted polynomial-size squared arithmetic circuits represent UST exactly, including on expanders.

**Fix:** title and state S as a theorem for “normalized squares/PSD forms of structured-decomposable amplitudes,” list the covered literature classes explicitly, and state that inner decomposability is essential. Remove or separately justify “any PSD factorisation”: matrix PSD rank, a quadratic form over shared structured features, and a general PSD circuit are different size models. This is a proof-coverage objection, not a claim that every possible PSD-rank extension is false.

### MAJOR 3 — K's tropicalization transfer does not cover arbitrary Boolean nonnegative PCs

Location: §5, Corollary K, especially “even a non-decomposable one.”

A monotone **formal arithmetic circuit for the arborescence polynomial** tropicalizes correctly: positivity prevents cancellation, so its formal monomial support is the required family. This applies whether or not that formal arithmetic circuit is syntactically decomposable.

But a nonnegative circuit agreeing with a distribution only on Boolean assignments need not compute that formal polynomial after replacing indicator leaves by variables. The simplest witness is the Boolean computation `x · x = x`: it computes the one-edge tree indicator, but its formal polynomial is x², and its tropicalization returns twice the edge weight. Boolean equality does not justify the asserted homogeneous multilinear polynomial identity.

**Fix:** either restrict the PC transfer to smooth decomposable indicator PCs, where the network-polynomial argument in M supplies the missing identity, or state the lower bound directly for monotone formal arithmetic circuits computing the arborescence polynomial. A lower bound for arbitrary nondecomposable Boolean PCs needs an additional theorem. Do not present it as the two-line tropicalization argument.

### MINOR 1 — C's denominator and n convention

The abstract rank statement needs `alpha_g = Pr(Good ≥ g) > 0`; for g larger than the number of pairs it currently divides by zero. The closed form similarly needs n > 0, with the n = 0 case declared vacuous. State `0 ≤ ε < 1/2` when claiming a positive asymptotic lower bound.

The table uses `floor(r0/2)` while the full corollary uses `floor(ceil(r0)/2)`. Both yield valid lower bounds, but they can differ by one. Use the latter consistently, or explain that the former is a deliberate weakening.

### MINOR 2 — X's range is necessary before squaring

The local theorem should explicitly require `0 ≤ ε < 1` (or use `max(0,1−ε)^2`). The proof only establishes a lower bound proportional to `1−ε` before squaring; squaring is invalid if this is negative.

For example, choose two different spanning trees T1,T2 of a fixed graph and use the fixed-size indicator circuit `f_t = (1+t) δ_T1 − t δ_T2`. It has λ = 1 and Z+ = 1+2t. For uniform P on τ trees, ε = t+1−1/τ. The displayed unrestricted-in-ε claim would require `s beta_*² (1+2t) ≥ (1−ε)²`, whose left side grows linearly in t and right side quadratically. It fails for sufficiently large t. This does not affect the intended ε < 1 version.

### MINOR 3 — Cite L1 as an adaptation, not the literal statement of de Colnet–Mengel

The cited KR 2021 paper is real and the theorem numbers are correct. Its Theorem 2 decomposes a nonnegative output of a smooth decomposable arithmetic circuit into at most |C| balanced products, with a common partition in the structured case. Lemma 13 is rank subadditivity for those products. Neither is literally the “every vtree node, arbitrary signed output” statement of L1.

The note's own induction supplies that signed, specified-cut extension once the model is fixed. Say “the same decomposition/rank argument as Theorem 2 and Lemma 13.” Also, their “induced matching” is a graph-theoretic matching used to build a matrix gadget; do not identify it literally with every permutation submatrix.

### MINOR 4 — B's deleted root and circuit accounting

Define every minor-table interior set after deletion of the Laplacian root: `U_v = V_v \ (boundary_v ∪ {r})`, and use `boundary_v \ {r}` for the selectable row/column sets. Otherwise the root table literally asks for a row that does not exist. This is an omitted convention, not an obstruction to the construction.

Make clear that signed off-diagonal minors are allowed, and that 0/1 leaf functions are implemented as `δ0+δ1` and `δ1`. These ensure smooth scope E_v without changing the function. Size refers to the resulting bounded-fanin arithmetic/indicator circuit, not merely the number of distinct table states.

### MINOR 5 — K's root and semiring conventions

The published KN definition sums over out-trees with **all possible roots**, and the proceedings formulate tropical evaluation in the **(max,+)** semiring over real valuations. The note says minimum-weight arborescence and seems to use a fixed-root determinant.

These are harmless polynomial-factor adjustments, but state them: negate weights to exchange max/min, and sum the rooted directed matrix-tree constructions over all roots to obtain the cited polynomial. Use the directed incoming-Laplacian convention appropriate for out-trees.

## Proof audit

### A and the geometric reduction

All six steps of A check out. The union of the selected cycle paths is a forest, so it extends to a tree. The fundamental cycle of each inserted e_i is exactly its original C_i. Prior swaps do not remove any edges from an unswapped cycle. If a row/column pair selects both distinguished edges of any cycle, that cycle is present; if none selects both but at least one selects neither, the edge count is too small. Thus the claimed submatrix is exactly `τ^−1` times a permutation matrix.

The argument needs no positivity of a competing circuit. The published mixed-cycle lemma supplies at least ceil(r0) vertex-disjoint mixed cycles at the balanced cut, giving A′. The ceiling is legitimate because the cycle count is an integer. The geometric hypotheses, bounded degree, and short-cycle length remain essential to the later asymptotic conclusion.

### C1: swaps and orbit size

A good pair has one cycle missing e and the other missing f. Simultaneously toggling both orientations is a valid composition of fundamental-cycle exchanges. The cycles are vertex-disjoint, so this cannot alter a different pair's cycle-edge pattern. Hence the good-pair set is invariant under all allowed toggles.

The resulting equivalence class has exactly `2^g` distinct trees. Each pair contributes a genuine independent binary choice; there is no stabilizer identifying two choices.

### C2(i): permutation blocks

For a good pair, each allowed row contains exactly one of its two A-distinguished edges, and each allowed column contains exactly one of its two B-distinguished edges. Matching choices give the original tree or its paired exchange. A mismatched choice completes one of the two cycles and is therefore not a tree. With g pairs the tensor product gives a `2^g × 2^g` permutation pattern. The trees all have probability 1/τ.

### C2(ii): disjointness includes the zero entries — HOLDS

Here is a more explicit invariant-based proof of the note's weakest step.

For a pair q, define a property of **any cell** `(row,column)`:

1. all nondistinguished edges of both cycles are present;
2. exactly one of the pair's two A-distinguished edges is in the row;
3. exactly one of the pair's two B-distinguished edges is in the column.

Every cell of a block with q good has this property, including its off-permutation zero cells.

For a block with q not good, all edges of q are fixed at their values in its representative spanning tree. That representative cannot satisfy the three properties: the four choices for the A/B pair either put one distinguished edge on each cycle, making q good, or put both on one cycle, making the representative cyclic. Both contradict the assumption that it is a tree with q not good.

Thus a cell shared by two block rectangles forces their good-pair sets to coincide. For the same good-pair set, a shared cell also fixes every edge outside its toggle set identically. The representative trees then lie in the same toggle orbit, so the blocks coincide.

This is **cell-disjointness**, not necessarily row-disjointness or column-disjointness. Cell-disjointness is sufficient for adding entrywise ℓ1 errors. No assumption about Q enters this argument.

### C3 and transfer of the global approximation error

Let a block have size k = 2^h and permutation matrix Π. Its target is `P_B = Π/τ`. Restriction of the whole split matrix Q to this block has rank at most r. Therefore

`||P_B − Q_B||_entrywise-1 ≥ ||P_B − Q_B||_* ≥ (k−r)_+/τ`.

The second inequality is the nuclear-norm low-rank approximation bound for a matrix with k singular values all equal to 1/τ. Equivalently, multiply by Πᵀ and use a rank-r kernel; the note's partial-diagonal nuclear-norm argument is valid even though the error matrix need not be symmetric.

The normalization is **τ, the total number of trees**, not the number of trees in a block. A block's target mass is k/τ. For h ≥ g,

`||P_B − Q_B||_1 ≥ mass(P_B) · (1 − r/2^g)`.

If r ≥ 2^g the desired rank inequality is already automatic; otherwise all these terms are nonnegative. Summing only over the cell-disjoint blocks with h ≥ g gives

`2ε ≥ ||P−Q||_1 ≥ alpha_g · (1 − r/2^g)`.

Here the first inequality uses the approximation hypothesis; the middle expression is bounded below by the error on the selected cells. There is no claim that each block has error at most its own mass times 2ε. Q may be signed and may have arbitrary total mass. This is exactly the correct global-to-block transfer.

### C4: domination, including adjacent cycles

The local completion lemma conditions outside **all incident edges** of a cycle. Previously revealed good-pair bits and the other cycle's missing-edge event depend only on cycle edges whose endpoints are disjoint from the current cycle. Consequently they are measurable with respect to that outside sigma-field, even when an edge joins the two cycles.

For either orientation, apply the local bound to the first cycle, then to the second conditional on the first orientation and the earlier good-pair history. The joint probability is at least p². The two orientations are disjoint, giving at least 2p² for the next good indicator under every possible past.

Such sequential Bernoulli variables stochastically dominate independent Bernoullis of parameter 2p²; independence of the actual indicators is not assumed. With mean `2p²n`, the lower-tail Chernoff exponent `exp(−p²n/4)` at half the mean is valid. Rounding g downward only weakens the failure event. This proves the claimed asymptotic bound for fixed bounded-degree/short-cycle parameters and fixed ε < 1/2.

### S: both linear-algebra steps check out

For the vector-sign lemma, given a row vector u and a stacked column vector v, group v's coordinates by original column and apply Cauchy–Schwarz to the amplitude index. The remaining bilinear form is bounded by `|u|ᵀ W b`, with `||b||₂ ≤ ||v||₂`. Nonnegativity of W then gives operator norm at most `||W||op`. The vector-valued signs are not a hidden scalar-sign assumption.

For the amplitude matrix D obtained by concatenating outputs horizontally, shared cut gates supply a factorization D = F H, so rank(D) is bounded by the number of cut gates. Its squared Frobenius norm is exactly the normalizer Z, including after PSD mixing. Trace duality and `||D||_* ≤ sqrt(rank D)||D||F` give `BC(P,Q) ≤ beta sqrt(s)`.

The statement `1−TV ≤ BC` requires P and Q to be distributions. This is satisfied by `Q = sum d_m²/Z`, **Z > 0**, and does not require each inner circuit to be nonnegative or separately normalized. For separate vtrees, normalize each nonzero squared component, use the pointwise square-root subadditivity, then Cauchy–Schwarz over its mixture weights. The resulting charged size is `sum s_k`.

Do not transfer C's permission for arbitrary signed, unnormalized Q to S. Nor does S prove that monotone and squared circuits have identical optimal sizes: it gives them the same lower-bound expression, which is weaker than such an equivalence.

### B1, B, B3

B1's complementary row/column Laplace expansion is correct. Once the root is placed on the boundary, the selectable boundary has k−1 vertices. The number of pairs of equal-size subsets is

`sum_j binom(k−1,j)^2 = binom(2k−2,k−1)`.

Signs may depend on the selected row and column sets, but are constants, so they do not obstruct a separated representation. Treat k = 0/trivial splits separately rather than evaluating a negative binomial parameter.

For B, in an expansion of a parent minor of the sum of the two children's Laplacian contributions, any child term omitting a child-interior row or column vanishes: the other child has a zero row or column there. Surviving terms are exactly the stored child minors with complementary selected indices and their Laplace signs. The children use disjoint edge variables. Thus every product is decomposable and respects the edge-partition tree. There are `2^{O(w)}` states and transitions per node; indicator smoothing preserves this bound. This supports the claimed `M 2^{O(w)}` size, up to routine polynomial conventions.

B3's inequalities also check: choose the smaller strict-interior vertex set U_A; its exterior neighbors are boundary vertices. Expansion gives `h|U_A| ≤ Dk`, and `M/3 ≤ |A| ≤ D(|U_A|+k)` yields the displayed lower bound on k.

I did not independently reconstruct or rerun the author's determinant implementation. The independent numerical work below targets A/C/S as requested; B's verdict is based on its algebraic construction.

### M

With unary inputs, replace a leaf by `f(0) bar_z_e + f(1) z_e`. Smoothness and decomposability yield the unique full-scope network polynomial. Since the output is exactly UST on every Boolean assignment, its coefficients are zero off trees and 1/τ on trees. Substituting `bar_z = 1` and rescaling gives ST_G. The degree-(N−1) extraction is harmless, though redundant once this coefficient identity is established.

Conversely, homogenization of a syntactically multilinear ST circuit and zero-padding missing variables at sums produces an indicator circuit: padding means “absent from the represented monomial,” not the constant-one function. Products remain disjoint. This is a polynomial-overhead equivalence, not a polynomial-size construction for ST_G. It leaves the stated unstructured problem open.

### X: λ and cancellation normalization

Let b be the absolute-value circuit and `Z+ = sum b`. Pointwise `b ≥ |f|`. Its normalized output `B = b/Z+` is a nonnegative decomposable distribution, assuming Z+ > 0. The published mixture/product reduction gives

`BC(P,B) ≤ sqrt(s) beta_*`.

For `q=f/λ` and the note's definition `λ = sum f > 0`,

`BC(P,B) ≥ sqrt(λ/Z+) sum_x sqrt(P(x)|q(x)|)`

`≥ sqrt(λ/Z+) sum_x min(P(x),q(x))`

`≥ sqrt(λ/Z+) (1−ε)`.

The middle inequality remains valid when q(x) is negative: its right side is then negative. Finally

`sum min(P,q) = (1 + sum q − ||P−q||_1)/2 = 1−ε`,

because `sum q = 1` exactly. This identity remains true for signed q. Thus the proof is valid with the normalization actually stated in the note. For 0 ≤ ε < 1, squaring yields exactly the claimed inequality. Zero-mass product terms may be discarded in the nonnegative mixture decomposition.

**Do not change λ to an arbitrary positive external scale.** That would break the identity and the theorem. For example, take any fixed nonnegative exact circuit for P and scale its output to `f=tP`, but incorrectly hold λ = 1 as t tends to zero. Then ε = (1−t)/2, whereas Z+/λ tends to zero. The claimed lower bound tends to the positive constant `1/(4s beta_*²)`. Under the note's actual definition λ = t, the normalized output is P, ε = 0, and this spurious counterexample disappears. The ratio Z+/λ is scale invariant and is genuinely an absolute-circuit-mass versus signed-mass ratio.

## Independent numerical verification

I wrote fresh code without importing or copying the author's numerical routines:

* `/home/chris/code/eth-research/toy/r3/signed-review/check.py`
* `/home/chris/code/eth-research/toy/r3/signed-review/results.json`
* `/home/chris/code/eth-research/toy/r3/signed-review/run.log`

Run: `uv run /home/chris/code/eth-research/toy/r3/signed-review/check.py`. Seed 20261104. BLAS/OpenMP threads were capped at one, below the six-core limit. The job ran in the background and completed. No browser or lingering worker was started.

The code independently enumerates spanning trees by union-find, simple cycles by DFS, and maximum vertex-disjoint mixed-cycle packings by exact search. Ranks are computed modulo 1,000,003. A modular rank is a **certified lower bound** on real rank, so passing `rank_mod ≥ 2^nu` proves the tested real-rank inequality; I am not claiming modular rank always equals real rank.

| Graph | Trees | Splits | Certified rank checks | Maximum ν | Maximum good pairs reached |
|---|---:|---:|---:|---:|---:|
| Triangle | 3 | 3 | 3 | 1 | 0 |
| K4 | 16 | 31 | 31 | 1 | 0 |
| Triangular prism | 75 | 255 | 255 | 2 | 1 |
| Cube | 384 | 26 | 26 | 2 | 1 |
| Wheel on 6 vertices | 121 | 26 | 26 | 1 | 0 |
| K3,3 | 81 | 255 | 255 | 1 | 0 |
| Two triangles, bridge | 9 | 63 | 63 | 2 | 1 |
| Two triangles, bridge plus two chords | 64 | 255 | 255 | 2 | 1 |
| Three triangles, bridges plus two chords | 375 | 26 | 26 | 3 | 1 |
| Four triangles, bridges | 81 | 26 | 26 | 4 | 2 |
| Four triangles, bridges plus two chords | 1,575 | 26 | 14 | 4 | 2 |

For M ≤ 9, all nontrivial unordered splits were enumerated. Larger graphs used seeded random splits plus deliberately cycle-fragmenting/connector-biased splits. Twelve larger compressed matrices were skipped for the rank computation because a dimension exceeded 400; their block checks still ran. The largest ranked dimension was 364. These are small-graph tests, not an expander asymptotic experiment.

All checks passed:

* **980** certified rank lower-bound checks over **992** splits.
* **1,413** block families, checking **228,382** rectangle cells for pairwise non-overlap and exact permutation support. This includes fresh two-pair examples, not just single-pair tests.
* Conditional good-pair success was checked after every observed earlier good/bad history against the coarse `2p²` bound; no violations.
* **63** whole-distribution low-rank approximation/block-mass checks, using truncated-SVD Q without requiring positivity or normalization. Smallest slack in the tested error inequality: 1/3.
* **195** permutation-versus-rank nuclear/error checks; smallest slack zero, as expected for r = 0.
* **700** normalized shared-basis multi-output squared/PSD checks, including random non-diagonal PSD output mixing. Maximum `BC/(beta sqrt(r))`: **0.943402361**. Here r is the charged shared cut-space dimension, so this directly tests the linear-algebra core of S, not a claim about the size of an arbitrary PSD architecture.
* **300** vector-sign operator-norm checks, including four amplitude channels. Maximum lifted/original norm ratio: **0.770685780**.

The script raises on a violation rather than silently skipping a failed assertion. The JSON and log retain the per-graph coverage and the explicit rank-size skips.

## Citation and novelty checks

Live deep searches and primary-source fetches were performed November 4, 2026. Absence from these searches is not proof of novelty.

### Verified citations

1. **Kluk–Nederlof exists, including the ICALP 2026 publication.** The proceedings' Theorem 1.4 states: for every k ≥ 1 there is a graph of pathwidth at most k on k^{O(1)} vertices such that every tropical circuit calculating DST_G uses at least `2^{Ω(k log log k)}` gates. The definition is the all-root directed out-tree polynomial. This confirms the important quantitative assertion, subject to the transfer and convention fixes above.
   * https://drops.dagstuhl.de/storage/00lipics/lipics-vol374-icalp2026/html/LIPIcs.ICALP.2026.130/LIPIcs.ICALP.2026.130.html
   * https://arxiv.org/abs/2512.23121
2. **de Colnet–Mengel attribution is substantively right, literally overstated.** “A Compilation of Succinctness Results for Arithmetic Circuits,” KR 2021, Theorem 2 and Lemma 13 give the balanced-product decomposition and rank method, as described above. A and L1 are adaptations, not a newly discovered general rank technique.
   * https://proceedings.kr.org/2021/20/kr2021-0020-de-colnet-et-al.pdf
3. **Loconte–Mengel–Vergari's SOS/PSD definitions match the restricted coverage.** The compatibility definitions, §5, and Proposition 2 distinguish the top-level PSD quadratic form from deeper PSD/μSOCS constructions. This distinction should appear in the claim, not only in a literature footnote.
   * https://arxiv.org/html/2408.11778v3
4. **Newer squaring work does not itself supply this UST lower bound.** “How to Square Tensor Networks and Circuits Without Squaring Them” discusses structured decomposability, compatibility, orthogonality, and marginalization. It is relevant model context, not a proof of the note's universal wording.
   * https://arxiv.org/html/2512.17090v2

### Novelty assessment

* **C:** I found no prior UST statement combining paired mixed-cycle toggle orbits, disjoint full rectangles, and entrywise-ℓ1 approximate-rank lower bounds. The permutation low-rank approximation inequality and adaptive Bernoulli domination are standard; the worthwhile new content is their application through the paired-cycle partition. “Apparently new UST application; novelty not established exhaustively” is warranted.
* **S:** Square-root rank and compatible-output/SOCS expressivity bounds predate this note. I found no prior identical UST approximation bound. The norm proof gives a useful general affinity-versus-cut-space inequality, but neither the square-root-rank method nor PSD-to-squares conversion is new. Describe the contribution as this inequality/application, not a universal separation or a characterization of all PSD models.
* **X:** I found no identical UST absolute-mass statement. It is a short consequence of the existing nonnegative circuit reduction and an overlap inequality, close in spirit to cancellation/ε-sensitive monotone-computation arguments. This is a useful obstruction, but weak grounds for a major standalone novelty claim.

Other inspected current records included the 2026 separation-rank/PC expressivity paper (`2605.12940`), the approximation-hardness paper (`2506.01281`), and the positive-unital model (`2506.01824`). None supplied the missing general-model transfers. In particular, literature search cannot replace the formal-circuit versus Boolean-function distinction in K.

## Recommended edit order

1. Add a single circuit-model paragraph fixing inputs, fanin, structure, normalization, and size charging.
2. Narrow S's headline and list its covered classes; preserve its proof.
3. Restrict K's Boolean-PC transfer, and align the cited root/semiring convention.
4. Add the denominator/range/root bookkeeping fixes.
5. Keep A and C. Replace C2(ii)'s prose with the explicit cell invariant above if a reader is likely to doubt it.

No failed numerical check was hidden. No modifications were made to the claim note, the author's scripts, or `/home/chris/code/eth-research/notes/paper-ust/`.

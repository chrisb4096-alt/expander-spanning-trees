# Adversarial review of R3 constants

Reviewed 2026-11-04, America/Chicago. Claim reviewed: `/home/chris/code/eth-research/notes/approaches/R3-constants.md`. The claim, R2, and the paper source were read but not edited. No commit was made.

## Verdict: HOLDS-WITH-FIXES — 8/10

**The two main improvements survive review.** The telescoped conditioning lemma, random-deletion lemma, distance-3 expansion and spectral lemmas, closed-form exponent, and spectral constant all check out. In particular, I reproduce **c = 8.011299619 × 10⁻⁸** for D = 6 and λ₂ = 6 − 2√5. I found no fatal issue in the affinity theorem or its appropriately scoped circuit corollary.

**The unconditional hub-method impossibility conclusion is not established.** The note acknowledges the missing expansion lower bound in its table and body, but presents the resulting impossibility as settled in the opening summary and again in §7.2. There is also a false intermediate Cheeger statement, although the lower bound actually needed downstream remains valid.

Severity summary: **0 FATAL; 2 MAJOR; 5 MINOR.** The score concerns the whole note, not just the valid positive theorem.

## Findings and concrete fixes

### MAJOR 1 — The hub obstruction is conditional, not a proved general impossibility

Locations: opening summary, line 17; §7.2, especially lines 281 and 289–291; status row L-c.

The proved statements are the caps in terms of **L**. To turn them into an obstruction to a polynomial in **1/h**, one needs a family with **h ≥ c₀/L**, or at least the stated alternative spectral lower bound uniformly along the relevant family. Neither is proved here. Finite examples with λ₂L² around 1.5–2.4 do not supply that assertion for growing cores or arbitrary L.

Furthermore, the written description of the gadget does not constrain core sizes relative to the number of attachments. As written, a uniform h ≥ c₀/L need not hold at all:

* Take two equal, arbitrarily large 2-connected cubic expander cores, each with m vertices.
* Use n ≥ 2 hubs, a = b ≥ 2, with globally distinct attachment vertices and paths of length L.
* For S equal to the A core, |S| = m ≤ N/2 and |∂S| = na. Hence **h ≤ na/m**.
* Let m/(naL) tend to infinity. Then hL tends to zero, despite the cores themselves being expanders and the resulting gadget being 2-connected.

The supplied implementation avoids this particular scaling problem by taking core size approximately max(a,b)n and checking core 2-connectivity. Those restrictions need to be in the mathematical construction. The generic phrase “two cubic expander cores” is insufficient. Likewise, a single hub would be an articulation between the cores, so n ≥ 2 or an equivalent connectivity assumption is needed.

**Fix:** either label the no-go statement conditional everywhere, including the opening summary, or specify the well-attached, size-controlled core family and prove its expansion uniformly. Keep the L-dependent upper caps as proved results. The “within a factor of about 5 of the limit” interpretation must retain the same qualification. No change to Theorems P or S is needed.

### MAJOR 2 — The displayed Cheeger upper inequality is false as stated

Location: §3, lines 107–113.

The assertion, for every simple H on at least three vertices,

> i(H) ≤ √(λ₂(2Δ − λ₂))

has the exact counterexample **H = K₃**:

* Δ = 2;
* λ₂ = 3;
* i(H) = 2;
* the claimed right side is √3 < 2.

This is relevant, not merely outside the application: deleting one vertex from K₄ gives such an H. The cited Mohar–Poljak survey indeed prints this same overbroad upper inequality; citation does not repair the counterexample.

**Importantly, the conclusion used in R3, λ₂ ≥ Δ − √(Δ² − i²), is still valid.** When λ₂ > Δ it is automatic because the right side is at most Δ. In the λ₂ ≤ Δ regime one can use the appropriate Cheeger–Mohar argument. H with two vertices, which can also occur after deletion, should be handled explicitly; the required lower bound holds there as well.

**Fix:** state/cite the valid lower-bound form directly, or separate the λ₂ ≤ Δ and λ₂ > Δ cases and the two-vertex case. Do not downgrade the optimized constants unnecessarily: this repair preserves μ(D,h′), not just the weaker h′²/(2D). The closed-form result already uses only the weaker bound.

### MINOR 1 — The final coarse arithmetic inequality misses D = 3

Location: Corollary P′, step 4.

The last displayed comparison requires

> 4 ln(D+1) + 8.3(D+1) + 1.3 ≤ 10(D+1).

At D = 3 its left side is **40.04517744**, greater than 40. Thus this particular chain of rounded inequalities does not prove the endpoint. The intended bound itself is comfortably true: evaluating the preceding, unrounded expression at D = 3 fixes it. My independent checks of the original closed-form expression also leave substantial slack; the smallest log-slack in the tested grid is 8.95 nats.

**Fix:** handle D = 3 directly and use the coarse comparison for D ≥ 4, or tighten the rounded intermediate constants. Do not change the final exponent merely to accommodate a loose rounding step.

### MINOR 2 — The stated ε threshold includes a point where the exponential disappears

Location: §5, line 223.

At 1 − ε = exp(−cN/2), the claimed size lower bound is exactly 1, not exponential. Even strict inequality without a linear-in-N margin is insufficient to imply exp(Ω(N)).

**Fix:** use the paper's already correct formulation: require

> (1 − ε)² ≥ exp(−(c − c′)N) for some fixed c′ > 0.

Then the lower bound is exp(c′N). For fixed ε < 1 there is no issue.

### MINOR 3 — Feasibility claims need narrower wording and updated caveats

Location: §6, especially “Typical values” and the h = 0.3, D = 3 / h = 0.67, D = 4 rows.

The old Monien–Preis ceilings 1/3 and 4/5 are legitimate upper ceilings from those results, **not attainability claims**. In particular, the 2023 Lichev–Mitsche random-cubic bisection upper bound 0.139822N already gives h ≤ 0.279644 + o(1) for random cubic graphs. Thus h = 0.3 is not an asymptotically realistic random-cubic example even without any November 2026 preprint.

I also read the relevant statements of Gaspers, arXiv:2611.05146, submitted November 3. It claims bisection bounds c₃N + o(N) for all maximum-degree-three graphs, with c₃ = (3/(4π)) arccos(5/6) ≈ 0.1398, and N/3 + o(N) for maximum degree four. **If accepted**, these make both rows above asymptotically vacuous: the respective expansion ceilings become approximately 0.279644 and 2/3. I did not audit that large new preprint's proof, so these should be labeled new, unreviewed constraints, not silently treated as settled.

**Fix:** call the displayed intervals certified lower/upper bounds, not typical attainable ranges; distinguish random families from arbitrary bounded-degree families; flag the two rows under the new preprint. The numerical entries themselves are correct.

### MINOR 4 — The asymptotic comparison suppresses a term outside the small-h regime

Location: opening summary, line 11.

For h ≤ 1, the stated R2 comparison Θ((D/h) log D) is appropriate. It is not uniform over all h. For h = Θ(D), R2's term

> D log(2 + 16D³/h²)

is Θ(D log D), whereas (D/h) log D is only Θ(log D). The literal P′ expression with h₁ = min(h,1) also gives O(D log(D/h₁)), not by itself O(D log(D/h)) for every h.

**Fix:** state the comparison for 0 < h ≤ 1, which is the intended small-expansion improvement, or retain the full R2 expression and h₁. Other regimes may then be discussed separately. The actual bounds are not affected.

### MINOR 5 — Numerical certification and bibliography should be stated more carefully

Locations: numerical table implementation and §9.

* The optimized table is evaluated with binary floating-point `lgamma` and a bisection pushed to the S(q) = 1/2 boundary. The tail formula is a valid conservative bound, but ordinary floating-point arithmetic is not a formal outward-rounded certificate. My independent 65-digit recomputation agrees with the reported rounded entries. **Fix:** back q away from the boundary and bound the computation's error, or use interval arithmetic if these values are advertised as rigorously certified constants. The closed form does not have this problem.
* “There are no 2024–2026 results” is stronger than a literature search establishes. **Fix:** say no matching result was found in the stated search.
* The UAI 2020 paper's authors are **Honghua Zhang, Steven Holtzen, and Guy Van den Broeck**, not Zhang–Juba–Van den Broeck.
* arXiv:0710.2296 is **Sonny Ben-Shimon and Michael Krivelevich**, not Krivelevich–Sudakov et al.

## Proof audit

### Conditioning, telescoping, and the circuit reduction

1. **Spatial Markov conditioning is used on full edge configurations.** F is the union of stars at W\{w}; independence of W ensures no edge incident to w is in F. A positive-probability configuration cannot leave one of those conditioned vertices with no selected incident edge.
2. **The comparison survives every configuration, not just an average.** Before contraction, the unconditioned component containing w is H with the surviving edges at w. The other W vertices are isolated. Identifying an isolated vertex with its selected neighbours identifies vertices of this component; it does not remove conducting edges. No such identification contains w. Thus contraction only lowers the relevant resistance. Parallel edges produced by contraction are harmless.
3. **Every deletion prefix is legitimate.** Since b ≥ 1 and H is connected, H plus w retains a B edge throughout A deletion. Contracting this connected subgraph cannot disconnect it. Every required conditional avoidance probability is positive, and the product identity applies.
4. **The flow is a real flow in that comparison graph.** All Y_j are distinct neighbours in H, and u_j is not in Y_j because G is simple. The zero-sum demand χ_j has squared norm 1 + 1/m_j. An H flow meeting this demand followed by 1/m_j on each remaining edge into w has energy at most κ(1 + 1/m_j) + 1/m_j. Thus the factor is at least m_j/((1+κ)(m_j+1)).
5. **The telescoping and degree reductions are correct.** The product of m_j/(m_j+1) is b/d. Multiplying the A and B bounds gives (ab/d²)(1+κ)^(−d), and ab ≥ d−1 together with d ≤ D gives the stated uniform bound. Nothing requires regularity.
6. **Passing to the chain past is valid.** Each past of the selected star bits is a coarsening of full configurations on the other selected stars. The uniform conditional lower bounds survive averaging over those configurations. Selection and alteration of W use G and auxiliary randomness, not the sampled tree. Fixing a successful W before sampling T preserves R2's argument unchanged.
7. **The circuit reduction has the correct factor of two.** It uses δ²/max β², so β² ≤ exp(−cN) gives δ²exp(cN), not an exponent divided by two. This assumes the paper's circuit class: nonnegative smooth decomposable circuits with leaf arity at most 2M/3 (in particular ordinary univariate leaves) and its stated size convention. Without a leaf-scope restriction, a single leaf containing the entire UST distribution defeats any such claim. Carry that qualification into any standalone theorem statement.

### Distance-3 deletion

* **Expansion:** every surviving vertex loses at most one neighbour; for any S of size at most half of G−W, it was also at most half of G. Therefore at most |S| cut edges are lost. The inequality h(G−W) ≥ h(G)−1 is in the correct direction.
* **Spectrum:** extend a mean-zero unit vector on G−W by zero. It remains mean-zero on G. Its energy in G exceeds its energy in G−W by at most 1. Hence every test vector on G−W has energy at least λ₂(G)−1. Taking the infimum proves the claimed lower bound. This is not a mistaken use of principal-submatrix eigenvalue interlacing.
* **Resistance use:** if λ₂(G)>1, the bound makes G−W connected, and its actual unnormalized Laplacian bounds χᵀL⁺χ for zero-sum χ. Replacing 1/λ₂(G−W) by 1/(λ₂(G)−1) is conservative.
* **Packing:** a radius-two ball contains at most 1+D+D(D−1) = D²+1 vertices. The greedy factor is correct.

### Connected-set counting and random deletion

The numbered dangerous-pair definition in §3.3 is sufficient as written. In particular, **it explicitly permits S ∩ R ≠ ∅**. This is the condition most likely to be accidentally lost in a rewrite.

* The BFS encoding has exactly k−1 ones and j zeros, even when distinct explored vertices share an outside neighbour: an outside vertex is recorded only once. The fixed graph and exploration order decode the string uniquely. The binomial bound follows.
* Expansion bounds k by floor(b/h); outer vertex boundary is at most the edge boundary b. The two hockey-stick sums give exactly binom(K+b+1,b+1)−1. The upper bound Λ^(b+1) follows from the elementary binomial estimate.
* Conditional on w ∈ R, the other candidate vertices are still independent Bernoulli(q). For fixed S, at most b−1 other outer neighbours could be selected. A union bound over (m(b)−1)-subsets gives the claimed factor, even when different dangerous sets overlap. Independence between different dangerous-set events is **not** required.
* b ≥ 2 uses 2-edge-connectivity, which follows from 2-connectivity. The separate one-neighbour argument genuinely needs vertex 2-connectivity. Since |S| ≤ N/2, V\(S∪{w}) is nonempty, including at N = 3. Zero surviving boundary in G−w would make w an articulation.
* For |S| ≤ 1/h′, the surviving boundary ≥ 1 suffices. For |S| > 1/h′ and h′ ≤ h/(D+1), (h−h′)|S| ≥ D. Hence the other lower bound hk−D suffices. This closes the one-neighbour case.
* In the many-neighbour case, |N(S)∩W*| ≥ m(b), and W* ⊆ R. Thus any surviving neighbour w would have been removed by the **original** dangerous-pair definition. S is allowed to contain restored vertices in R\W*. No monotonicity of expansion under restoration is assumed.
* Linearity of expectation gives E|R*| ≥ q|W₀|/2; no concentration inequality is needed to obtain one sufficiently large outcome. R* remains independent because it is a subset of W₀. The surviving graph has at least two vertices: any deleted vertex has at least two distinct neighbours, all outside the independent deletion set.

### Closed form and assembled constants

Apart from the D = 3 rounding issue above, the proof is correct:

* m(b)−1 ≥ θb/2 and binom(b−1,m−1) ≤ 2^(b−1).
* With x = 2Λq^(θ/2), S(q) ≤ (DΛ/2)x²/(1−x).
* At q₀, θ = 1/(D+1), x² = 1/(2DΛ), and x ≤ 1/4 for D ≥ 3. Consequently S(q₀) ≤ 1/2.
* The h powers are 1 from the mixed-vertex density, 3(D+1) from q₀, and 2D from the spectral/flow factor: total **5D+4**, with no missing factor of two.
* M ≥ N follows from minimum degree at least two. The |Z| bound, independent-set factor D+1, alteration factor 2, and degree-quality factor are all accounted for.
* Substituting h = λ₂/2 in the mixed-vertex bound and using the distance-3 packing gives exactly

  c ≥ λ₂ / [3D(D+λ₂/2)(D²+1)] · (D−1)/D² · (1−1/λ₂)^D.

* Theorems H and I make the advertised conservative substitutions. No change of normalized versus unnormalized Laplacian occurs.

### Bare-path upper caps

Lemma L-a is valid: on deleting the first edge of a bare path, its remaining L−1 edges must carry any unit flow from the newly exposed endpoint to w. The resulting avoidance factor is at most 1/L. If an intermediate deletion disconnects the graph, the total avoidance probability is zero, which only strengthens the claim.

Lemma L-b is valid for the fixed star-bit chain. A row or column norm bounds the operator norm below; a conditional marginal is at most its unconditional average for at least one positive-probability past. The two sides need not use the same past because γ is a maximum over all pasts. The convex-secant inequality f(x) ≤ 4 ln(2)x on [0,1/4] is also correct.

For L ≥ 2 the stated caps in terms of L follow, provided the construction really has only the hubs mixed and the required balanced split. Equal a and b alone is not sufficient for balance if the two unspecified core sizes differ wildly; choose equal cores as the implementation does. The missing step is the uniform conversion from L to 1/h, not either probability lemma. The status-bit discussion remains a conjectural research direction, not a consequence of this review.

## Independent numerical verification

New code and outputs are confined to `/home/chris/code/eth-research/toy/r3/constants-review/`. I did not run or import the author's R3 scripts. The new tree enumerator uses its own union–find enumeration. Electrical calculations use **(L + J/n)⁻¹** on the mean-zero subspace, not a generic pseudoinverse with an unsafe cutoff. All commands used one BLAS/OpenMP thread and no parallel worker pool, well below the six-core limit. The first check run was launched in the background; all runs have finished.

### Exact finite graph checks

Results in `checks.json`:

* **538** 2-connected graph-atlas graphs on 3–7 vertices, including degree-two cases as extra boundary tests.
* **33,044** connected-set counting buckets checked; maximum count/binomial bound = **1**.
* **3,911** distance-3 deletion sets checked for both expansion and λ₂. Worst residuals were −2.22×10⁻¹⁶ and −4.44×10⁻¹⁵ respectively: numerical roundoff at equality, not substantive violations.
* **9,097** independent deletion outcomes checked using exactly the note's threshold-based dangerous-pair alteration. Every resulting graph retained h′ = h/(D+1).
* **17,537** dangerous-set witnesses included a vertex restored from R to G−R*. Thus the subtle restored-vertex case was actually exercised.
* **87** exact finite expectation/union-bound calculations across selected atlas graphs: q = 1/8, q = 1/2, and the admissible q₀. All union bounds passed; all **29 admissible-q₀** cases retained the asserted expected size. The larger-q tests are not claimed to satisfy S(q) ≤ 1/2.
* **23** tree-enumerated graphs, **4,193** spanning trees total. Graphs include K₂,₄, K₃,₃, K₆, wheels, a circular ladder, an unequal-arm theta graph, and three freshly seeded random cubic graphs.
* **20,620** positive-probability full-star conditioning checks and **30,416** positive-probability sequential-prefix checks. Minimum exact/closed-bound ratio = **1** to rounding; minimum exact/flow-bound ratio = **1** to rounding. No violations.

These finite checks supplement, not replace, the proof. They deliberately test independent W that need not be distance-3 separated, shared neighbours, contractions, rare conditioning configurations, and unequal degrees.

### Independent arithmetic

`arithmetic.py` uses exact rational h and h′ for floors and 65-digit arithmetic for the constants. It independently recomputes the Peierls polynomial through b = 300 and includes a conservative geometric tail. Its largest tail ratio among the following optima is below 0.38, so this cutoff is ample for the displayed digits. It searches the same 12-point h′ grid to make the table comparison like-for-like.

| D | h | R2 ln(1/c), recomputed | P ln(1/c), recomputed | Claim rounded correctly? |
|---|---:|---:|---:|---|
| 3 | 0.05 | 306.560218 | 73.646048 | Yes |
| 3 | 0.30 | 74.375440 | 47.965338 | Yes; feasibility caveat above |
| 4 | 0.67 | 69.323968 | 50.373946 | Yes; feasibility caveat above |
| 6 | 0.76 | 115.693197 | 75.075860 | Yes |

| D | λ₂ = D−2√(D−1) | Recomputed c_S | ln(1/c_S) |
|---|---:|---:|---:|
| 6 | 1.527864 | 8.01130×10⁻⁸ | 16.339828 |
| 8 | 2.708497 | 5.08862×10⁻⁷ | 14.491089 |
| 10 | 4 | 5.57560×10⁻⁷ | 14.399697 |
| 16 | 8.254033 | 2.46704×10⁻⁷ | 15.215075 |
| 32 | 20.864471 | 3.14200×10⁻⁸ | 17.275821 |
| 64 | 48.125492 | 2.78716×10⁻⁹ | 19.698242 |

For the D = 6 Ramanujan row, the recomputed R2 log reciprocal is **115.626714**, matching the claimed comparison. I also checked the final P′ expression on D ∈ {3,4,6,10,64,256}, h ∈ {0.001,0.05,0.5,1}; it is conservative throughout. This grid is not offered as a proof of an infinite parameter statement.

Reproduction:

```bash
cd /home/chris/code/eth-research
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  uv run --with numpy --with networkx toy/r3/constants-review/review_checks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  uv run --with numpy --with networkx --with mpmath toy/r3/constants-review/arithmetic.py
```

## Brief novelty check — November 4, 2026

I used deep web searches and read the relevant source passages, not just search snippets. No exact prior statement of the two new positive constants was found. That supports a tentative novelty assessment, not a proof of priority.

* **Connected-set counting is known.** Fomin–Villanger, *Treewidth computation and extremal combinatorics*, supplies the fixed-root connected-set/neighbourhood binomial count. The note's proof is a valid self-contained encoding argument. [Paper](https://fedorvf.github.io/articles/2012/2012n.pdf), [arXiv](https://arxiv.org/abs/0803.1321).
* **Percolation/alteration on expanders is established territory.** Ben-Shimon–Krivelevich, *Vertex Percolation on Expander Graphs*, studies deletion probability n^(−α) and the giant expanding component. It does not immediately give the same deterministic retained-set statement with the present explicit constant. [arXiv](https://arxiv.org/abs/0710.2296).
* **The flow lemma is elementary but potentially useful as packaged here.** Spreading electrical flow, resistance monotonicity, and sequential determinant/Schur-complement identities are standard. I found no identical b/d telescoped conditional-star statement in the checked sources; avoid claiming those component tools are new.
* **Exact UST circuit lower bounds already exist.** Zhang–Holtzen–Van den Broeck, UAI 2020, is directly relevant prior art, but does not provide this bounded-degree, approximate-TV theorem. [Proceedings and paper](https://proceedings.mlr.press/v124/zhang20c.html).
* **Recent approximation work is relevant but different.** *On the Hardness of Approximating Distributions with Tractable Probabilistic Models* studies bounded f-divergence hardness and an approximation separation between decomposable and deterministic-decomposable PCs. It is not this UST-expander constant result. [arXiv:2506.01281](https://arxiv.org/abs/2506.01281).
* **The October 2026 UST product-approximation paper is not a contradiction.** Mészáros–Peled couples a UST with a random 1-out subgraph with expected edge disagreement O(n log k/√k). Small expected edit distance is not a small total-variation guarantee and does not remove this circuit lower bound. [arXiv:2610.11554](https://arxiv.org/abs/2610.11554).
* **Feasibility references:** [Lichev–Mitsche, 2023](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v30i2p40); [Gaspers, November 3, 2026, new preprint](https://arxiv.org/html/2611.05146v1). Only the relevant stated bounds of the new preprint were checked, not its full proof.

## Files written

* `/home/chris/code/eth-research/notes/approaches/R3-constants-review.md`
* `/home/chris/code/eth-research/toy/r3/constants-review/review_checks.py`
* `/home/chris/code/eth-research/toy/r3/constants-review/arithmetic.py`
* Their outputs: `checks.json`, `checks.log`, `arithmetic.json`, and `arithmetic.log`, all in that same review-only directory.

I made no edits to the paper source, claim note, or original numerical scripts. The final working-tree check showed concurrent edits to original R3 numerical scripts and signed-circuit files; those edits are not mine and were not reverted. The claim note and paper source still match the versions reviewed. No browser was started and no background review process remains running.

# Adversarial review of R3-matroid

Reviewed November 4, 2026, US Central time.

## Verdict

**HOLDS-WITH-FIXES. Score: 8/10 as stated.**

I did not find a counterexample to Theorems A, B or C. The difficult conditional step in A works for every positive-probability past, not merely on average. The restriction to balanced matroids in B is justified, and S8 really violates the unrestricted conditioning comparison. The Fano construction also survives review, including the counting argument at balance M/4.

There are **0 FATAL, 1 MAJOR and 5 MINOR findings**. The major finding concerns the description of prior work, not a failure of the new inequalities. Several universal sentences in the definitions and duality discussion are literally false and should be corrected before this becomes a paper.

I read the claim, its numerical code, the imported tools and geometry arguments, and the R2 statement and conditioning argument. All new numerical checks used independent code, without importing the claim's toolkit or R2 code. The original code was used to identify test matrices and to compare results, not as the implementation of the tests.

## Findings requiring changes

### MAJOR 1. The description of arXiv:2109.06941 identifies the wrong polynomial

Section 12 calls its lower bound an exact lower bound for the **graphic basis-generating polynomial**. Theorem 1.1 of that paper concerns the **directed rooted spanning-tree polynomial** on the bidirected version of an expander. Its variables distinguish x_(u,v) and x_(v,u). Its displayed definition is set-multilinear by the non-root vertices.

Identifying the two directed variables produces the undirected polynomial, but a lower bound does not survive variable identification automatically. A hard polynomial can become easier after identification. Do not cite this theorem as the claimed undirected graphic-basis lower bound without an additional argument.

**Fix:** identify it as a directed, exact monotone-arithmetic lower bound. Distinguish it from both the undirected basis polynomial and TV approximation. Its discrepancy and epsilon-sensitive results remain relevant prior art, but neither is the theorem proved here.

Source read: [arXiv:2109.06941, definitions preceding Theorem 1.1 and Theorems 1.1–1.2](https://arxiv.org/html/2109.06941).

### MINOR 1. Not every circuit or cocircuit has a frame

The construction in Section 1 fails for a singleton loop K={e}. Here pi(e)=e and e belongs to cl(empty), so the proposed cut does not contain e. No cut contains a loop at all. Dually, a singleton coloop cannot have the proposed cycle frame.

**Fix:** say that every circuit and cocircuit **of size at least two** admits a frame. Singletons cannot cross a split, so this does not change A, B, C or the applications.

### MINOR 2. The graph interpretation confuses unions with binary sums

A matroid cycle in the stated definition is a union of circuits. It need not be an even subgraph. Two triangles sharing an edge have a union with odd-degree vertices. Their symmetric difference, not their union, is an even subgraph.

Likewise, a cut as defined here is a union of cocircuits and need not be a single graph cut delta(U). For example, all three edges of a triangle form such a union but are not delta(U) for any U.

**Fix:** describe graph cycles here as unions of ordinary cycles, and graph cuts as unions of bonds. Explain that binary cycle/cocycle spaces use symmetric difference and are narrower notions. Ordinary vertex stars used in the proof remain valid cuts.

### MINOR 3. The exponential corollary needs a lower bound on w0

Section 0 says that m >= c n and w0, Gamma = O(1) give exp(Omega(n)). The numerical exponent supplied by C is w0 m/(Gamma+1). An upper bound on w0 does not give a positive constant lower bound on that exponent divided by n.

**Fix:** require c>0 fixed, w0 >= c0>0 fixed, and Gamma=O(1), or say that w0 and Gamma are fixed positive parameters independent of n. The displayed finite-size inequality is unaffected.

### MINOR 4. One explanatory duality sentence reverses the wrong object

Section 9 says that a basis containing no circuit is equivalent to a cobasis meeting every cocircuit. The latter is false. In the triangle matroid, a cobasis is a singleton and misses the cocircuit consisting of the other two edges.

**Fix:** use the two correct statements: a basis contains no circuit and meets every cocircuit; its complementary cobasis contains no cocircuit and meets every circuit. Proposition 1 and the actual bit definitions already use the correct duality.

### MINOR 5. State the length hypothesis in the signed-circuit resistance example

After Lemma 7, the bound R <= (ell-1)/m requires that each of the m signed circuits have size at most ell. That restriction is not stated in the example itself.

**Fix:** add it. More generally the same disjoint-support construction gives R <= sum_j(|C_j|-1)/m^2. TU circuit coefficients can be normalized to +/-1 at the distinguished element, and averaging the resulting flows proves this directly.

## Mathematical audit

### Proposition 1 and the imported tools

Complementing a basis gives a basis of the dual. Complementing the row and column labels separately gives permutation-equivalent affinity matrices. Substitution at every leaf preserves scopes and hence smoothness and decomposability, and a bijection preserves TV. This is an actual same-size PC transformation, not just an equality of supports.

The imported data-processing inequality has the right direction. Coarsening each side can only increase beta. Its two maps act separately on the sides, including for a mixture of circuit and cocircuit pieces.

The adaptive chain-rule proof is valid. At the final step, the norms of the slices of f and g depend only on their respective earlier coordinates. The kernel may depend on the joint past. Zero-probability pasts contribute nothing. No independence between pieces is needed.

For the forbidden-cell kernel, trace one and determinant squared pq give the larger squared singular value

    (1 + sqrt(1-4pq))/2 <= 1-pq.

The asymmetric p and q used in R3 follow from the same calculation, even though the imported lemma phrases its final bound with equal lower bounds.

The circuit reduction applies to an arbitrary target distribution. Binarization, selection of a medium scope, nonnegative product decomposition, normalization of terms, and the overlap inequality do not use graph structure. It supplies splits with n/3 < |X| <= 2n/3. These are included in R3's closed balanced interval. There is no mismatch in balance or leaf arity. Use the stated nodes-plus-wires convention, or the imported node-only remark.

### Lemma 2: conditioning really gives a uniform minor

Positive probability makes sigma independent. A witnessing basis also shows that deleting tau after contracting sigma leaves rank r(M)-|sigma|. The forward and reverse basis maps are inverse bijections. Consequently all surviving bases have the same probability, rather than a multiplicity-dependent weight.

This qualification matters: an arbitrary infeasible assignment on F would not define a conditional law. The proof only uses feasible assignments.

### Lemma 3: the exchange bound holds for every matroid

The flat E-D can be extended to a hyperplane avoiding e, since e is outside its closure. Its complement is a cocircuit D' contained in D and containing e. A loop cannot occur here.

For each basis avoiding e, the fundamental circuit meets D' at another element f. Exchanging f for e gives a basis containing e. Given the image, a preimage is determined by f in D'-e, so there are at most |D'|-1 preimages. This proves the counting inequality and the stated probability bound. No representability or negative dependence is hidden in this argument.

The dual statement follows by complementing bases. This is a sound elementary lemma; I would not claim originality for it without a more focused literature check.

### Lemma 4: contraction and deletion both work

Let E_N=E-(sigma union tau), H=E-D and Z=H-(sigma union tau). For x in E_N, the contraction/deletion closure identity is

    x in cl_N(Z) iff x in cl_M(Z union sigma).

Deletion restricts the ground set; it does not change ranks of the remaining subsets. Since sigma is contained in H, the set Z union sigma is contained in H. Its closure therefore misses e. Also, cl_N(Z) contains all the remaining elements of H. Thus

    D_N = E_N - cl_N(Z)

is a cut containing e and contained in D-tau. Applying Lemma 3 gives the bound with denominator |D-tau|. The proof remains valid even for contractions that are not feasible basis conditionings; feasible conditionings are enough for the application.

The exact dual statement, worth writing explicitly, is:

    If C is a cycle, e in C and C intersect tau is empty,
    P_{U(M/sigma\tau)}[e not in B] >= 1/|C-sigma|.

In the dual, tau is contracted and sigma is deleted. The avoidance condition belongs to tau, not sigma. This verifies the deletion side as well as the contraction side.

### Lemma 5: every full past gets both bounds

For a proper subset S of a circuit, the induced runs in its cyclic order are paths. Ordering each run forward ensures that the frame for the current element avoids all earlier selected elements. The runs need not share a starting point or be independent random objects.

After conditioning on B intersect F=sigma and on inclusion of a prefix of S, the contracted elements are sigma plus that prefix. The frame misses all of them. Lemma 4 applies to this particular minor. Its positive lower bound makes the next prefix feasible, which closes the induction. The cocircuit case follows from the explicit dual version above.

This establishes a pointwise lower bound for **each** feasible assignment on F. It is stronger than an unconditional or average bound.

### Theorem A: the two probabilities are simultaneous

For a given piece, both halves are proper and nonempty because it crosses the split. Compatibility makes every frame of the new piece avoid the whole union of earlier piece supports. Overlapping territories are allowed; they do not cause a problem. The triangular-prism checks below deliberately include compatible pieces with overlapping territories.

For every full earlier configuration sigma, Lemma 5 bounds the first-half event and the second-half event under that same conditional law. Since their intersection is impossible, these are exactly the probabilities of cells (1,0) and (0,1). Both lower bounds therefore hold simultaneously.

A bit past is a coarsening of these full configurations. Averaging either pointwise inequality over the conditional distribution of sigma preserves it. In particular, one is not averaging a good bound for one event over a different past than the other event.

Data processing, the forbidden-cell calculation and the adaptive chain rule now give the claimed product. The argument never asserts that different pieces are conditionally independent. The pieces and their order are fixed before sampling the basis.

### Theorem C: greedy selection and PC reduction match

The conflict graph is undirected because compatibility requires avoidance in both directions. Restricting to the mixed pieces cannot increase its maximum degree. Each greedy choice removes at most Gamma+1 vertices, so the independent set has size at least ceil(m/(Gamma+1)). This is the correct rounding.

Apply A separately for each balanced split. The selected family can depend on the split; it cannot depend on the random basis. Uniformity over all balanced splits is exactly what the PC reduction requires. The exponential estimate follows from -log(1-w0) >= w0. The only correction needed here is the asymptotic wording in MINOR 3.

Under duality, cuts and cycles swap, supports and territories stay the same, and weights, conflict degrees and crossing counts stay the same. Together with Proposition 1, this proves the self-dual assertion. The intrinsic conflict estimate is a safe upper bound, although not optimized.

The direct-sum obstruction is also correct. Bases of a direct sum are independent uniform bases of its components. A split along a balanced union of components therefore has beta=1. Some condition ruling out such splits is indispensable.

### Lemma 6 and Theorem B: the direction of negative correlation is correct

For a nonloop, noncoloop f in any minor N', negative correlation gives

    P[e in B | f in B] <= P[e in B | f not in B].

Equivalently, contracting f lowers the inclusion probability relative to deleting f. Replacing deleted elements by contracted elements therefore gives a lower bound for circuit inclusion events. This uses negative correlation in each intervening minor, which is why balance, rather than negative correlation of U(M) alone, is required.

Loops and coloops give identical deletion and contraction matroids on the remaining elements. If the right-hand product is positive, induction ensures that the left-hand conditional prefixes exist. If it is zero, the requested inequality is trivial. The proof explicitly handles this issue.

Duality gives the deletion-all comparison for avoidance events. Applying these bounds at every full past proves B exactly as it proves A. No negative-association theorem is needed. The hypothesis that the pieces are disjoint ensures that the old supports do not contain the new half being exposed.

The S8 obstruction is real. In the binary representation [I4|A] with rows of A equal to 0011, 0101, 1001, 1111, use zero-based labels. There are 48 bases, with counts

    # containing 3 = 20; # containing 7 = 28; # containing both = 12.

Take F={7}, S={3} and observe 7 absent. Then

    P[3 in B | 7 absent] = 2/5 < 3/7 = P_{U(M/7)}[3 in B].

The ratio is 14/15. This reproduces the reported worst ratio, not merely a numerical near-violation. Our enumeration found four violating conditioned comparisons. It counts 5,600 comparisons because it counts each feasible assignment separately; the original 1,624 counts aggregate over assignments. These totals are not in conflict.

### Lemma 7: the resistance formula is valid, including degeneracies

For a full-row-rank TU representation, Cauchy–Binet assigns determinant squared one to each basis and zero to each nonbasis. The determinantal law is therefore uniform, not a nonuniform representable-matroid law.

After deletion, choose an independent subset of rows for the remaining columns. This retains TU and imposes an equivalent linear system on the coefficient vector x. If e is a coloop, that system without e is infeasible and both sides of the claimed identity are zero. Otherwise the columns excluding e span the row space, and Sherman–Morrison plus the least-norm solution gives the identity. A loop has a_e=0, R=0 and avoidance probability one, also correctly.

I checked every deletion subset of R10 independently below. The general formula does not extend to arbitrary real representations with uniform basis weights, and the note correctly restricts it to regular matroids.

## Applications and the Fano construction

The graphic specialization is sound. A cyclically ordered graph cycle is framed by its forward vertex stars, with weight the product of inverse vertex degrees. The uniform bound improves D^(-2(ell-1)) to D^(-ell), a factor D^(ell-2). For simple graph cycles ell>=3 this is a genuine improvement. The note correctly warns that its per-cycle comparison need not dominate every degree-sensitive version of v1.

For R2, 2-connectivity makes each vertex star a bond. Independent vertices have disjoint stars. B supplies the avoidance factors; the extra deletions in the R2 comparison can only increase resistance. Thus the claimed reduction to R2 is valid in the stated 2-connected setting. It does not claim to recover every arbitrary-graph star statement using cocircuits.

Cographic transfer is exactly Proposition 1 and should not be sold as an independent lower-bound mechanism.

### Section 7.3 survives the structural and counting checks

A 2-sum with F7 deletes both copies of the distinguished point and leaves six elements per core edge. A p-avoiding circuit or cocircuit of a gadget survives. Circuits through the identified points splice as claimed, and the cocircuit rule is its dual. A nonbridge core edge can be retained in a two-element circuit minor, giving F7 back after the 2-sum. Every core edge is a nonbridge by the linked-graph hypothesis. This justifies the asserted F7 minor and nonregularity.

The four p-avoiding Fano lines cover the six remaining points and intersect pairwise. If every line were monochromatic, their common intersections would force one color everywhere. Consequently every nontrivial gadget split has a mixed local line. Its size-four frames stay inside the gadget.

For a lifted core cycle, the local frame at a_e and the lifted star frame at b_e have exactly the required two-point intersections. If a graph star is not a bond, it is a union of bonds. Lifting these bonds gives a union of cocircuits, which is an allowed cut in A. Alternatively use the bond corresponding to the component of G-v containing the two cycle edges. Thus the distinction in MINOR 2 does not invalidate this construction.

Every lifted-star frame has size at most 3D, and the size-four local frames also obey that bound because D>=2. The product bound (3D)^(-2ell) follows. Vertex-disjoint core cycles give compatible lifts, including when their territories share a block on an edge between the cycles.

For the balance calculation, if c gadgets are split, each whole-red gadget contributes six red elements and each split gadget at most six. Since there are at least 2M red elements, there are at least M/3-c whole-red gadgets; the same holds for blue. Arbitrarily coloring the split gadgets cannot remove these. With c<epsilon M and epsilon<=1/12, each color has at least M/4 core edges.

The imported mixed-cycles lemma explicitly permits arbitrary m0>0. Substituting m0=M/4 and dividing by ell*kappa gives

    min{hN/(4D ell^2 kappa), M/(2 ell^2 kappa D(1+D/h))}

vertex-disjoint mixed core cycles. The epsilon printed in Section 7.3 makes this at least 2 epsilon M. At most c of these vertex-disjoint cycles contain a split gadget, since any such gadget corresponds to an edge and an edge cannot belong to two selected cycles. At least epsilon M clean mixed lifts remain. The alternative case already gives epsilon M compatible local pieces.

The displayed epsilon is uniformly positive for a bounded-degree linked expander family, since N/M>=2/D. This completes the asymptotic lower bound in the actual element count n=6M. I found no missing balance constant.

F7 is balanced, and binary balanced matroids are closed under 2-sums through their Rayleigh characterization. This example is genuinely nonregular, but it is not a witness that the all-matroid local theorem goes beyond the balanced class. The note's separate S8-gadget conjecture should remain separate.

### Uniform matroids are a valid limitation example

For U_(k,2k), the affinity matrix decomposes by j=|B intersect X|. Each nonzero block is constant, with squared norm

    binom(|X|,j) binom(2k-|X|,k-j) / binom(2k,k).

Thus beta squared is the largest hypergeometric atom. Across balanced splits it is Theta(k^(-1/2)), rather than exponentially small. This agrees with the smooth decomposable O(n^2) counting DP for the exact distribution. High Tutte connectivity alone cannot imply the desired exponential PC bound.

## Independent numerical checks

The checks ran with uv Python and NumPy/SciPy, one BLAS thread and at most one active numerical process at a time. Working directory was `/home/chris/code/eth-research/toy/r3/matroid-review/`. No script in the original numerical directory was imported or modified.

The small-matroid tester built bases with its own GF(2) rank routine, derived circuits and cocircuits independently, chose valid small frames in two cyclic orders, and used SVD of the exact enumerated probability matrix. It tested **every nontrivial split up to exchanging X and Y**, not just random or balanced splits. All generated compatible singleton and two-piece families were checked. These are finite checks, not exhaustive searches over all possible frames.

| Matroid | Elements | Bases | A comparisons | Compatible two-piece comparisons |
|---|---:|---:|---:|---:|
| Graphic K4 | 6 | 16 | 720 | 0 |
| Graphic triangular prism | 9 | 75 | 16,432 | 576 |
| Cographic triangular prism | 9 | 75 | 15,328 | 288 |
| R10 | 10 | 162 | 55,744 | 0 |
| F7 | 7 | 28 | 2,192 | 0 |
| F7* | 7 | 28 | 2,240 | 0 |
| AG(3,2) | 8 | 56 | 6,160 | 0 |
| S8 | 8 | 48 | 5,888 | 0 |
| Total | | | **104,704** | **864** |

No violation occurred. It is important that several small non-graphic matroids admit no compatible pair among these frames; their single-piece tests alone do not validate multiplication. The prism tests, the mixed-type test and the four-piece decorated test address that gap.

Additional results:

- **885,934 exchange/minor-certificate inequalities**, checking all cuts. Feasible conditionings were exhaustive for the tested matroids with at most eight elements and randomly chosen for the nine- and ten-element examples. Conditional probabilities were compared by integer arithmetic. No violations.
- **648 direct minor-law comparisons** against independently constructed contraction/deletion basis sets. No mismatched law.
- **71,273 frame inequalities** over proper subsets, admissible outside conditionings and both piece types. No violations.
- **32,256 full-past chain checks** for compatible prism pieces. Counts include repeated halves arising from different splits. No violations.
- A separate **three-piece circuit/cocircuit/circuit** test on F7 direct-sum F7* direct-sum a triangle checked all six orders on ten usable test splits, with every feasible full past. Its 60 bit-matrix bounds passed.
- Lemma 6 had no violations on the tested balanced examples for |F|,|S|<=2. S8 had the exact violation displayed above.
- All **5,120** pairs (deleted subset Z, surviving element e) of R10 satisfied Lemma 7. A separate search found and checked a TU signing. Maximum absolute numerical error was 3.28e-15.
- U_(k,2k), k=2 through 7, matched the hypergeometric formula at every balanced split size. For |X|=k the beta-squared values were approximately 0.666667, 0.450000, 0.514286, 0.396825, 0.432900 and 0.356935. The odd/even oscillation is expected and does not contradict the asymptotic statement.

### Beyond the 18-element decorated triangle

I constructed an independent binary representation of the **F7-decorated four-cycle**, with 24 elements and rank 11. For a core column v and fresh block coordinates a,c, its six columns are

    a, a+v, c, c+v, a+c, a+c+v       over GF(2).

Exhaustive rank testing of all 2,496,144 candidate 11-subsets found **196,608 bases**. The four possible block-size patterns are permutations of (3,3,3,2), each with 49,152 bases. This agrees with 4 * 16^3 * 12 from the 2-sum rule. It also independently checks which interface basis choice produces size two versus size three.

I checked all 16 local framed lines and one framed lift of the core cycle, including the cut property of every frame. Across all **7,852,453 balanced splits modulo complements**, this fixed 17-piece family always had a mixed piece; its minimum crossing count was one. This finite example is not evidence of expander growth by itself.

For 27 selected splits, including whole-block and within-block adversaries, sparse spectral calculations on the complete basis law checked 343 greedy compatible-family bounds. No violations; the largest beta-squared/RHS ratio was 0.500002. Four compatible local pieces also passed 2,400 full-past frame checks, with the smallest actual-probability/bound ratio 12/7.

One initial scratch-test assertion failed because **my tester reversed the two 2-sum interface choices** and predicted 4 * 12^3 * 16. The direct rank enumeration exposed that mistake. Correcting the tester gave the count above. This was not a failure of the claim's construction.

Temporary implementations and logs are `/tmp/r3-matroid-adversary.py`, `/tmp/r3-matroid-adversary.log`, `/tmp/r3-fano-review.py` and `/tmp/r3-fano-review.log`. The mixed-type and uniform checks ran from stdin. No long-running numerical process remains. No repository code or paper source was changed.

## Novelty check

I ran several rounds of deep web searches on November 4, 2026 and read the sources below. I found no direct prior theorem giving this framed circuit/cocircuit TV lower bound, or the nonregular Fano-decorated family under this model of approximation. This is a bounded search result, not a certification of novelty.

- [arXiv:2512.19515, Negations are powerful even in small depth](https://arxiv.org/html/2512.19515v1), especially Theorem 1.5 and its matroid discussion. It constructs an explicit matroid-related monotone Boolean function with a superpolynomial lower bound, alongside arithmetic separations. The function concerns whether an available set contains a basis, not approximation of the probability mass of the uniform basis law. It is important related work, but I found no theorem there that subsumes A/C. Avoid claiming the first matroid circuit lower bound of any kind.
- [arXiv:2109.06941](https://arxiv.org/html/2109.06941). Relevant directed exact arithmetic and discrepancy/epsilon-sensitive lower bounds, with the correction in MAJOR 1.
- [Jukna–Seiwert, Greedy can beat pure dynamic programming](https://arxiv.org/abs/1803.05380). Relevant earlier exact spanning-tree and tropical/monotone lower bounds. Again, an exact optimization lower bound is not automatically a constant-TV distribution lower bound.
- [Zhang–Holtzen–Van den Broeck, On the Relationship Between Probabilistic Circuits and Determinantal Point Processes](https://proceedings.mlr.press/v124/zhang20c.html). Relevant exact representational separation and spanning-tree context. The new claim should distinguish approximate representation from these exact DPP barriers.
- [arXiv:2506.01281, On the Hardness of Approximating Distributions with Tractable Probabilistic Models](https://arxiv.org/html/2506.01281). It addresses approximation notions and PC restrictions. It should be part of the distributional positioning; I found no general uniform-matroid-basis result there replacing the current argument.
- [Choe–Wagner, Rayleigh matroids](https://arxiv.org/pdf/math/0307096). This directly supports the binary no-S8 characterization, balance and the dual/minor/2-sum closure facts. It also explains why the Fano example does not leave the balanced class.
- [Mészáros–Peled, arXiv:2610.11554](https://arxiv.org/abs/2610.11554). The recent coupling controls the expected number of changed edges, not total variation between the unmodified distributions. It is not an approximate-PC upper bound that contradicts R3. I also verified the cited October identifiers [2610.11829](https://arxiv.org/abs/2610.11829) and [2610.33918](https://arxiv.org/abs/2610.33918); their stated topics do not supply the missing PC theorem.

The OpenReview copy of the approximation paper was blocked by a browser-verification page. I read the accessible arXiv copy instead; no conclusion relies on an unread search snippet.

The defensible novelty claim is the **all-matroid framed criterion for TV/affinity hardness, with a balanced-matroid conditioning version and the explicit nonregular construction**. The affinity toolkit, graphic/cographic duality, exchange counting and electrical identity should be presented as ingredients, not separately new results.

## Recommended disposition

Keep A, B, C and the Fano application. Correct the directed-polynomial citation and the five minor statements above. Preserve the explicit distinction between local all-matroid bounds, global balanced-matroid bounds, and the still-open problem of obtaining useful global resistance bounds for non-graphic regular families.

Only this review was written inside the repository. The paper-source tree was read only, and no commit was made.

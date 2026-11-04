# R3: the complete graph K_n (unbounded degree)

Written 2026-11-04 (CST). Helper for the "complete graph" direction. Scripts and logs are in `~/code/eth-research/toy/r3/complete/`.

## Verdict: PROVED. UST(K_n) costs 2^{Θ(n)} for TV-approximation by decomposable PCs (leaf arity ≤ 2M/3; reviewed in R3-complete-review.md, holds with fixes)

Notation: μ_n = UST(K_n) is the uniform labelled tree (Cayley/Prüfer) on the M = n(n−1)/2 edge variables. A *balanced split* is any edge split A ⊔ B = E(K_n) with |A|, |B| ≥ M/3. The split is arbitrary and need not be induced by a vertex partition. β(μ; A|B) = max over products q = q_A ⊗ q_B of BC(μ, q) = ‖√μ‖_op, as in `lem:LP`.

**Theorem K (PROVED, §3).** Let n ≥ 24 and r = ⌊n/24⌋. Every balanced split satisfies

  β(μ_n; A|B)² ≤ λ^{r/2} + (e/6)^{r/4} ≤ 2·e^{−r/150},

where λ = (1+√(1−4p))/2 and p = e^{−1}·(1/23)·e^{−4/23} ≥ 0.01344, so λ ≤ 0.98642.

**Corollary K (PROVED, from `cor:reduction`).** Let C be a smooth decomposable PC with nonnegative weights, leaves of arity at most 2M/3 and size counted as nodes plus wires. If TV(Q_C, μ_n) ≤ 1 − δ, then

  |C| ≥ (δ²/2)·e^{⌊n/24⌋/150} = 2^{Ω(n)}.

- **Error range.** This holds for every fixed TV < 1, and more generally for TV up to 1 − √2·e^{−⌊n/24⌋/300}.
- **No extra assumptions.** Determinism, exactness and any particular split structure are not assumed. Different product terms may use different balanced splits.

**Proposition U (PROVED, §3.6; CHECKED NUMERICALLY for n ≤ 12).** μ_n has an *exact* smooth, decomposable and deterministic PC with literal leaves and O(n⁴·3ⁿ) nodes plus wires.

**Corollary T (PROVED).** Within any fixed TV < 1, the smallest decomposable PC for UST(K_n) has size 2^{Θ(n)} = 2^{Θ(√M)}. Equivalently, max over balanced splits of β(μ_n; A|B)² equals e^{−Θ(n)}.
- The lower side is Theorem K.
- The upper side is Proposition U plus `cor:reduction` with δ = 1, which forces max β² ≥ 3^{−n}/O(n⁴).
- So the brief's target "2^{Ω(n)}" is the truth up to the constant. A bound exponential in the number of variables (2^{Ω(n²)}) is false.

**Constant.** The proved rate, β² ≤ 2e^{−n/3600} roughly, is very weak. Exact computation at n = 6–9 (§4) gives −2 ln β/n between 0.12 and 0.29 for every split tried, including local-search worst cases. At n = 9 the smallest value is 0.17 (star split). The ceiling from Proposition U is about ln 3 ≈ 1.1 per vertex.

**Key new idea.**
- Use *side-of-attachment bits* of vertices instead of edge or cycle bits:
  - a_w = [w has no A-edge to U] and b_w = [w has no B-edge to U];
  - w ranges over a small set W of "mixed" vertices (A-degree and B-degree both Θ(n)) that are isolated inside T[W].
  - Each bit aggregates about n/2 edge variables and is Θ(1)-balanced, unlike single edges (probability 2/n).
  - a_w·b_w = 0 because an isolated w must attach to U.
- Condition on the rectangle z = T ∩ E(W), using the averaged Lemma 1 below.
- The degree loss of R2 (D^{−D}) disappears by electrical networks:
  - In K_n the effective resistance between any two neighbours of w is O(1/n).
  - Building w's star one edge at a time multiplies the tree count by (j/(j−1))(1+1/m) per edge.
  - So Pr[w is a leaf attached through a B-edge | anything about the other W-vertices] ≥ d_B/(e·m), whatever the degree.

## 1. Literature (searched 2026-11-04: deep web search ×3 rounds, plus fetches)

| Work | What it proves | Relation to us |
|---|---|---|
| Martens & Medabalimi 2014, [arXiv:1411.7717](https://arxiv.org/abs/1411.7717), Thm 37 | Decomposable+complete SPNs computing, or converging to, the K_m spanning-tree indicator need size ≥ 2^{m/30240} (m ≥ 20). Proof: s² weak (balanced-product) functions, each nonzero on a 2^{−Ω(m)} fraction of trees. | **Support / exact-limit only.** "Approximated arbitrarily well" means a sequence converging pointwise. No fixed-TV statement, and a TV-approximator can put mass on non-trees. Exponent is in vertices, like ours. They give no upper bound. |
| Zhang, Holtzen & Van den Broeck, UAI 2020, [PDF](https://starai.cs.ucla.edu/papers/ZhangUAI20.pdf) | UST(K_n) is a DPP; restates MM for smooth decomposable PCs. | Exact only. |
| Sharir & Shashua, SPQN, [PMLR v84](https://proceedings.mlr.press/v84/sharir18a/sharir18a.pdf) | Adapts MM to a complete-graph distribution. | Exact only. |
| Chattopadhyay, Datta, Ghosal & Mukhopadhyay, ITCS 2022, [arXiv:2109.06941](https://arxiv.org/abs/2109.06941) | Monotone arithmetic circuits: 2^{Ω(n)} for expander spanning-tree polynomials. ε-sensitive bounds for F_{n,m} ± ε·(spanning-tree polynomial of K_n) via a universal low-discrepancy distribution over nearly balanced *vertex* partitions (IP gadget). | Different model: polynomial-level perturbation, monotone circuits, discrepancy. Not TV approximation of the tree law. Their "balanced" is over vertices of a rooted functional representation. |
| Leland & Choi, NeurIPS 2025, [arXiv:2506.01281](https://arxiv.org/abs/2506.01281) | NP-hardness of bounded-f-divergence approximation for any tractable-marginal model. 2^{Ω(n)} TV-approximation bound for *deterministic* decomposable PCs (Sauerhoff function, via de Colnet–Mengel weak approximation and Bova et al. rectangle partitions). | No unconditional bound for *non-deterministic* decomposable PCs on an explicit distribution. No overlap. |
| de Colnet & Mengel, IJCAI 2020, [arXiv:2011.13721](https://arxiv.org/abs/2011.13721); Bova, Capelli, Mengel & Slivovsky (DNNF vs d-DNNF) | Weak/strong approximation lower bounds for d-DNNF; rectangle-partition (determinism) arguments. | Need determinism. Our argument does not. |
| Chubarian & Turán ([PMLR v138](https://proceedings.mlr.press/v138/chubarian20a.html)) | *Upper* bounds: OBDD approximations of bounded-treewidth BN classifiers. | Context only. |
| Mészáros & Peled, Oct 2026, [arXiv:2610.11554](https://arxiv.org/abs/2610.11554) | Couples UST(G), min degree k, with the random 1-out graph so that E\|symmetric difference\| = O(n log k/√k). | A "product over vertices" approximation in edit distance only; a coupling bound says nothing about TV by itself. TV is ≈ 1 for a separate reason: the 1-out graph is a spanning tree with probability n^{n−2}/(n−1)^{n−1} ~ e/n, and is uniform given that, so TV = 1 − O(1/n). No conflict with Theorem K. It does say vertex-local structure is what UST(K_n) looks like, which is what our bits use. |
| Anari–Liu–Oveis Gharan, Chen–Eldan, entropic independence (e.g. [arXiv:2106.04105](https://arxiv.org/abs/2106.04105), [FOCS 2022](https://ieee-focs.org/FOCS-2022-Papers/pdfs/FOCS2022-4Bu7jGV9xIcveUWYj3oWoi/551900a123/551900a123.pdf)) | Correlation decay / MLSI for spanning trees under pinnings. | Not needed. Rayleigh monotonicity under pinning plays the same role and is elementary. |

Also checked and not relevant: [2609.18707](https://arxiv.org/abs/2609.18707) (FPRAS for TV between product mixtures), 2605.12940 (PC vs LLM expressivity), sum-of-squares circuits (Loconte–Mengel–Vergari, [AAAI 2025](https://ojs.aaai.org/index.php/AAAI/article/view/34100); no spanning trees), and quantum spanning-tree sampling [2609.40314](https://arxiv.org/abs/2609.40314).

**Novelty status.**
- **Theorem K.** I found no prior fixed-TV (or any non-vanishing-error) lower bound for UST(K_n) against non-deterministic decomposable PCs. The ingredients are standard electrical-network facts. The combination is new as far as these searches go: side bits plus the averaged rectangle conditioning plus the leaf-via-side bound.
- **Proposition U.** I did not find the 3ⁿ subset-DP circuit stated in the PC literature. It is probably folklore (subset DP for rooted trees), so cite it as "easy".

## 2. Setup

Fix a balanced split A ⊔ B. Write θ_v = deg_A(v)/(n−1), where deg_A counts A-edges of K_n at v. Choose W ⊆ V with |W| = r and put U = V∖W, m = n − r. For w ∈ W:
- S_A(w) = {wu ∈ A : u ∈ U}, S_B(w) = {wu ∈ B : u ∈ U};
- d_A = |S_A(w)|, d_B = |S_B(w)|, so d_A + d_B = m;
- a_w = 1[T ∩ S_A(w) = ∅] is a function of x_A, and b_w = 1[T ∩ S_B(w) = ∅] is a function of x_B.

Let z = z(T) = T ∩ E(W), where E(W) is the set of edges with both ends in W. Then z = (z ∩ A, z ∩ B), with each part a function of one side, so {z(T) = z} is a rectangle. Let W'_z be the set of w ∈ W with no edge of z at w. For w ∈ W'_z, the tree must connect w to U, so a_w·b_w = 0 (one forbidden cell).

Standard facts used:
- **(F1) Spatial Markov.** Conditioning UST(G) on a consistent in/out status of an edge set F gives UST of G' = (G − out)/in on the undecided edges. Loops are removed; parallel edges are kept.
- **(F2) Edge identities.** For e not a loop or bridge: τ(G) = τ(G−e) + τ(G/e), τ(G/e)/τ(G−e) = R_{G−e}(e) (Kirchhoff), hence Pr_G[e ∉ T] = 1/(1+R_{G−e}(e)). Also τ(H+e) = τ(H)(1 + R_H(e)).
- **(F3) Rayleigh monotonicity.** Effective resistance does not increase when edges are added or vertices are identified (shorted).
- **(F4) Exact resistance.** Let H_d be K_m on U plus a vertex w joined to d ≥ 1 vertices D ⊆ U, and let u ∈ U∖D. Then R_{H_d}(w,u) = 1/d + (d+1)/(d·m).
  Proof: inject a unit current at w and extract it at u. By symmetry, D-vertices sit at potential φ_D and the other m−d−1 vertices at φ_O, with u at 0. Kirchhoff at an O-vertex gives φ_O = dφ_D/(d+1). At a D-vertex: 1/d = (m−d−1)(φ_D−φ_O) + φ_D = φ_D·m/(d+1). Finally φ_w = φ_D + 1/d. ∎

## 3. Proofs

### 3.1 Lemma 1: averaged rectangle conditioning (PROVED)

Let μ be a law on X×Y, and let z = (φ(x), ψ(y)) for maps φ, ψ. Write μ_z = μ(·|z). Then

  β(μ; X|Y)² ≤ Σ_z μ(z)·β(μ_z; X|Y)².

*Proof.*
1. Let q = q_X ⊗ q_Y. The cell R_z = φ^{−1}(z₁) × ψ^{−1}(z₂) is a rectangle, and the cells partition X×Y.
2. q restricted to R_z and renormalised is a product q_z. So BC(μ,q) = Σ_z √(μ(z) q(R_z))·BC(μ_z, q_z) ≤ Σ_z √(μ(z) q(R_z))·β(μ_z).
3. Cauchy–Schwarz and Σ_z q(R_z) = 1 give the claim. Maximise over q. ∎

This lemma is not in `sec-tools.tex`. It is needed because the set of usable coordinates W'_z depends on z. Making z the first coordinate of `lem:chain` would be wrong: the chain maximises over pasts, and some pasts make every later kernel trivial.

### 3.2 Lemma 2: mixed vertices (PROVED)

Let 0 < c < 1/9 and Mid = {v : c ≤ θ_v ≤ 1−c}. For every balanced split, |Mid| > (1/3 − 3c)·n.

*Proof.*
1. Let L = {θ_v < c} and H = {θ_v > 1−c}. They are disjoint and |E(L,H)| = |L||H|, since the graph is complete.
2. Each L–H edge is an A-edge at its L end or a B-edge at its H end. So |L||H| ≤ Σ_{l∈L} deg_A(l) + Σ_{h∈H} deg_B(h) < c(n−1)(|L|+|H|) (when L∪H ≠ ∅). Hence min(|L|,|H|) < 2cn.
3. Case |H| < 2cn: n(n−1)/3 ≤ 2|A| = Σ_v deg_A(v) < (n−1)(c|L| + |H| + |Mid|) < (n−1)(3cn + |Mid|).
4. Case |L| < 2cn: the same with B, using |B| ≥ M/3 and deg_B < c(n−1) on H. ∎

### 3.3 Lemma 3: per-vertex bounds under conditioning (PROVED)

Fix a possible value z with w ∈ W'_z. Fix any possible value σ of T on the stars ∪_{x∈W'_z∖{w}} (S_A(x) ∪ S_B(x)). Write P = Pr[· | z(T)=z, σ]. Then

- (i) leaf bound: P[a_w = 1] ≥ (d_B/m)(1+1/m)^{−(m−1)} > d_B/(e·m);
- (ii) deletion bound, if d_B ≥ 1: P[a_w = 1] ≥ (1 + 1/d_B + (d_B+1)/(d_B m))^{−d_A} ≥ exp(−d_A/d_B − 2d_A/m).

The same bounds hold with A and B swapped.

*Proof.* By (F1), P is UST(G') for the minor G' of K_n obtained by deciding E(W) according to z and the stars in σ. Two facts about G' are all we use:
- (α) w's edges to W are decided out, since w ∈ W'_z. w's m edges to U are undecided, since w's star is not in σ. So w is never identified with another vertex, and its incident edges in G' are exactly its m U-edges.
- (β) All U–U edges are undecided. So G' is obtained from a supergraph of K_m ∪ (w's U-star) by deleting other edges and identifying vertices. (F3) therefore lets us bound resistances in G' by resistances in K_m ∪ (part of w's star).

*(i).*
1. Let G'' = G' − w. It is connected: the U-images form a connected image of K_m; every other isolated W-vertex is glued to U by σ (it has a T-edge to U); every non-isolated W-component keeps its undecided U-edges.
2. List w's edges as e_1, …, e_m with endpoints u_1, …, u_m, which are distinct vertices of U in K_n. For each B-edge e, the event "T ∩ star(w) = {e}" has probability τ(G'')/τ(G'), and it implies a_w = 1 by (α). So P[a_w=1] ≥ d_B·τ(G'')/τ(G').
3. By (F2), τ(G') = τ(G''+e_1)·Π_{j=2}^{m} (1 + R_j), where R_j = R_{G''+e_1..e_{j−1}}(w, u_j), and τ(G''+e_1) = τ(G'') because w is a leaf there.
4. By (β) and (F3), R_j ≤ R_{H_{j−1}}(w,u_j). Here H_{j−1} is K_m plus w joined to u_1..u_{j−1} in K_n, where these vertices are distinct and u_j is not among them. By (F4) this equals 1/(j−1) + j/((j−1)m).
5. So 1 + R_j ≤ (j/(j−1))(1 + 1/m), and τ(G')/τ(G'') ≤ m(1+1/m)^{m−1} < e·m. ∎

*(ii).*
1. Let S_A(w) = {e_1..e_{d_A}}. Then P[a_w=1] = Π_i P[e_i ∉ T | e_1..e_{i−1} ∉ T].
2. By (F1) and (F2), the i-th factor is 1/(1 + R_{G'_i}(w,u_i)), where G'_i = G' − e_1 − … − e_i. It still contains the image of K_m ∪ star(w, N_B(w)) (all d_B B-edges remain).
3. By (F3) and (F4), R ≤ 1/d_B + (d_B+1)/(d_B m) ≤ 1/d_B + 2/m.
4. Then use 1+x ≤ e^x. ∎

**General form (PROVED; not needed for K_n).** For any graph G, let ρ = max_{u,u' ∈ N_U(w)} R_{G[U]}(u,u') and d = deg_U(w). The same argument gives P[w is a leaf through a B-edge] ≥ (d_B/d)(1+ρ)^{−(d−1)}.
- Thomson's principle bounds R_j by 1/(j−1) + ρ: spread the unit flow evenly over the j−1 star edges, then route each share by the u_i→u_j unit flow in G[U], and use convexity of energy.
- The degree cost is e^{−ρd}, which is O(1) whenever neighbour resistances are O(1/d).
- This is exactly where R2's D^{−D} loss goes away for dense graphs. For K_m it gives the weaker (1+2/m)^{−(m−1)}.

### 3.4 Lemma 5: few tree edges inside W (PROVED)

For S ⊆ E(K_n), Pr[S ⊆ T] ≤ (2/n)^{|S|}. Two proofs:
- by transfer-current determinants and Hadamard's inequality, since diagonal entries are Pr[e∈T] = 2/n;
- directly from the forest formula: the number of trees containing a forest with component sizes c_1..c_t is (Π c_i)·n^{t−2}, and Π c_i ≤ 2^{|S|}.

Hence, for j ≥ 1, Pr[|T ∩ E(W)| ≥ j] ≤ C(C(r,2), j)(2/n)^j ≤ (e r²/(jn))^j. ∎

### 3.5 Proof of Theorem K

1. Take c = 1/12 in Lemma 2, so |Mid| > n/12 ≥ 2r. Let W be any r vertices of Mid, with r = ⌊n/24⌋ and ε := 1/24.
2. For w ∈ W: d_A ≥ deg_A(w) − (r−1) ≥ c(n−1) − r + 1 ≥ cn − r, and likewise d_B ≥ cn − r. With m = n − r ≥ (1−ε)n, θ_A := d_A/m ∈ [(c−ε)/(1−ε), 1 − (c−ε)/(1−ε)] = [1/23, 22/23].
3. **Per-vertex kernel.** By symmetry assume θ_A ≤ 1/2. Lemma 3(ii) gives p_A ≥ exp(−θ_A/θ_B − 2θ_A) ≥ e^{−4θ_A}, and Lemma 3(i) with roles swapped gives p_B ≥ θ_A/e. So p_A·p_B ≥ θ_A e^{−1−4θ_A}. On [1/23, 1/2] the minimum is at θ_A = 1/23 (f(θ) = θe^{−4θ} rises to θ = 1/4 and then falls, and f(1/2) = 0.0677 > f(1/23) = 0.0365). This gives p_A·p_B ≥ p = 0.01344.
4. **Chain, for fixed z.** Order W'_z arbitrarily; it is fixed once z is fixed. By `lem:dpi`, β(μ_z; A|B) ≤ ‖√(law of ((a_w)_w, (b_w)_w) under μ_z)‖_op. Apply `lem:chain` to the pairs (a_w, b_w), w ∈ W'_z.
   - Each past (the bits of earlier vertices) is a function of σ for the current w. So by Lemma 3, averaged over σ, the current kernel K satisfies K(1,0) = P[a_w=1 | past] ≥ p_A(w) and K(0,1) = P[b_w = 1 | past] ≥ p_B(w), with K(1,1) = 0.
   - `lem:two-by-two` (formula form) gives γ_w² ≤ (1+√(1−4p_A p_B))/2 ≤ λ. Hence β(μ_z)² ≤ λ^{|W'_z|}.
5. **Averaging.** By Lemma 1, β² ≤ Σ_z μ(z)·λ^{|W'_z|}. Note |W'_z| ≥ r − 2|z|. If |z| < ⌈r/4⌉ then |W'_z| > r/2. By Lemma 5 with j = ⌈r/4⌉ ≥ r/4, the remaining z have total mass ≤ (4e r/n)^{j} ≤ (e/6)^{r/4}. So β² ≤ λ^{r/2} + (e/6)^{r/4}.
6. **Numbers.** λ^{1/2} = e^{−0.00684} and (e/6)^{1/4} = e^{−0.198}, so both terms are ≤ e^{−r/150}. ∎

Corollary K is `cor:reduction`: |C| ≥ δ²/max β².

**Other balance thresholds.** For splits with |A|, |B| ≥ αM, the same proof runs with c < α/3 and |Mid| > (α − 3c)n, so any constant α > 0 still gives 2^{Ω_α(n)}.

### 3.6 Proposition U: exact 3ⁿ circuit (PROVED; CHECKED NUMERICALLY)

For S ⊆ V and v ∈ S, let f(S,v) be the indicator, over the variables E(S), of "the edges inside S form a spanning tree of S" (rooted at v for the recursion). Put f({v},v) = 1, and for |S| ≥ 2 let s_0 = min(S∖{v}). Then

  f(S,v) = Σ_{S' ⊆ S∖{v}, s_0 ∈ S'} Σ_{u ∈ S'} f(S',u) · [x_{vu}=1] · f(S∖S', v) · Π_{e ∈ E(S',S∖S')∖{vu}} [x_e = 0].

- **Exactness and determinism.** Removing v's edge to the subtree containing s_0 splits the tree uniquely into (S', u, rest). So the sum is a partition of the support: the circuit is exact and deterministic.
- **Decomposable and smooth.** The four factors have disjoint scopes E(S'), {vu}, E(S∖S') and the cut, whose union is E(S). Every summand has scope E(S).
- **Size.** Σ_s C(n,s)·s·Σ_t C(s−2,t−1)·t·O(1 + t(s−t)) = O(n⁴3ⁿ). The root is f(V, 0) scaled by n^{−(n−2)}.

`upper.py` checks the following:
- for n ≤ 7, the generated trees are exactly all labelled trees, each generated once;
- for n ≤ 12, the count equals n^{n−2};
- the size, which is about 3.85·10⁷ at n = 12 (size^{1/n} = 4.29, falling toward 3).

## 4. Numerics (CHECKED NUMERICALLY)

Scripts are in `toy/r3/complete/`. All runners exit non-zero if a worker fails. Single-threaded BLAS, at most 5 jobs.

- **`kn_exact.py`** (runner `run_exact.sh`, logs `logs/exact-n*.jsonl`). Enumerates all n^{n−2} trees using `toy/r2/common.py` `enum_trees`, imported and not modified. It computes the following exactly:
  - β = ‖√μ‖_op;
  - B1, the averaged-conditioning bound (Lemma 1);
  - B2, the exact bit norms after conditioning;
  - B3, the proved product bound from Lemma 3.
  
  It asserts β ≤ B1 ≤ B2 ≤ B3 ≤ 1 and that the forbidden cell is never hit.
- **`kn_lemma.py`** (runner `run_lemma.sh`, logs `logs/lemma-n*.jsonl`). For n = 40…200, samples uniform trees by Prüfer code. For each sample it picks a random w ∈ W'_z, conditions on z and on the actual σ, and computes P[a_w=1 | z,σ] and P[b_w=1 | z,σ] exactly by matrix-tree determinants on the minor. It compares them with Lemma 3(i) and 3(ii).
- **`upper.py`**: Proposition U checks.
- **`constants.py`**: optimises the constants of Theorem K.

Exact β for UST(K_n). All rows have chain_ok = True. B1–B3 use r = 3.

| n | split | \|A\|/M | β | B1 | B2 | B3 | −ln β / n |
|---|---|---|---|---|---|---|---|
| 6 | random | 7/15 | 0.4408 | 0.5342 | 0.9050 | 0.9365 | 0.137 |
| 6 | star (B = edges at S) | 10/15 | 0.6944 | 0.6944 | 0.9044 | 0.9365 | 0.061 |
| 7 | local-search worst | 9/21 | 0.5240 | 0.5365 | 0.8828 | 0.9365 | 0.092 |
| 8 | random | 14/28 | 0.3438 | 0.3709 | 0.8843 | 0.9343 | 0.134 |
| 8 | local-search worst | 13/28 | 0.4293 | 0.4667 | 0.8841 | 0.9343 | 0.106 |
| 9 | random | 18/36 | 0.2772 | 0.2859 | 0.8906 | 0.9317 | 0.143 |
| 9 | clique (A = E(S)) | 15/36 | 0.3745 | 0.4062 | 0.8744 | 0.9337 | 0.109 |
| 9 | star | 21/36 | 0.4705 | 0.4705 | 0.8906 | 0.9317 | 0.084 |
| 9 | vertex-bipartite cut | 14/36 | 0.4335 | 0.4338 | 0.8906 | 0.9317 | 0.093 |

The smallest rate over everything tried comes from the star and vertex-bipartite splits: 0.061 per vertex for β at n = 6, and 0.084–0.096 for n = 7–9 (0.17–0.19 for β²). At n ≤ 9 the proved chain B3 is weak, since r = 3 and the constants are small. Its role is to confirm the chain of inequalities, not to show the rate.

Lemma 3 under real conditioning:
- There were about 2,000 logged conditional checks across n ∈ {40, 60, 100, 150, 200}, with random, clique, star, bipartite and skewed splits and r up to 12.
- There were **0 violations**.
- The minimum of (exact probability)/(bound) is 1.21 for the leaf bound and 1.025 for the deletion bound.
- The asymptotically exact value is P[a_w=1] ≈ θ_B·e^{−θ_A} (degree = 1 + Poisson(1)). The leaf bound θ_B/e is therefore loose by e^{θ_B} ∈ [1, e], which matches the observed ratios of 1.2–1.6.

## 5. For a reviewer: weakest steps and failure modes

1. **Lemma 3 under the conditioning σ. This is the step to check hardest.** The claim is that conditioning on the other isolated vertices' stars can only help. It rests on (α) and (β) in §3.3:
   - w's U-star is entirely undecided, and w is never merged with anything;
   - U–U edges are never decided, so contractions only short vertices (F3).
   
   If one tried to use non-isolated W-vertices (w with z-edges), (α) fails: w could already be attached, and the leaf event would no longer imply a_w = 1. That is why W'_z and the averaged Lemma 1 are needed.
2. **Lemma 1 is new (easy) and not in sec-tools.** `lem:chain` alone cannot absorb z (see §3.1).
3. **Adaptivity.** The coordinates (W'_z and its order) are fixed after conditioning on z, for each z separately, so `rem:adaptive` is respected. W itself depends on the split, which is allowed because β is bounded split by split.
4. **Constant.** κ ≈ 1/3600 per vertex, against a numerical ≈ 0.17. Most of the loss comes from r = n/24, the factor 1/2 from requiring |W'_z| > r/2, and p ≈ 0.013. A sharper version would choose W ⊆ Mid as large as Lemma 5 allows and use P[a_w=1] ≈ θ_B e^{−θ_A}. This is not needed for the 2^{Ω(n)} claim.
5. **What is not covered.**
   - Splits with |A| = o(M). The constant depends on the balance α (§3.5).
   - Size measured against a circuit whose leaves have arity above 2M/3; this is the same restriction as `cor:reduction`.
   - The exponent is Θ(n) = Θ(√M), not linear in the number of variables, and by Proposition U it cannot be.

## 6. Extensions and leads (CONJECTURED unless marked)

- **Dense expanders and other high-degree families** (open item 1 in `sec-open`, on degree dependence). The general leaf lemma (§3.3, PROVED) costs only e^{−ρd}, where ρ bounds resistances between w's neighbours inside G[U]. Suppose a family has (a) a Lemma-2-type supply of mixed vertices, (b) R_{G[U]}(u,u') = O(1/d) between neighbours, and (c) few tree edges inside W. Then the same proof should give β² ≤ e^{−Ω(n)} with a constant independent of the degree. Plausible examples are G(n, 1/2), K_{n,n} and powers of expanders.
  - Where it can fail: Lemma 2 used completeness (|E(L,H)| = |L||H|). Two expanders joined by a sparse cut, with A and B being the two sides, have no mixed vertices. That case needs R2-style cycle or expansion arguments instead.
- **R2's degree loss.** For bounded degree D with ρ = O(1/h²)-type resistance bounds, the leaf lemma gives per-vertex p ≥ (1/D)e^{−O(D·ρ)}. This might beat R2's D^{−4D/h}. Not checked. Worth passing to the R2/high-girth helper.
- **Rooted/directed variant.** For arborescences rooted at a fixed vertex, the out-edge (parent) variables form a Prüfer-like functional representation. The same side bits apply to the in-star of w. Not worked out.

## Files

- `notes/approaches/R3-complete.md` (this note)
- `toy/r3/complete/kn_exact.py`, `kn_lemma.py`, `upper.py`, `constants.py`, `run_exact.sh`, `run_lemma.sh`
- `toy/r3/complete/logs/` (`exact-n{6,7,8,9}.jsonl`, `lemma-n{40,60,100,150,200}.jsonl`, `*.status` = DONE)

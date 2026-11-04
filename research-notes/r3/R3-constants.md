# R3: better constants for the vertex-star UST bound

Direction: BETTER CONSTANTS for R2 (`notes/approaches/R2-highgirth.md`, Theorem B). Written 2026-11-04.
Reviewed (Sol 6.1, R3-constants-review.md): holds with fixes, 8/10, no fatal; fixes applied. The 1/h-power lower limit is conditional.
Numerics are in `toy/r3/constants/`, and all logs are in `toy/r3/constants/logs/`.

**Verdict: PARTIAL.**
* **Achieved and PROVED (fallback target).** The constant is now singly exponential in D, with no 1/h in the exponent. For every 2-connected G with max degree D ≥ 3 and edge expansion h, and every balanced split:

  β(UST(G); A|B)² ≤ exp(−c N),  with  c ≥ h₁^{5D+4} · (D+1)^{−4(D+1)} · e^{−10(D+1)},  h₁ := min(h, 1)   (Cor. P′).

  For 0 < h ≤ 1, R2 gives c ≈ D^{−4D/h}·(16D³/h²)^{−D}. In this regime, in nats, ln(1/c) drops from Θ((D/h)·log D) to O(D·log(D/h)). Optimising the same proof numerically gives a much better c (Theorem P and §6).
* **Achieved and PROVED (main target, in a regime).** If λ₂(L_G) > 1, then

  c ≥ λ₂ / (3D(D+λ₂/2)(D²+1)) · (D−1)/D² · (1 − 1/λ₂)^D   (Theorem S).

  This is polynomial in D, ≈ D^{−5}·e^{−D/(λ₂−1)}, so it is polynomial whenever λ₂ − 1 ≥ Ω(D). Near-Ramanujan D-regular graphs with D ≥ 6 qualify. For D = 6 and λ₂ = 1.53 it gives c ≈ 8·10⁻⁸, versus R2's e^{−115.6} ≈ 10^{−50}.
* **Main target in general: OPEN, with a conditional obstruction (Prop. L).** For the size-controlled, 2-connected hub family in §7.2, every edge at a mixed hub starts a bare path of length L. Every vertex-star chain certificate has c ≤ O(L^{−D/2}/(DL)), and the Theorem-A certificate has c ≤ O(L^{−D}/(DL)) (PROVED). If h ≥ c₀/L uniformly along this family (CONJECTURED), the degree in 1/h must grow linearly with D, ruling out c polynomial in D and 1/h for vertex-star certificates. Conditional on that expansion bound, Cor. P′ is within a factor of about 5 of the Theorem-A limit in that exponent. The best lead past it is replacing star bits by "status" bits (§7.3).

---

## 0. Statements and status

| # | Statement | Status |
|---|---|---|
| L2+ | **Telescoped key lemma.** For W ⊆ Z independent, H = G−W connected, κ = 1/λ₂(L_H), and every positive-probability σ on the stars of W∖w: Pr[T∩S_A(w)=∅ ∣ σ] ≥ (b/d)(1+κ)^{−a}, where a = \|S_A(w)\|, b = \|S_B(w)\| and d = a+b. Hence p_A p_B ≥ (ab/d²)(1+κ)^{−d} ≥ ((D−1)/D²)(1+κ)^{−D} | PROVED (§2); CHECKED NUMERICALLY over all σ, 65 cuts, 9 graphs |
| Ch | Cheeger–Mohar: λ₂(L_H) ≥ μ(Δ,i) := Δ − √(Δ²−i²) ≥ i²/(2Δ) for H with max degree Δ and edge expansion i | PROVED by citation (Mohar 1989) |
| L-H | W 3-separated (pairwise distance ≥ 3) ⇒ h(G−W) ≥ h − 1 | PROVED (§3.1) |
| L-S | W 3-separated ⇒ λ₂(G−W) ≥ λ₂(G) − 1 | PROVED (§3.1); CHECKED NUMERICALLY (N=300, D=6..16) |
| L-C | Number of connected S ∋ u with \|S\| = k and \|N(S)\| = j is ≤ C(k−1+j, j) | PROVED (§3.2), = Fomin–Villanger bound; CHECKED exhaustively on 8 graphs |
| L-P | **Peierls alteration.** If h′ ≤ h/(D+1) and S(q) ≤ 1/2 (series in §3.3), then there is an independent W* ⊆ Z with \|W*\| ≥ q\|Z\|/(2(D+1)) and h(G−W*) ≥ h′. Needs 2-connectivity | PROVED (§3.3) |
| T-P | **Theorem P.** c ≥ [2h/(3D(D+h))]·[q/(2(D+1))]·((D−1)/D²)·(1+1/μ(D,h′))^{−D}, for any admissible (h′, q) | PROVED |
| C-P′ | Closed form: q = (8DΛ³)^{−(D+1)}, Λ = e(1+1/h), h′ = h/(D+1) is admissible. Hence c ≥ h^{5D+4}(D+1)^{−4(D+1)}e^{−10(D+1)} for h ≤ 1 | PROVED (§4) |
| T-H | **Theorem H** (h > 1): c ≥ [2h/(3D(D+h))]·(D²+1)^{−1}·((D−1)/D²)·(1+1/μ(D,h−1))^{−D} | PROVED, but vacuous for large cubic and 4-regular graphs (§6) |
| T-S | **Theorem S** (λ₂ > 1): the formula above | PROVED |
| T-I | R2's packing and Lemma 4, with L2+ in place of L2: c ≥ [2h/(3D(D+h))]·\|B(s−1)\|^{−1}·((D−1)/D²)·(1+1/μ(D,h/(2D)))^{−D} | PROVED (uses R2's L4); never the best regime at realistic h |
| E | Circuit corollary: \|C\| ≥ (1−ε)²·e^{cN} for TV ≤ ε < 1. ε enters only through the prefactor | PROVED (paper cor:reduction) |
| L-a | Bare-path hubs: if each A-edge at w starts a path whose L−1 interior vertices have degree 2, then Pr[a_w = 1] ≤ L^{−a} | PROVED (§7) |
| L-b | For every vertex-star chain certificate, −ln(∏γ_k²) ≤ Σ_{w∈W∩Z} −ln(1−π_w), with π_w = min(Pr[a_w=1], Pr[b_w=1]). Any Theorem-A certificate has exponent ≤ 4ln2·Σ_w Pr[a_w=1]Pr[b_w=1] | PROVED (§7) |
| L-c | If the size-controlled hub family G_{L,n} has h ≥ c₀/L uniformly, vertex-star methods cannot give c = poly(1/h) of D-independent degree | CONDITIONAL; upper bound h ≤ 2/(L−1) and L-dependent caps PROVED; lower bound h ≳ 1/L CONJECTURED, not proved; λ₂L² ∈ [1.5, 2.4] CHECKED NUMERICALLY (L ≤ 16, tiny n), not a uniform bound |
| L-d | For the gadget split, the true β is e^{−Ω(n)}: the loss is the method's, not the truth's | CONJECTURED (heuristic in §7.3) |

---

## 1. Setup reused from R2

G is simple, 2-connected, with N ≥ 3 vertices, M edges, max degree D ≥ 3 and edge expansion

  h = min_{0<|S|≤N/2} e(S, V∖S)/|S|.

λ₂ = λ₂(L_G). T ~ UST(G). A balanced split is E = A ⊔ B with |A|, |B| ≥ M/3.

S_A(w) and S_B(w) are the A- and B-edges at w. N_A(w) and N_B(w) are the corresponding neighbours, with a = |S_A(w)| and b = |S_B(w)|. Vertex w is *mixed* (w ∈ Z) if a, b ≥ 1. The star bits are a_w = 1[T∩S_A(w)=∅] and b_w = 1[T∩S_B(w)=∅]; (a_w, b_w) = (1,1) is impossible.

Used as is from R2 (all reviewed in `R2-review.md`):
* **L1 (spatial Markov).** Given T∩F (F any edge set), T∖F is the UST of (V, E∖F)/(T∩F).
* **K2.** Pr[e∉T] = 1/(1+R_{K−e}(e)), where R is effective resistance with unit conductances.
* **R2 (Rayleigh).** Identifying vertices does not increase effective resistance.
* **L3.** |Z| ≥ 2hM/(3D(D+h)) ≥ 2hN/(3D(D+h)). The last step uses M ≥ N, since min degree is ≥ 2 under 2-connectivity.
* **Theorem A.** If W ⊆ Z is independent and, for every w ∈ W and every positive-probability configuration σ of T on ∪_{x∈W∖w} S(x), Pr[a_w=1∣σ] ≥ p_A(w) and Pr[b_w=1∣σ] ≥ p_B(w), then

  β² ≤ ∏_{w∈W} (1+√(1−4p_A p_B))/2 ≤ exp(−Σ_{w∈W} p_A(w)p_B(w)).

  The chain rule conditions on the past bits, which are coarser than σ, so lower bounds that hold for every σ survive. The last step is lem:two-by-two.

Hence, for any W as in Theorem A with a uniform bound p_A p_B ≥ x:

  **c ≥ (|W|/N)·x.**   (★)

Every regime below is an instance of (★) with a particular W and x.

## 2. Lemma 2+ (telescoped key lemma)

**Lemma 2+.** Let W ⊆ Z be independent, with H := G − W connected, and put κ := 1/λ₂(L_H). Let w ∈ W have A-edges e_1, …, e_a, with e_j = w u_j, and b ≥ 1 B-edges, d = a+b. Then for every configuration σ of T on F := ∪_{x∈W∖w} S(x) with Pr[σ] > 0,

  Pr[T∩S_A(w)=∅ ∣ σ] ≥ (b/d)·(1+κ)^{−a},

and symmetrically Pr[T∩S_B(w)=∅ ∣ σ] ≥ (a/d)(1+κ)^{−b}. So p_A p_B ≥ (ab/d²)(1+κ)^{−d} ≥ ((D−1)/D²)(1+κ)^{−D}.

*Proof.* By L1, conditioned on σ, T∖F is the UST of K := (V, E∖F)/σ_in. Let K_j := K − {e_1,…,e_j}. Conditioning on {e_1,…,e_{j−1} ∉ T} is again a UST, this time of K_{j−1}. Then K2 gives

  Pr[T∩S_A(w)=∅ ∣ σ] = ∏_{j=1}^a Pr[e_j∉T ∣ σ, e_1..e_{j−1}∉T] = ∏_{j=1}^a 1/(1+R_{K_j}(u_j,w)).

*Comparison graph.* Let G_j := (V, E∖F) − {e_1,…,e_j}. In G_j every x ∈ W∖w is isolated, and K_j is obtained from G_j by identifying each such x with its σ-neighbours. Those neighbours are not w, by independence. By R2, R_{K_j}(u_j,w) ≤ R_{G_j}(u_j,w). The component of w in G_j is H together with w and the m_j := d − j ≥ b ≥ 1 remaining edges from w to the set Y_j ⊆ V(H) of remaining neighbours. Note u_j ∉ Y_j because G is simple.

*Spreading flow.* Send a unit flow from u_j to the uniform measure on Y_j inside H, then 1/m_j along each edge yw (y ∈ Y_j). By Thomson's principle,

  R_{G_j}(u_j,w) ≤ χ_jᵀ L_H⁺ χ_j + m_j·(1/m_j)² ≤ κ‖χ_j‖² + 1/m_j = κ(1+1/m_j) + 1/m_j,

with χ_j = δ_{u_j} − (1/m_j)1_{Y_j} ⊥ 1. So 1 + R ≤ (1+κ)(1+1/m_j) = (1+κ)(m_j+1)/m_j.

*Telescoping.* ∏_{j=1}^a m_j/(m_j+1) = ∏_{j=1}^a (d−j)/(d−j+1) = (d−a)/d = b/d.

*Final bound.* Since a, b ≥ 1, ab ≥ d−1. The function (d−1)/d² is decreasing for d ≥ 2, and d ≤ D. ∎

**What is new compared with R2's Lemma 2.** R2 bounds each factor by 1/(1+R_u) and routes through a single neighbour y, so 1 + R_u ≤ 2 + R_H(u,y) ≤ 2(1+κ), since ‖δ_u−δ_y‖² = 2. Spreading the flow over all remaining neighbours replaces this by (1+κ)(m+1)/m. That is better by 2m/(m+1) ∈ [1, 2) per factor, and the 1/m losses telescope to b/d instead of compounding to 2^{−a}. This is the "conditional DPP kernel" idea made concrete: the product is the sequential Schur-complement expansion of Pr[T avoids S_A(w)] = det(I − K_{S_A(w)}) for the conditioned transfer-current kernel K, and each factor is bounded by one explicit flow.

**Numerics (CHECKED NUMERICALLY).** `toy/r3/constants/exact_lemma2plus.py` enumerates all spanning trees of small graphs. For each w it computes the exact minimum over σ on the other stars of W of Pr[a_w=1∣σ]. It compares that minimum with (i) the exact-energy flow bound ∏ 1/(1+E_flow) and (ii) the closed form (b/d)(1+1/λ₂(H))^{−a}. It also checks β ≤ ‖√μ‖ ≤ chain ≤ Theorem-A bound with the closed form. Results (`run_exact.sh`, logs `exact-*.jsonl`):

* 65 balanced cuts on 9 graphs (K33, cube, Petersen, Heawood, Möbius–Kantor, Pappus, rr3_14, rr3_16, rr4_12), with τ up to 1.26·10⁶, W both independent and 3-separated: **0 violations**.
* min(exact/flow) = 1.0000. The flow bound is attained for some σ, so the spreading flow is not the slack.
* min(exact/closed) = 1.000 on K33 and ≥ 1.07 elsewhere.
* Every chain ≤ bound check passed.
* Per instance, R2's exact-resistance product can exceed the closed form (e.g. Heawood with sep 2). The gain is at the level of the worst-case closed forms, which is what Theorem B uses.

## 3. Deleting W while keeping expansion

**Ch (cited).** For a connected simple graph H with ≥ 2 vertices, max degree Δ and edge expansion i,

  λ₂ ≥ μ(Δ,i) := Δ − √(Δ²−i²) ≥ i²/(2Δ).

For ≥ 3 vertices, if λ₂ > Δ this is automatic, since i ≤ Δ and μ(Δ,i) ≤ Δ. If λ₂ ≤ Δ, the Cheeger–Mohar argument gives i² ≤ λ₂(2Δ−λ₂); solving the quadratic gives the stated lower bound. The upper inequality is not asserted when λ₂ > Δ (K₃ is a counterexample). For two vertices, H = K₂, Δ = i = 1 and λ₂ = 2 ≥ μ(1,1) = 1.

See Mohar, *Isoperimetric numbers of graphs*, JCTB 47 (1989) 274–291 (https://www.sciencedirect.com/science/article/pii/0095895689900294), and Thm 2.9 of Mohar–Poljak 1993 (https://people.orie.cornell.edu/dpw/orie6334/Fall2016/MoharPoljak93.pdf), with the above case distinction rather than the overbroad upper inequality printed there.

μ is decreasing in Δ, so Δ may be replaced by D. The optimized constants retain μ(D,h′); Cor. P′ uses only the weaker λ₂ ≥ i²/(2D).

### 3.1 Dense packings when h > 1 or λ₂ > 1

W is *3-separated* if distinct w, w′ ∈ W have distance ≥ 3: W is independent and no vertex has two W-neighbours. Greedy selection in Z gives |W| ≥ |Z|/(D²+1), since |B₂(w)| ≤ 1 + D + D(D−1).

**Lemma H.** W 3-separated ⇒ h(G−W) ≥ h−1.

*Proof.* Let S ⊆ V∖W with |S| ≤ (N−|W|)/2 ≤ N/2. Each v ∈ S has ≤ 1 neighbour in W, so e(S,W) ≤ |S| and e_{G−W}(S, V∖W∖S) = e_G(S, V∖S) − e(S,W) ≥ (h−1)|S|. ∎

**Lemma S.** W 3-separated ⇒ λ₂(L_{G−W}) ≥ λ₂(L_G) − 1.

*Proof.* Take f on V∖W with f ⊥ 1 and ‖f‖ = 1, and extend it by 0 on W. The extension g satisfies g ⊥ 1_V and ‖g‖ = 1. Then

  gᵀL_G g = fᵀL_{G−W} f + Σ_{v∉W} f_v²·|N(v)∩W| ≤ fᵀL_{G−W} f + 1,

so fᵀL_{G−W} f ≥ λ₂(L_G) − 1. ∎

The same proof gives λ₂(G−W) ≥ λ₂ − t for independent W in which every vertex has ≤ t W-neighbours. This could shave the 1/(D²+1) packing loss toward 1/D when λ₂ ≫ 1. That variant is not used or worked out here.

### 3.2 Counting connected sets

**Lemma C.** In any graph, the number of connected S ∋ u with |S| = k and outer vertex boundary |N(S)| = j is ≤ C(k−1+j, j).

*Proof.* Fix orders on V and on adjacency lists. Explore S from u in BFS fashion. Keep a set of discovered vertices (initially {u}) and a set of known-outside vertices. When processing a discovered vertex v, examine each neighbour x of v in order. If x is already discovered or known outside, skip it without writing anything. Otherwise write the bit 1[x ∈ S]: on 1, mark x discovered and queue it; on 0, mark it outside.

Every vertex of S∖u yields exactly one 1. Every vertex of N(S) yields exactly one 0, since it is adjacent to a processed vertex of S. Given G, the string determines S, because the decoder replays the exploration. So the strings have k−1 ones and j zeros, and S ↦ string is injective. ∎

This is the bound of Fomin–Villanger, *Treewidth computation and extremal combinatorics*, Combinatorica 2012 (arXiv 0803.1321), restated with vertex boundary.

**Corollary.** For b ≥ 1, K := ⌊b/h⌋ and Λ := e(1+1/h), the number of connected S ∋ u with |S| ≤ N/2 and |∂_E S| = b is at most

  F(b) := Σ_{k=1}^{K} Σ_{j=0}^{b} C(k−1+j, j) = C(K+b+1, b+1) − 1 ≤ Λ^{b+1}.

*Proof.* Such S has |S| ≤ b/h, by expansion, and |N(S)| ≤ b. Apply the hockey-stick identity twice. For the last step, C(K+b+1, b+1) ≤ (e(K+b+1)/(b+1))^{b+1} ≤ (e(1+1/h))^{b+1}.

### 3.3 Lemma P (Peierls alteration)

**Lemma P.** Let G be 2-connected with N ≥ 3. Fix 0 < h′ ≤ h/(D+1) and q ∈ (0,1], and put

  θ := (1−h′/h)/D,   m(b) := max(2, ⌊θb⌋ + 1),
  S(q) := D·Σ_{b≥2} F(b)·C(b−1, m(b)−1)·q^{m(b)−1}.

If S(q) ≤ 1/2, then there is an independent W* ⊆ Z with |W*| ≥ q|Z|/(2(D+1)) and h(G − W*) ≥ h′. In particular G − W* is connected.

*Proof.*

1. *Random subset.* Greedily pick an independent W₀ ⊆ Z with |W₀| ≥ |Z|/(D+1). Let R ⊆ W₀ include each vertex independently with probability q.
2. *Dangerous sets.* Call a pair (w, S) dangerous if w ∈ R, S is connected in G, w ∉ S, w is adjacent to S, |S| ≤ N/2, and |N(S)∩R| ≥ m(b) where b = |∂_E S|. S may intersect R.
3. *Alteration.* Let R* := R minus every w that lies in some dangerous pair.
4. *Expectation.* Fix w ∈ W₀ and condition on w ∈ R; the other inclusions are unaffected. Choose a neighbour u ∈ S of w (≤ D ways). By the Corollary there are ≤ F(b) sets S ∋ u with boundary b. Other than w, S has ≤ b−1 boundary neighbours, and at least m(b)−1 of them must lie in R. By a union bound over (m(b)−1)-subsets this has probability ≤ C(b−1, m(b)−1)q^{m(b)−1}. Also b ≥ 2, because 2-connected implies 2-edge-connected. So Pr[w removed ∣ w ∈ R] ≤ S(q) ≤ 1/2 and E|R*| ≥ q|W₀|/2. Fix an outcome with |R*| ≥ q|W₀|/2, and put W* := R*. It is independent and contained in Z.
5. *Expansion of H* := G − W*.* Suppose some S ⊆ V(H*) with |S| ≤ |V(H*)|/2 has e_{H*}(S) < h′|S|. Some component of H*[S] has the same property, because the boundaries of the components add up. So assume S is connected. Let b = |∂_E^G S|. Then |S| ≤ N/2 gives |S| ≤ b/h, and e_{H*}(S) = b − e(S,W*).
   * *Many W*-neighbours.* e(S,W*) > b − h′|S| ≥ b(1−h′/h). Each W* vertex has ≤ D edges into S, so |N(S)∩W*| ≥ ⌊θb⌋ + 1.
   * *Zero W*-neighbours.* Then e_{H*}(S) = b ≥ h|S|, a contradiction.
   * *One W*-neighbour w.* Then e_{H*}(S) ≥ b − D ≥ h|S| − D. Also e_{H*}(S) ≥ 1: otherwise w separates S from V∖(S∪w), which is nonempty since |S| ≤ N/2 and N ≥ 3, contradicting 2-connectivity. If |S| ≤ 1/h′, then 1 ≥ h′|S|. If |S| > 1/h′ ≥ (D+1)/h, then h|S| − D − h′|S| ≥ (D+1)(1−h′/h) − D ≥ 0, using h′ ≤ h/(D+1). Either way this is a contradiction, so |N(S)∩W*| ≥ m(b).
   * *Contradiction.* Pick w ∈ N(S)∩W*. Then w ∈ R, w ∉ S, and |N(S)∩R| ≥ |N(S)∩W*| ≥ m(b). So (w, S) is dangerous and w ∉ R*, which contradicts w ∈ W*.

   Hence h(H*) ≥ h′. ∎

*Where 2-connectivity enters.* It gives b ≥ 2 and the one-neighbour case. This matches R2's R1: expansion alone is not enough.

## 4. The new constants

Each theorem below combines (★), L3 and Lemma 2+, with κ bounded by Ch, Lemma H or Lemma S.

**Theorem P (any 2-connected G).** For every 0 < h′ ≤ h/(D+1) and q with S(q) ≤ 1/2,

  c ≥ c_P := [2h/(3D(D+h))] · [q/(2(D+1))] · ((D−1)/D²) · (1 + 1/μ(D,h′))^{−D}.

*Proof.* Take W = W* from Lemma P in (★). Lemma 2+ applies with κ = 1/λ₂(G−W*) ≤ 1/μ(D,h′) by Ch, since G−W* has max degree ≤ D. ∎

**Corollary P′ (closed form).** Take h′ = h/(D+1), so θ = 1/(D+1), and q₀ := (8DΛ³)^{−(D+1)} with Λ = e(1+1/h). Then S(q₀) ≤ 1/2, and for 0 < h ≤ 1 and D ≥ 3,

  ln(1/c) ≤ (5D+4)·ln(1/h) + 4(D+1)·ln(D+1) + 10(D+1).

For h > 1, apply this with h replaced by 1.

*Proof.*

1. *Bounding m(b).* m(b)−1 = max(1, ⌊θb⌋) ≥ θb/2. Check the two cases θb < 2 and θb ≥ 2. So q^{m−1} ≤ q^{θb/2}, and C(b−1, m−1) ≤ 2^{b−1}.
2. *Geometric series.* Each term is ≤ Λ^{b+1}2^{b−1}q^{θb/2} = (Λ/2)x^b, with x := 2Λq^{θ/2}. Hence S ≤ (DΛ/2)·x²/(1−x) ≤ DΛx² when x ≤ 1/2. At q₀, q₀^θ = 1/(8DΛ³), so x² = 4Λ²/(8DΛ³) = 1/(2DΛ). That gives S ≤ 1/2, and also x ≤ 1/4.
3. *Size of κ.* Ch gives κ ≤ 2D/h′² = 2D(D+1)²/h². For h ≤ 1, 1 + 2D(D+1)²/h² ≤ 3D(D+1)²/h².
4. *Collecting the logarithm.* The coefficient of ln(1/h) is 1 (from L3), plus 3(D+1) (from Λ³ in q₀), plus 2D (from κ^D). The rest is ln(3D(D+1)/2) + ln(2(D+1)) + (D+1)(ln 8D + 3 + 3ln2) + ln(D²/(D−1)) + D·ln(3D(D+1)²). For D = 3, the unrounded rest is 54.41332 or less, below 16ln4 + 40 = 62.18070 or more. For D ≥ 4, it is ≤ 4(D+1)ln(D+1) + 4ln(D+1) + 8.3(D+1) + 1.3 ≤ 4(D+1)ln(D+1) + 10(D+1), using 4ln(D+1) ≤ 1.39(D+1) and 1.3 ≤ 0.26(D+1).

The inequality was also checked on a grid with D ∈ [3,200] and h ∈ (0,1]: the worst slack is −8.98 nats. ∎

So c ≥ h^{5D+4}·(D+1)^{−4(D+1)}·e^{−10(D+1)}. R2 has c = D^{−(4D/h+O(1))}·(2+16D³/h²)^{−D}.

**Theorem H (h > 1).** Take W 3-separated and greedy in Z, with κ = 1/μ(D, h−1) from Lemma H and Ch:

  c ≥ [2h/(3D(D+h))]·(D²+1)^{−1}·((D−1)/D²)·(1+1/μ(D,h−1))^{−D}.

**Theorem S (λ₂ > 1).** Take W 3-separated and greedy in Z, with κ = 1/(λ₂−1) from Lemma S. L3 is used with h ≥ λ₂/2, the easy Cheeger direction: e(S,S^c) ≥ λ₂|S||S^c|/N:

  c ≥ [λ₂/(3D(D+λ₂/2))]·(D²+1)^{−1}·((D−1)/D²)·(1−1/λ₂)^D.

If h > λ₂/2 is known, h can be used in the first factor.

**Theorem I.** This is R2's s-packing and Lemma 4 with Lemma 2+, as in the table in §0. It is never the best regime at realistic h. The final constant therefore does **not** depend on R2's Lemma 4.

**Best constant.** c_new(D,h[,λ₂]) = max(c_P, c_H, c_S, c_I). c_P is evaluated with the exact series and optimised over (h′, q); the other three are closed forms.

*Evaluation of c_P.* `toy/r3/constants/constants_table.py` computes ln S(q) in log space. It sums b ≤ 4000 exactly, using C(K+b+1, b+1) ≥ F(b). It bounds the tail b > 4000 by the geometric series Λq^{−1}r^{4001}/(1−r), with r = 2Λq^θ, using F ≤ Λ^{b+1}, C(b−1, m−1) ≤ 2^b and m−1 ≥ θb−1. It then bisects on ln q (60 steps) for each h′ = (h/(D+1))·i/12, i = 1..12. This is a floating-point evaluation, not an interval-arithmetic certificate: the reported lo satisfies S ≤ 1/2 as computed, but rounding error at the boundary is not certified. The reviewer independently reproduced the rounded table entries with 65-digit arithmetic. Cor. P′ is the proved numerical-free guarantee.

## 5. Dependence on ε (δ)

cor:reduction gives |C| ≥ δ²/max_{A|B}β² with δ = 1 − TV(P, Q_C). So

  |C| ≥ (1−ε)²·e^{c_new·N}   for TV ≤ ε < 1.

The constant c does not depend on ε. The bound stays exponential if (1−ε)² ≥ exp(−(c−c′)N) for some fixed 0 < c′ ≤ c, giving |C| ≥ exp(c′N). At 1−ε = exp(−cN/2) the bound is only 1. For fixed ε < 1 there is no issue. Nothing in R2 or R3 loses more in ε.

## 6. Old vs new: c(D,h)

From `logs/constants_table.txt`, floating-point evaluated and independently reproduced to 65 digits by the reviewer, not interval-arithmetic certified. Entries are ln(1/c), so lower is better. "Best new" is max over regimes.

**Certified bounds on h, not attainable ranges.** For large N, cubic graphs have h ≤ 1/3 + o(1) and 4-regular graphs have h ≤ 4/5 + o(1). Both follow from h ≤ 2·bw/N and the Monien–Preis bisection bounds bw ≤ (1/6+o(1))N and (2/5+o(1))N (J. Discrete Algorithms 4 (2006) 475–498, https://doi.org/10.1016/j.jda.2005.12.009). Rows above these ceilings are omitted below and kept in the log, marked.

For random cubic graphs, the certified asymptotic lower/upper bounds are approximately 0.086 and 0.2796 (up to o(1)): the lower bound uses Cheeger with λ₂ ≈ 3 − 2√2; the upper uses Lichev–Mitsche (2023), bw ≤ 0.139822N + o(N). For random 4-regular graphs the corresponding certified bounds are approximately 0.27 and 0.8; Ramanujan 6-regular graphs have h ≥ 0.76. These are bounds, not claims that all intervening values are attained.

**New, unreviewed constraint.** Gaspers, arXiv:2611.05146 (2026-11-03), claims bisection bounds giving expansion ceilings ≈ 0.2796 for maximum degree 3 and 2/3 for maximum degree 4. If correct, both rows marked † below are asymptotically vacuous for arbitrary such families. The D = 3, h = 0.3 row is already asymptotically vacuous for random cubic graphs by Lichev–Mitsche.

| D | h | R2 stated: ln(1/c) | c_R2 | new: ln(1/c) | c_new | regime |
|---|---|---|---|---|---|---|
| 3 | 0.05 | 306.6 | 10^−133 | 73.6 | 10^−32.0 | P (h′=0.0125, ln q=−32.8) |
| 3 | 0.1 | 169.9 | 10^−73.8 | 63.5 | 10^−27.6 | P |
| 3 | 0.2 | 99.2 | 10^−43.1 | 53.6 | 10^−23.3 | P |
| 3 | 0.3 † | 74.4 | 10^−32.3 | 48.0 | 10^−20.8 | P |
| 4 | 0.1 | 274.8 | 10^−119 | 83.8 | 10^−36.4 | P |
| 4 | 0.2 | 157.7 | 10^−68.5 | 71.2 | 10^−30.9 | P |
| 4 | 0.4 | 96.1 | 10^−41.7 | 59.0 | 10^−25.6 | P |
| 4 | 0.67 † | 69.3 | 10^−30.1 | 50.4 | 10^−21.9 | P |
| 6 | 0.2 | 290.6 | 10^−126 | 107.6 | 10^−46.7 | P |
| 6 | 0.5 | 149.8 | 10^−65.1 | 84.9 | 10^−36.9 | P |
| 6 | 0.76 | 115.7 | 10^−50.2 | 75.1 | 10^−32.6 | P |
| 6 | 1.0 | 97.8 | 10^−42.5 | 69.1 | 10^−30.0 | P |
| 6 | λ₂=1.53 (Ramanujan), h=λ₂/2 | 115.6 | 10^−50.2 | **16.3** | **8.0·10^−8** | S |
| 6 | λ₂=3 | 78.3 | 10^−34 | 11.8 | 7.3·10^−6 | S |
| 8 | λ₂=2.71 (Ramanujan) | 123.6 | 10^−53.7 | 14.5 | 5.1·10^−7 | S |
| 10 | λ₂=4.00 (Ramanujan) | 135.8 | 10^−59.0 | 14.4 | 5.6·10^−7 | S |
| 16 | λ₂=8.25 (Ramanujan) | 184.0 | 10^−79.9 | 15.2 | 2.5·10^−7 | S |
| 32 | λ₂=20.9 (Ramanujan) | 328.6 | 10^−143 | 17.3 | 3.1·10^−8 | S |
| 64 | λ₂=48.1 (Ramanujan) | 628.8 | 10^−273 | 19.7 | 2.8·10^−9 | S |

Further columns in the log:
* R2 + Lemma 2+ alone (Theorem I) gains 10–100 nats.
* The closed form P′ is 25–40 nats worse than the optimised P.
* Theorem H wins whenever h > 1, but no D ≤ 4 graph family has h > 1 for large N. For D = 6, I do not know an infinite family with h ≥ 1.5; Alon (1997) gives h ≤ D/2 − c√D.

*Shape.* In regime P at D = 3, ln(1/c) grows by ≈ 10 nats per halving of h. That is polynomial in h of degree ≈ 14.6, against 5D+4 = 19 in the closed form. R2 grows ∝ 1/h.

*Medium-graph sanity check (CHECKED NUMERICALLY).* Script `spectral_regime.py`, logs `spec-*.jsonl`. On random D-regular graphs with N = 300 and D = 6, 8, 10, 16, we used random, Fiedler-bisection and BFS-ball splits with greedy 3-separated W.
* λ₂(G−W) ≥ λ₂(G) − 1 held in all 12 cases. The observed drop was only 0.24–0.41.
* Wilson-sampled σ with exact conditional probabilities gave 0 violations of Lemma 2+ in 480 checks (240 samples × both sides), with min ratio 1.05.
* The per-split certificate with κ = 1/(λ₂−1) is 3·10⁻⁵ to 7·10⁻⁴ per N, against the proved c_S of 3·10⁻⁷ to 3.5·10⁻⁶.

## 7. Limits of vertex-star certificates and the next lead

### 7.1 Bare paths (PROVED)

**Lemma L-a.** Suppose each A-edge e_j = wu_{j,1} at w starts a path w, u_{j,1}, …, u_{j,L} whose interior vertices u_{j,1..L−1} have degree 2 in G. Then Pr[a_w = 1] ≤ L^{−a}.

*Proof.* As in §2, the probability is ∏_j 1/(1+R_{G−{e_1..e_j}}(u_{j,1}, w)). After e_j is deleted, any unit flow from u_{j,1} to w crosses the L−1 path edges in series at full strength, so R ≥ L−1 and each factor is ≤ 1/L. ∎

**Lemma L-b.** Take any ordered W and the chain rule (lem:chain) with bits (a_w, b_w). Each conditional 2×2 law has a zero (1,1)-cell. Its √-matrix has a first row of norm √(1−Pr[a=1∣past]) and a first column of norm √(1−Pr[b=1∣past]). Some past has Pr[a_k=1∣past] ≤ Pr[a_k=1], and likewise for b. Hence γ_k² ≥ 1 − π_k, with π_k = min(Pr[a_k=1], Pr[b_k=1]).

For Theorem A, any admissible p_A ≤ inf_σ Pr[a_w=1∣σ] ≤ Pr[a_w=1]. The function f(x) = −ln((1+√(1−4x))/2) is convex with f(1/4) = ln 2, so the certified exponent is ≤ 4ln2·Σ_w Pr[a_w=1]Pr[b_w=1]. ∎

### 7.2 The hub gadget

G_{L,n} has two 2-connected cubic expander cores X_A and X_B, each of size Θ(max(a,b)·n), with fixed size-comparison constants, and n ≥ 2 hubs. Each hub has a ≥ 2 A-paths and b ≥ 2 B-paths of length L ≥ 2 to globally distinct attachment vertices in the respective cores. The split is A = E(X_A) ∪ A-paths. It is balanced for a = b. It is 2-connected, Z = {hubs}, max degree max(a+b, 4), and N ≥ n(a+b)(L−1).

By L-a and L-b:
* every vertex-star chain certificate has c ≤ 2L^{−min(a,b)}/((a+b)(L−1));
* the Theorem-A certificate has c ≤ 2.78·L^{−(a+b)}/((a+b)(L−1)).

With a = b = D/2 these are O(L^{−D/2}/(DL)) and O(L^{−D}/(DL)).

On expansion, h ≤ 2/(L−1) (PROVED: the interior of one path). If h ≥ c₀/L, as expected for subdivided expanders (CONJECTURED, not proved here), then c_method ≤ O(h^{D/2+1}/D). Alternatively, a uniform λ₂ ≳ L^{−2} along the family would give c_method ≤ O(h^{D/4}) via Cheeger. This spectral lower bound is also unproved; the finite observations below do not establish it.

Thus **CONDITIONAL on the uniform expansion lower bound, no vertex-star certificate gives c = poly(D, 1/h) with a D-independent degree.** Under that same condition, Cor. P′, at h^{5D+4}, is within a factor (5D+4)/(D+1) ≈ 5 in that exponent of the Theorem-A limit h^{D+1}. The optimised Theorem P at D = 3 numerically behaves like h^{≈14.6}, a factor ≈ 3.7 under the same conditional comparison.

Numerics (`hub_gadget.py`, logs `hub-*.jsonl`): exact Kirchhoff probabilities for a = b = 2 and a = b = 3, L up to 16.
* Max Pr[a_w=1] is 0.29–0.40 of L^{−a} for a=2, and 0.22–0.35 of L^{−a} for a=3. The ratio slowly increases with L, consistent with L-a.
* λ₂·L² ranges over 1.53–2.33.
* The chain-cap per N falls from 3.2·10⁻² at L=2 to 9.5·10⁻⁵ at L=16 for a=b=2.

### 7.3 Why this is the method's loss, not the truth's (CONJECTURED), and the lead

In G_{L,n}, T∩A and T∩B interact only through the hubs. Each hub has an A-*status*: is it joined to X_A by a full A-path in T? It also has a B-status. Path interiors are monochromatic, so (not joined, not joined) is impossible. When both cores are internally connected in T, (joined, joined) can occur at only one hub without creating a cycle.

Both statuses have probability Θ(1), so the status pairs are strongly anti-correlated with Θ(1) information per hub. This suggests β ≤ e^{−Ω(n)} = e^{−Ω(N/(DL))}, which would be polynomial in h under the conjectured h ≳ 1/L. Star bits only see the L^{−a}-rare event "no A-edge at w".

**Best lead for the main target.**
* Replace a_w by a deep status bit: the T∩A-component of w leaves B_r(w), with r ≈ the local resistance scale.
* The (1,1)-impossibility then fails in general. It needs a 2×2 lemma for a *rare* cell, Pr[(1,1) ∣ past] ≤ η, in place of lem:two-by-two.
* It also needs the analogue of Lemma 2+ for status events, which are increasing and decreasing events of a DPP restricted to a ball.

Not attempted within the time box.

## 8. Weakest steps for a reviewer

1. **Lemma P.** Check three things.
   * Dangerous sets must be allowed to meet R∖R*; the proof does allow it.
   * The one-neighbour case uses 2-connectivity, N ≥ 3 and h′ ≤ h/(D+1).
   * The union bound C(b−1, m−1)q^{m−1} uses the independence of inclusions given w ∈ R.

   The numerical S(q) is floating-point log-space with an explicit tail. Cor. P′ avoids numerics entirely.
2. **Lemma 2+ comparison graph.** K_j is G_j with identifications, which needs W independent so that no identification touches w. Also u_j ∉ Y_j. The exact checks attain ratio 1.000, so the bound is tight in places.
3. **Uniformity over σ versus chain conditioning.** Theorem A conditions on past bits. Lemma 2+ is proved for every σ on all of W∖w, which is a finer conditioning, so averaging is fine.
4. **Vacuous rows.** Any comparison at h above the Monien–Preis limits says nothing about large graphs; §6 also flags the random-cubic ceiling and the stronger claims of the new, unreviewed Gaspers preprint. Theorem H is only meaningful for D ≳ 6–8.
5. **Prop. L(c).** The expansion lower bound of the gadget is not proved. Only the qualitative "degree grows with D" conclusion depends on it, and the Cheeger route needs λ₂ ≳ L^{−2} for large n, also unproved.

## 9. Literature and novelty (searched 2026-11-04, deep web search)

* No matching result on 2024–2026 approximate-PC lower bounds for UST/DPPs or constants of this kind was found in the stated search.
  * arXiv 2506.01281 (NeurIPS 2025) proves general hardness and separations, not UST.
  * UAI 2020, arXiv 2006.15233 (Honghua Zhang, Steven Holtzen, Guy Van den Broeck), proves exact smooth-SPN lower bounds for UST(K_n), with size ≥ 2^{n/30240}.
  * arXiv 2109.06941 covers monotone circuits for the spanning-tree polynomial on expanders.
  * Martens–Medabalimi, arXiv 1411.7717.
* **Lemma 2+.** Elementary: a spreading flow plus telescoping. Not found as such. It is the explicit-flow form of the sequential Schur expansion of det(I − K_S).
* **Lemma P.** Standard Peierls/alteration technique, with the counting bound of Fomin–Villanger (arXiv 0803.1321). Related is Sonny Ben-Shimon and Michael Krivelevich, *Vertex percolation on expanders* (arXiv 0710.2296), which covers random deletion at rate n^{−α}, with high-probability statements and no explicit constants. The deterministic explicit-constant statement used here was not found.
* **Mohar 1989** gives the Cheeger form used.
* **Context.** Mészáros–Peled, arXiv 2610.11554 (Oct 2026), couple UST with a random 1-out subgraph with O(n log k/√k) expected disagreement at min degree k. Per-vertex dependence fades with degree, consistent with c ~ D^{−5} in Theorem S. It is not a TV or PC statement.

## 10. Files

* `toy/r3/constants/constants_table.py`: all closed forms and the optimised S(q). Output in `logs/constants_table.txt`.
* `toy/r3/constants/exact_lemma2plus.py` with `run_exact.sh` (≤ 6 jobs, loud failures): exact Lemma 2+, Lemma S and Theorem-A chain checks. Logs `logs/exact-*.jsonl`.
* `toy/r3/constants/spectral_regime.py`: N = 300 spectral-regime checks. `hub_gadget.py`: the limit gadget. Both are run by `run_medium.sh`; logs `logs/spec-*.jsonl` and `logs/hub-*.jsonl`.
* `toy/r3/constants/count_check.py`: exhaustive check of Lemma C on 8 graphs (max count/bound = 1.000). Log `logs/count_check.txt`.

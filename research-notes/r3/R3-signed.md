# R3: nonnegativity vs decomposability for the uniform spanning tree

Written Nov 4 2026, about 01:30–03:15 CST. Scope: direction R3. Can small *decomposable* circuits that may *subtract* represent UST(G) on expanders, exactly or approximately? The v1 paper (`notes/paper-ust/`) was read and not modified. Scripts and logs are in `toy/r3/signed/`.

Labels: **[PROVED]** means a full proof is in this note, or a proof that reduces to a cited v1 lemma or a cited paper. **[CHECKED]** means exact small-graph numerics (§9). **[CONJ]** means conjectured. **[OPEN]** means open.

---

## Verdict

Reviewed (Sol 6.1, R3-signed-review.md): A, B, C hold; S, K, M, X hold with fixes; 8/10; fixes applied.

**PARTIAL.** For the structured-decomposable classes specified below, the answer is decomposability, not nonnegativity. For unstructured signed decomposable circuits, the question is open and equivalent to a known hard problem.

1. **Structured signed circuits, exact and approximate: outcome (b) holds [PROVED].** On linked expanders (the hypotheses of v1 `thm:general`), take any smooth structured-decomposable circuit with real weights and any vtree. If it computes Q with ½‖UST − Q‖₁ ≤ ε < ½, it has size 2^{Ω(N)}. The exponent matches v1's nonnegative bound up to a constant factor (Corollary C). Subtraction buys nothing once the decomposition is structured. The threshold ½ is sharp for unnormalised Q, since Q = 0 achieves ε = ½.
2. **Normalized squares/PSD forms of structured-decomposable amplitudes obey the same lower bound [PROVED].** An explicit sum of squares of smooth structured-decomposable signed circuits, each with its own vtree, charged total size S = Σ_k s_k, satisfies S ≥ (1 − TV)²/β_*². The shared structured-output PSD model is also covered (Theorem S). This matches the nonnegative lower-bound expression, not necessarily optimal sizes; inner decomposability is essential.
3. **Positive route [PROVED, CHECKED].** There is an explicit smooth structured-decomposable signed PC with indicator leaves for UST, of size poly(N)·2^{O(w)}, where w is the edge-order width (Theorem B). It is the determinant (Laplace-expansion) dynamic program of Bodlaender–Cygan–Kratsch–Nederlof, read as a PC. On expanders it is exponential, because every balanced split has Θ(N) boundary vertices (Lemma B3).
   - **Directed spanning trees, bounded pathwidth.** Combined with Kluk–Nederlof (ICALP 2026), this gives a small separation from nonnegativity alone: signed structured size 2^{O(k)}·poly(k) against monotone size 2^{Ω(k log log k)} (Corollary K).
4. **Unstructured signed decomposable: OPEN, with a reduction and a necessary condition.**
   - Exact signed smooth decomposable PCs for UST are, up to polynomial factors, the same as syntactically multilinear circuits for the spanning-tree polynomial ST_G (Proposition M).
   - So outcome (a) in its exact unstructured form says "ST_G has poly-size syntactically multilinear circuits." Outcome (b) says "it does not."
   - Both are open. The determinant analogue is a named open problem (Saha's SIGACT column), and the best lower bound for any explicit polynomial is Ω(n²/log² n) (Alon–Kumar–Volk).
   - **Best lead (Theorem X, cancellation) [PROVED].** A signed decomposable circuit of size s that approximates UST(G) to ε < 1 must cancel exponentially: Z⁺/λ ≥ (1 − ε)²/(s·β_*²) = e^{Ω(N)}/s, where Z⁺ is the mass of its absolute-value circuit. Any outcome-(a) circuit must live in that regime, and any outcome-(b) proof must handle it.

---

## 0. Statements in one place

Notation:
- G = (V, E) is connected, with N vertices, M edges and maximum degree D.
- P = UST(G) is viewed as a function on {0,1}^E: P(x) = 1[x is a spanning tree]/τ, with τ = τ(G).
- A split is E = A ⊔ B. It is balanced if M/3 ≤ |A| ≤ 2M/3.
- M_{A,B}(f) is the 2^{|A|} × 2^{|B|} value matrix of f.
- β(P; A|B) = ‖√P‖_op (v1 `lem:LP`), and β_* = max over balanced splits.
- A cycle is *mixed* if it has an edge in A and an edge in B.
- ν(A) is the maximum number of pairwise vertex-disjoint mixed cycles.
- ∂(A) is the set of *boundary vertices*, those incident to edges of both A and B, and k := |∂(A)|.

| # | Statement | Status |
|---|---|---|
| L1 | A smooth structured-decomposable circuit with real weights and size s, on vtree V, satisfies rank M_{Var(v), rest}(f_C) ≤ s for every vtree node v. Some node is balanced if leaf arity ≤ 2M/3. | PROVED (adaptation of de Colnet–Mengel 2021 Thm 2 / Lemma 13) |
| A | rank M_{A,B}(UST) ≥ 2^{ν(A)} for every split. | PROVED, CHECKED |
| A′ | Linked expanders: every structured signed circuit computing UST exactly, on any vtree, has size ≥ 2^{⌈r₀⌉}. | PROVED |
| C | If rank Q = r on the split and ½‖P − Q‖₁ ≤ ε, then r ≥ 2^g(1 − 2ε/(1 − δ(g))) for every g with 1 − δ(g) > 0. Here δ(g) = Pr[fewer than g good pairs] (§3). | PROVED, CHECKED |
| C4 | Under UST the number of good pairs dominates Bin(n, 2p²). | PROVED |
| Cor C | Linked expanders: an ε-approximate structured signed circuit, ε < ½, has size ≥ 2^{⌊p² n⌋}(1 − 2ε/(1 − e^{−p² n/4})) with n = ⌊⌈r₀⌉/2⌋ > 0 (otherwise vacuous). | PROVED |
| S | A normalized explicit sum of structured-decomposable squares (or shared structured-output PSD form) of charged total size S satisfies 1 − TV ≤ BC ≤ β_* √S. | PROVED, CHECKED (K = 1) |
| B | Theorem B items: | |
| | • rank M_{A,B}(UST) ≤ C(2k−2, k−1). | PROVED, CHECKED |
| | • Minor-table circuit: smooth structured-decomposable signed PC with indicator leaves, size poly(N)·2^{O(w)}. | PROVED, CHECKED |
| B3 | Every balanced split of an h-expander has k ≥ 2M/(3D(1 + D/h)). | PROVED |
| K | Arborescences on Kluk–Nederlof graphs (pathwidth k): signed structured 2^{O(k)} poly vs monotone formal arithmetic circuits (hence smooth decomposable indicator PCs) 2^{Ω(k log log k)}. | PROVED modulo [KN25, Thm 1.4] |
| M | Exact signed smooth decomposable indicator PC for UST(G) ⇔ syntactically multilinear circuit for ST_G, up to poly factors. | PROVED (standard) |
| X | A signed decomposable circuit of size s with λ = Σ f > 0 and ½‖P − f/λ‖₁ ≤ ε < 1 satisfies s ≥ (1 − ε)²(λ/Z⁺)/β_*². | PROVED, CHECKED (Z⁺ numbers) |

---

**Input model.** Unless stated otherwise, circuits have unary/indicator inputs and real weights; size counts bounded-fanin arithmetic/input gates (equivalently wires up to constants). Smooth sums have equal child scopes. Structured-decomposable products are binary and genuinely vtree-aligned: their child scopes are the two child-node scopes, not merely arbitrary disjoint vtree-node scopes. Inputs occur only at vtree leaves. One may instead use disjoint variable blocks as vtree leaves, charging their input representation cost; whenever a balanced cut is used, each block has arity ≤ 2M/3. Proposition M specifically requires indicator leaves (or explicitly supplied multilinear block representations, with conversion cost charged), never free arbitrary multivariate inputs. Theorem X allows explicitly charged multivariate inputs with its stated arity cap.

## 1. The structured rank lemma

**Lemma L1 [PROVED; adapted standard argument].** Let C be a smooth circuit with real weights over x ∈ {0,1}^E that is structured-decomposable with respect to a vtree V in the input model of §0. In particular every gate's scope equals Var(v) for some vtree node v, products follow the two vtree children, and inputs occur only at leaves. This holds after smoothing along V, by Shih–Van den Broeck–Beame–Amarilli 2019 (polynomial cost). Fix a node v, and let G_v be the set of gates with scope Var(v). Then
f_C = Σ_{g ∈ G_v} g(x_{Var v}) · h_g(x_{E∖Var v}), so rank M_{Var v, rest}(f_C) ≤ |G_v| ≤ |C|.

*Proof.* Let u range over v and its ancestors. By induction in topological order, every gate with scope Var(u) computes Σ_{g∈G_v} g·h_{g,u}, with h_{g,u} a function of x_{Var u ∖ Var v}.
- **Base.** For u = v this is the gate itself.
- **Sum gates.** Smoothness gives every child the same scope Var(u), so take linear combinations of the h's.
- **Product gates.** Children have scopes Var(u₁) and Var(u₂) for the two vtree children of u. If v lies below u₁, the second factor is a function of Var(u₂) ⊆ Var(u) ∖ Var(v), and it is absorbed into the h's.

The root has scope E. Each term g·h_g has rank one. ∎

If every vtree leaf holds at most 2M/3 variables, walk from the root to the larger child. The first node with |Var| ≤ 2M/3 has |Var| > M/3, so it is balanced. Hence a structured circuit has size ≥ min over balanced splits A of the rank of M_{A,B} for its own vtree. In particular this holds for the split at that balanced node.

This adapts the same decomposition/rank argument as de Colnet–Mengel's Theorem 2 and Lemma 13 to signed outputs at a specified vtree cut ([KR 2021 §5.2](https://proceedings.kr.org/2021/20/kr2021-0020-de-colnet-et-al.pdf)). It is also the vtree-respecting form of paper-dimers Prop 9.1.

---

## 2. Theorem A: exact rank from disjoint mixed cycles

**Theorem A [PROVED; CHECKED].** For every split A|B of any graph G, rank M_{A,B}(UST(G)) ≥ 2^{ν(A)}.

*Proof.*
1. **Setup.** Let C₁, …, C_m be pairwise vertex-disjoint mixed cycles. Pick e_i ∈ E(C_i) ∩ A and f_i ∈ E(C_i) ∩ B.
2. **A tree through all the paths.** The paths C_i − e_i are vertex-disjoint, so their union is a forest. Extend it to a spanning tree T₀ of G. No e_i is added, because e_i would close C_i. So f_i ∈ T₀ for all i.
3. **The common forest.** Put F := T₀ ∖ {f₁, …, f_m}. It has N − 1 − m edges and contains every C_i − e_i − f_i.
4. **Rows and columns.** For s, t ∈ {0,1}^m let X_s := (F ∩ A) ∪ {e_i : s_i = 1} (a row) and Y_t := (F ∩ B) ∪ {f_i : t_i = 1} (a column).
5. **Claim: X_s ∪ Y_t is a spanning tree iff t = 1 − s.**
   - **If.** Start from T₀ and, for each i with s_i = 1, swap f_i for e_i. At each step the current tree still contains C_j − e_j for every unswapped j: the cycles are vertex-disjoint, so the swapped edges are not on C_j. So the fundamental cycle of e_j is C_j, which contains f_j, and the swap gives a tree.
   - **Only if.** If s_i = t_i = 1 for some i, the set contains all of C_i. Otherwise s_i + t_i ≤ 1 for all i, and some s_i + t_i = 0, so the set has fewer than N − 1 edges.
6. **Conclusion.** The 2^m × 2^m submatrix on rows X_s and columns Y_t is (1/τ) times a permutation matrix. ∎

**Corollary A′ [PROVED].** Let G satisfy the hypotheses of v1 `thm:general`: linked with parameters ℓ, κ; edge expansion h; and r₀ as defined there. Then every smooth structured-decomposable circuit with real weights computing UST(G) exactly has size ≥ 2^{⌈r₀⌉}, whatever its vtree (leaf arity ≤ 2M/3).

*Proof.* Take the vtree's balanced node from L1. v1 `lem:mixed` with m₀ = M/3 gives ⌈t/(ℓκ)⌉ ≥ r₀ pairwise vertex-disjoint mixed cycles. Then apply Theorem A. ∎

**Remarks.**
- Vertex-disjointness is only used to make ∪(C_i − e_i) a forest, and to keep the swaps independent. Edge-disjoint cycles whose union minus one edge each is acyclic suffice.
- The permutation submatrix plays the same rank-witness role as de Colnet–Mengel's graph-theoretic induced matching gadget; it is not literally that matching. Their gadget is a vertex-variable function F_G, and ours is the cycle gadget for edge variables.

---

## 3. Theorem C: approximate rank (signed, ε < ½)

Exact rank says nothing about approximation: a small perturbation of P can have tiny rank in principle. Theorem C upgrades A to ℓ₁-approximation, using many *disjoint* permutation blocks.

**Pair gadget.**
- Group the m vertex-disjoint mixed cycles into n = ⌊m/2⌋ pairs (C_i, C_j).
- On each cycle fix e ∈ A and f ∈ B. The *arcs* of a cycle are its other edges.
- For a spanning tree t, call pair q **good** if t contains the arcs of both cycles and either of two patterns holds:
  - pattern 1: e_i, f_j ∈ t and f_i, e_j ∉ t, that is, t ∩ E(C_i) = C_i − f_i and t ∩ E(C_j) = C_j − e_j;
  - pattern 2: f_i, e_j ∈ t and e_i, f_j ∉ t, the mirror image.
- Let G(t) be the set of good pairs, g(t) = |G(t)|, and F(t) := t minus the four special edges of every good pair.
- The **block** of t is the set of trees with the same key (G(t), F(t)).

**Lemma C1 (blocks are cubes) [PROVED; CHECKED].** A block with g good pairs has exactly 2^g trees. The flip at a good pair is pattern 1 ↔ pattern 2, that is, t ↦ t − e_i + f_i − f_j + e_j. The flips at the g good pairs act freely and transitively on the block.

*Proof.*
- **The flip gives a tree.** In pattern 1, t ⊇ C_i − f_i. The fundamental cycle of f_i is C_i, which contains e_i, so t − e_i + f_i is a tree. It still contains C_j − e_j, because the cycles are vertex-disjoint. So adding e_j and removing f_j gives a tree.
- **The flip preserves the key.** No edge of another pair changes, so the pair stays good in the other pattern and every other pair keeps its status. So G and F are unchanged.
- **Freeness.** Distinct flip sets give distinct trees.
- **Transitivity.** Inside a block, a tree is determined by F and the pattern at each good pair. ∎

**Lemma C2 (block rectangles are disjoint induced permutations) [PROVED; CHECKED].**
- *Rectangle.* For block b = (G, F), rows are X = (F ∩ A) ∪ {one of e_i, e_j per q ∈ G}, and columns are Y = (F ∩ B) ∪ {one of f_i, f_j per q ∈ G}. This is a 2^g × 2^g rectangle R_b.
- *Claim (i).* P restricted to R_b is (1/τ)·(permutation).
- *Claim (ii).* The rectangles R_b of distinct blocks are disjoint.

*Proof.*
- **(i) Trees in the rectangle.**
  - If row and column choose (e_i, f_j) or (e_j, f_i) at every q ∈ G, then X ∪ Y is a tree of the block (pattern 1 or 2).
  - If at some q they choose (e_i, f_i), then X ∪ Y contains C_i, because F contains its arcs. Likewise for (e_j, f_j).
  - So exactly one column per row gives a tree.
- **(ii) Disjointness.** Recover the block from any cell (X, Y). Let G′ be the set of pairs q whose arcs lie in X ∪ Y, with exactly one of e_i, e_j in X and exactly one of f_i, f_j in Y.
  - Every q ∈ G is in G′.
  - For q ∉ G, the edges of q in X ∪ Y are those of t, for any t in the block. If the arcs of q lie in t, acyclicity of t leaves each cycle of q with e only, f only, or neither. Of the resulting patterns, only (e_i, f_j) and (f_i, e_j) give one e and one f, and those make q good, a contradiction.
  - So G′ = G. Then F = (X ∪ Y) minus the special edges of G. ∎

**Lemma C3 (ℓ₁ to rank) [PROVED].** For an n × n permutation matrix Π and any Q with rank Q ≤ r ≤ n, ‖Π − Q‖₁ ≥ n − r, where ‖·‖₁ is the entrywise sum.

*Proof.*
1. Put E := I − Πᵀ Q. Then I − E = ΠᵀQ has rank ≤ r, so K := ker(I − E) has dimension ≥ n − r.
2. With P_K the orthogonal projector onto K, E P_K = P_K. Hence n − r ≤ tr(P_K) = tr(E P_K) ≤ ‖E‖_* ‖P_K‖_op = ‖E‖_*.
3. Finally ‖E‖_* = ‖Π − Q‖_* ≤ ‖Π − Q‖₁, by the triangle inequality over unit-entry matrices. ∎

The bound is tight: take Q = Π restricted to r diagonal cells.

**Theorem C [PROVED; CHECKED].** Let Q be any real function on {0,1}^E (not necessarily nonnegative or normalised) with rank M_{A,B}(Q) ≤ r and ½‖P − Q‖₁ ≤ ε. Put δ(g) := Pr_{T∼UST}[g(T) < g]. Then for every g with α_g := 1 − δ(g) = Pr[g(T) ≥ g] > 0,
r ≥ 2^g · (1 − 2ε/(1 − δ(g))).

*Proof.*
1. Every restriction Q_b of Q to R_b has rank ≤ r. By C2, the rectangles are disjoint.
2. By C3 with scale 1/τ: 2ε ≥ Σ_b ‖P_b − Q_b‖₁ ≥ (1/τ) Σ_b (2^{g_b} − r)⁺.
3. Restrict to the blocks with g_b ≥ g. Each such block has 2^{g_b} − r ≥ 2^{g_b}(1 − r 2^{−g}).
4. Every tree lies in exactly one block. So Σ_{b: g_b ≥ g} 2^{g_b}/τ = Pr[g(T) ≥ g] = 1 − δ(g). ∎

For ε = 0 and g = n with all pairs good, this recovers rank ≥ 2^n, a weaker version of Theorem A (2^{m/2} rather than 2^m). The pairing is what makes the rectangles disjoint. With single cycles, a cell where a cycle has "neither e nor f" belongs to several blocks.

**Lemma C4 (good pairs are frequent) [PROVED].** Let p_k be the v1 `lem:local` constant of cycle C_k. Conditionally on any event determined by the edge states T ∩ E(C_l) of the other pairs, pair q = (i, j) is good with probability ≥ 2 p_i p_j. Hence g(T) stochastically dominates Bin(n, 2p²) with p = min_k p_k.

*Proof.*
1. **Pattern 1 for a fixed outside event.** Let Π be any such event of positive probability, and put F_i := T ∖ I(C_i). Vertex-disjointness gives E(C_l) ∩ I(C_i) = ∅ for l ≠ i. So both Π and {T ∩ E(C_j) = C_j − e_j} are functions of F_i (exactly as in v1 `thm:cut`). Then
   Pr[pattern 1, Π] = E[1_Π 1{T∩E(C_j) = C_j − e_j} Pr(T∩E(C_i) = C_i − f_i | F_i)] ≥ p_i Pr[T∩E(C_j) = C_j − e_j, Π].
   The same argument with F_j gives ≥ p_i p_j Pr[Π].
2. **Pattern 2.** Same argument. The two patterns are disjoint, so the pair is good with conditional probability ≥ 2 p_i p_j.
3. **Domination.** Apply step 1 with Π = the past goodness indicators. This gives Pr[good_q | good_{<q}] ≥ 2p², which yields the domination by the standard sequential coupling. ∎

**Corollary C (linked expanders) [PROVED].** Let G satisfy the hypotheses of v1 `thm:general`, with p = D^{−(ℓ−1)}, r₀ as there, and n := ⌊⌈r₀⌉/2⌋ > 0 (the n = 0 case is vacuous). Let C be any smooth structured-decomposable circuit with real weights (any vtree, leaf arity ≤ 2M/3) computing Q with ½‖UST(G) − Q‖₁ ≤ ε. Then
|C| ≥ 2^{⌊p² n⌋} · (1 − 2ε/(1 − e^{−p² n/4})).

*Proof.*
1. Lemma L1 gives a balanced node, and |C| ≥ rank there.
2. `lem:mixed` gives ⌈r₀⌉ vertex-disjoint mixed cycles of length ≤ ℓ for that split, each with p_k ≥ p (v1 `lem:local`).
3. Take g = ⌊p² n⌋ ≤ (2p²)n/2. Chernoff gives δ(g) ≤ Pr[Bin(n, 2p²) < p² n] ≤ exp(−2p² n/8).
4. Apply Theorem C. ∎

So for 0 ≤ ε ≤ ½ − η (η > 0) and N large, |C| ≥ η·2^{p² r₀/2 − 2} = e^{Ω(N)}. For the prism family of v1 (p = 1/64, r₀ ≥ N/11520) the exponent is p² r₀ ln2/2 ≥ N/(1.36·10⁸). v1's nonnegative bound is exp(p² r₀) ≥ exp(N/(4.7·10⁷)). The rates are the same up to the constant ln2/2.

**What C does not give.**
- **ε ≥ ½.** For unnormalised Q this regime is trivial (Q = 0). For Q constrained to be a nonnegative distribution, it needs a different argument; see §8 (OPEN).
- **Small graphs.** Pair counts there are 1–2 (§9), so δ(g) is large and the numerics can only confirm the lemmas, not show the rate.

---

## 4. Theorem S: normalized squares/PSD forms of structured-decomposable amplitudes

**Perron sign lemma [PROVED; CHECKED].** For W ≥ 0 entrywise and any matrix S with |S_ab| ≤ 1, ‖W ∘ S‖_op ≤ ‖W‖_op.

*Proof.* |xᵀ(W∘S)y| ≤ Σ W_ab |x_a||y_b| ≤ ‖W‖_op‖x‖‖y‖. ∎

The vector version also holds. Let W̃ := [W∘S₁ | … | W∘S_K] with Σ_m S_m(a,b)² ≤ 1. Then ‖W̃‖_op ≤ ‖W‖_op: the same computation, with Cauchy–Schwarz over m.

**Theorem S [PROVED].**
- **Setting.** Let d = (d₁, …, d_K) be full-scope outputs of one smooth structured-decomposable circuit in the §0 input model, with real weights of size s on vtree V and leaf arity ≤ 2M/3. Let Q = ‖d(x)‖²/Z, where Z := Σ_x ‖d(x)‖² > 0.
- **Coverage.** This covers the normalized square of a smooth structured-decomposable circuit (K = 1); an explicit sum of such squares on separate vtrees, charged S = Σ_k s_k (not the DAG union once); and d(x)ᵀKd(x)/Z with K ⪰ 0 and d from one shared smooth structured circuit. For the latter, write K = LLᵀ and replace d by Lᵀd: mixing preserves the same cut space, so no fictitious expanded circuit is charged. It also covers Loconte–Mengel–Vergari's compatible Σ²_cmp class and their top-level psd class (Proposition 2). It does not automatically cover nested PSD units, compact μSOCS with exponentially many implicit squares, or positive-unital/quantum architectures without a separate interface-to-size argument.
- **Claim.** For the balanced split A|B of V, BC(P, Q) ≤ β(P; A|B)·√s, so 1 − TV(P, Q) ≤ β_*√s.
- **Several vtrees.** If instead Q = Σ_k c_k²/Z with each c_k structured on its own vtree and of size s_k, then BC ≤ β_* √(Σ_k s_k).

*Proof.*
1. **Stacked rank.** Stack the value matrices D_m of the d_m side by side: D̃ = [D₁ | … | D_K]. By L1, each D_m = U H_mᵀ, where U lists the gates with scope Var(v). So rank D̃ ≤ s.
2. **Signs.** Let S_m(x) := d_m(x)/‖d(x)‖, with S_m(x) := 0 when ‖d(x)‖ = 0.
3. **Bound.** Σ_x √P(x)‖d(x)‖ = ⟨W̃, D̃⟩ with W = √P. This is ≤ β‖D̃‖_* ≤ β√s‖D̃‖_F = β√s√Z. Divide by √Z.
4. **Several vtrees.** Use √(Σc_k²) ≤ Σ|c_k|, apply the K = 1 bound to each c_k at its own balanced split, and finish with Cauchy–Schwarz. ∎

**Consequence.** On linked expanders the covered normalized squares/PSD forms of structured-decomposable amplitudes have charged size ≥ (1 − TV)²(1 − p²)^{−r₀}. This is v1's nonnegative lower-bound expression, for every TV < 1, not a claim of identical optimal sizes.

**Inner decomposability is essential.** On Boolean edge indicators let H(x) := [z^{N−1}] Π_e (1 − x_e + x_e z) = [|x| = N−1]. Both H and det L_r(x) have polynomial-size ordinary arithmetic circuits. Then d(x) := det L_r(x)·H(x)/√τ has d(x)² = UST(x) exactly: the determinant counts trees contained in x, and on N−1 edges it is one precisely for a tree. This polynomial-size amplitude is not decomposable. Unrestricted squaring therefore can represent UST even on expanders.

Thus the expressivity gains that Loconte–Mengel–Vergari (AAAI 2025) prove for the covered SOCS classes over monotone PCs do not evade this UST lower bound. Their gains come from functions of small square-root rank but large nonnegative rank. UST has *both* exponential at every balanced split, because β is exponentially small.

---

## 5. Theorem B: the positive route (determinant as a structured signed PC)

**Lemma B1 (rank upper bound) [PROVED; CHECKED].** rank M_{A,B}(UST) ≤ C(2k − 2, k − 1), where k = |∂(A)| ≥ 1. For a trivial split (k = 0), rank is 1 instead.

*Proof.*
1. **Matrix tree theorem.** τ·P has the same value matrix as the coefficient matrix of ST_G(x) = det L(x)[V − r] (coefficients 0/1), with root r ∈ ∂. Split L = L_A + L_B, where L_A is supported on V_A × V_A, the vertices touched by A, and similarly L_B.
2. **Laplace expansion.** det(X + Y) = Σ_{|R|=|C|} ± det X[R, C] det Y[R̄, C̄].
3. **Which terms survive.** A term is nonzero only if R, C ⊆ V_A and R̄, C̄ ⊆ V_B. So R = U_A ∪ R′ and C = U_A ∪ C′ with R′, C′ ⊆ ∂ − r and |R′| = |C′|. Here U_A = V_A ∖ ∂.
4. **Count.** That is Σ_j C(k−1, j)² = C(2k−2, k−1) terms, each a product (polynomial in x_A) × (polynomial in x_B). The coefficient matrix is then a sum of that many rank-one matrices. ∎

**Theorem B (minor-table circuit) [PROVED; CHECKED].**
- **Setting.** Fix an edge order and the caterpillar or balanced vtree over it, and let w be the maximum boundary size of a vtree node.
- **Tables.** Delete the Laplacian root r throughout. Put U_v := V_v ∖ (∂_v ∪ {r}). Every vtree node v holds the table D_v[R, C] := det L_v[U_v ∪ R, U_v ∪ C] over R, C ⊆ ∂_v ∖ {r} with |R| = |C|. Signed off-diagonal minors are allowed.
- **Combination.** A node's table is computed from its children's tables by the generalized Laplace expansion. Every summand is a product of one entry of each child.
- **Leaves.** The leaf table of edge e holds the order-0 minor 1 and the order-1 minors ±x_e (or 0).

Replace the order-0 leaf entry by the indicator [x_e = 0] and x_e by [x_e = 1]. The resulting circuit is a smooth structured-decomposable PC with signed weights (±1) and indicator leaves. It computes τ·UST(G) exactly and has size M·2^{O(w)}, crudely ≤ M·4^{5w/2}.

*Proof of the Boolean form.*
1. Let y⁰_e, y¹_e be the two leaf symbols. Every leaf entry is linear in exactly one of them, and every gate multiplies one entry from each child. So the computed polynomial Q(y⁰, y¹) is homogeneous of degree 1 in each pair {y⁰_e, y¹_e}.
2. Q(1, x) = ST_G(x) is the Laplace identity.
3. A polynomial of that homogeneity is determined by its dehomogenisation, so Q = Σ_S 1[S is a tree] Π_{e∈S} y¹_e Π_{e∉S} y⁰_e. Setting y¹ = [x = 1] and y⁰ = [x = 0] gives 1[x is a tree].
4. **Smoothness and structure.** Each product combines the two child-node scopes, and all summands have scope Var(v). When a scope-carrying leaf function is 1 or x_e (rather than the homogenizing order-0 entry above), implement it as δ0+δ1 or δ1 respectively; zero is 0·(δ0+δ1). Thus such constants do not silently drop edge scopes.
5. **Size.** Size counts the resulting bounded-fanin arithmetic/indicator circuit, including summand products, binary sum trees, and leaf implementations, not just distinct table states. Table sizes are ≤ 4^w. A combination step sums, per entry, over row and column subsets of ∂_a ∪ ∂_b, so it costs ≤ 4^{|∂_a ∪ ∂_b|} per entry. Every vertex of ∂_a ∪ ∂_b ∪ ∂_v lies in at least two of the three sets: a vertex of ∂_a has an edge outside a, which is in b or outside v. Hence |∂_a ∪ ∂_b| ≤ 3w/2. ∎

**Corollary (bounded width).** On graphs of bounded branchwidth or pathwidth, UST has a poly-size signed structured PC. Nonnegative structured PCs of size N·Bell(w)^{O(1)} = N·w^{O(w)} also exist there (the standard monotone connectivity DP over boundary partitions; not re-derived here), so on bounded width there is no exponential separation. Corollary K below finds a superpolynomial one.

**Lemma B3 (expanders have large boundary) [PROVED].** If G has edge expansion h and A|B is balanced, then k ≥ 2M/(3D(1 + D/h)).

*Proof.*
1. Let S be the set of vertices with all their edges in A, T the set with all their edges in B, and ∂ the rest.
2. No edge joins S to T, so every edge leaving S or T ends in ∂. Hence |∂_E S| ≤ D k and |∂_E T| ≤ D k.
3. Since |A| ≥ M/3 and every A-edge has both endpoints in S ∪ ∂: |S| + k ≥ 2|A|/D ≥ 2M/(3D). Likewise |T| + k ≥ 2M/(3D).
4. One of S, T has at most N/2 vertices, say S. Then h|S| ≤ Dk, so (1 + D/h)k ≥ |S| + k ≥ 2M/(3D). ∎

So on expanders the B1 upper bound and the A lower bound are both 2^{Θ(N)}, and the determinant circuit is exponential for every vtree. This is not a limitation of the construction alone: Corollary A′ says *every* structured signed circuit is exponential.

**Corollary K (directed, bounded pathwidth) [PROVED modulo the cited theorem].**
- **Directed circuit.** Theorem B works verbatim for the out-arborescence polynomial with a fixed root, via the directed matrix-tree theorem, L_in = Σ_{u→v} x_{uv}(e_v − e_u)e_vᵀ. Summing over the N roots gives a smooth structured-decomposable signed PC for the uniform law on spanning arborescences (any root), of size poly(N)·2^{O(k)} on pathwidth-k digraphs.
- **Monotone lower bound.** Kluk–Nederlof [arXiv 2512.23121, ICALP 2026, Thm 1.4] give digraphs of pathwidth ≤ k on k^{O(1)} vertices on which every tropical circuit for DST_G has 2^{Ω(k log log k)} gates. Their polynomial sums out-trees over all roots, as in the rooted-determinant sum above, and their semiring is (max,+) over real valuations; negating weights exchanges maximum and minimum weight.
- **From formal circuits to tropical circuits.** A monotone formal arithmetic circuit for the arborescence polynomial tropicalises (Jerrum–Snir) to DST_G at the same size. This includes a smooth decomposable nonnegative indicator PC for the law: [x = 0] → 1 and [x = 1] → x_e gives that formal polynomial, up to normalization, by the network-polynomial identity in M. It does not include arbitrary nondecomposable Boolean PCs: x·x = x on Booleans becomes x² formally, and [x_e = 0][x_e = 1] vanishes on Booleans but becomes x_e after substitution.
- **Conclusion.** Signed structured 2^{O(k)} versus monotone formal arithmetic circuits for the arborescence polynomial (hence smooth decomposable nonnegative indicator PCs) of size 2^{Ω(k log log k)}. This is a genuine but small gap from nonnegativity alone within the decomposable PC model, of order S^{Ω(log log log S)}.

Checks: `bw_directed.py` matches the weighted determinant and brute-force arborescence counts (§9). The Boolean leaf version is argued, not run, for the directed case.

Novelty: KN state the gap against "algorithms that rely on algebra" themselves. The only addition here is that the algebraic side is a structured-decomposable PC.

---

## 6. Unstructured signed circuits: equivalence and the cancellation theorem

**Proposition M [PROVED; standard, cf. Martens–Medabalimi 2014 §2].** Up to polynomial factors, the following are equivalent:
- (i) a smooth decomposable PC with real weights and indicator leaves of size s computing UST(G) exactly;
- (ii) a syntactically multilinear arithmetic circuit of size s computing ST_G.

*Proof.*
- **(i) ⇒ (ii).** The formal polynomial in the leaf symbols is homogeneous of degree 1 in every pair, as in Theorem B. On Boolean points it equals UST, so it is Σ_S UST(S) Π y¹ Π y⁰. Set y⁰ = 1 and y¹ = x, and multiply by τ.
- **(ii) ⇒ (i).** Smooth the circuit (polynomial cost; Shih et al. 2019). At a sum gate, multiply each child by Π y⁰_e over the variables of the gate's scope missing from that child's scope. The result has the pair-homogeneity above, and its dehomogenisation is still ST_G. So it computes τ·Σ_S 1[S tree] Π y¹ Π y⁰ by the uniqueness argument. ∎

**Consequence.**
- Outcome (a) in its exact unstructured form is the statement that ST_G has polynomial-size syntactically multilinear circuits for an expander family. That is open.
- The determinant analogue is an explicit open question: [Saha, SIGACT News guest column](https://www.csa.iisc.ac.in/~chandan/research/sigact.pdf), and the discussion in [Alon–Kumar–Volk, CCC 2018](https://drops.dagstuhl.de/storage/00lipics/lipics-vol102-ccc2018/LIPIcs.CCC.2018.11/LIPIcs.CCC.2018.11.pdf).
- Outcome (b) unstructured exact would be a superpolynomial syntactically multilinear lower bound for an explicit polynomial in VP. The best known lower bound for *any* explicit polynomial is Ω(n²/log² n).
- Approximate (b) implies exact (b), so it is at least as hard.
- Kush 2026 ([ECCC TR26-043](https://eccc.weizmann.ac.il/report/2026/043/)) shows that full min-partition rank is achieved by poly-size multilinear ABPs. So no rank measure alone can do it.

**A route that fails, recorded so nobody retries it.**
- **The idea.** The signed version of v1 `lem:decomp` holds: f_C = Σ_{u≤s} f_u h_u over balanced splits A_u|B_u, with signed factors. Raz's random-partition rank argument applied to this two-factor form gives, against a random equipartition (Y, Z), rank(f_u h_u) ≤ 2^{M/2 − |d_u|}, where d_u = |A_u ∩ Y| − |A_u ∩ Z|.
- **Why it fails.** d_u is typically about √M, and Pr[|d_u| = j] = O(1/√M). So E Σ_u rank ≤ s·2^{M/2}·O(1/√M), which gives only s = Ω(√M)·rank(f)/2^{M/2}. That is polynomial at best. Raz needs many factors, as in formulas or log-depth, not two.

**Theorem X (cancellation is necessary) [PROVED].**
- **Setting.** Let C be a smooth decomposable circuit with real weights and real-valued leaves of arity ≤ 2M/3, computing f with λ := Σ_x f(x) > 0. Assume 0 ≤ ε := ½‖P − f/λ‖₁ < 1 (needed before squaring 1 − ε).
- **Absolute circuit.** Let C⁺ be C with every weight and leaf function replaced by its absolute value, and Z⁺ := Σ_x f_{C⁺}(x). This is computable in one pass by marginalising C⁺; with indicator leaves, set both indicators to 1.
- **Claim.** |C| ≥ (1 − ε)² · (λ/Z⁺) / β_*².

*Proof.*
1. **Pointwise domination.** By induction, f_{C⁺} ≥ |f| pointwise, since |Σ w g| ≤ Σ|w||g| and |gh| = |g||h|. Also Z⁺ ≥ Σ|f| ≥ λ.
2. **Overlap.** Q⁺ := f_{C⁺}/Z⁺ is a nonnegative smooth decomposable PC of the same size. Using √(ab) ≥ min(a, b) and min(a, κb) ≥ κ min(a, b) for κ ≤ 1:
   BC(P, Q⁺) ≥ Σ_x √(P |f|/Z⁺) ≥ √(λ/Z⁺) Σ_x min(P, |f|/λ) ≥ √(λ/Z⁺) Σ_x min(P, f/λ) = √(λ/Z⁺)(1 − ε).
3. **v1 reduction.** v1 `lem:decomp` and the proof of `lem:overlap` bound BC(P, Q⁺) itself, not just 1 − TV, by √|C|·β_*. ∎

**Corollary X.** On linked expanders, any signed decomposable circuit for UST, exact or ε-approximate with ε < 1, has cancellation ratio Z⁺/λ ≥ (1 − ε)² (1 − p²)^{−r₀} / |C|. A polynomial-size outcome-(a) circuit must therefore have *exponentially* more absolute mass than value.

This is the precise sense in which nonnegativity matters. Signed circuits can only beat v1 by cancelling e^{Ω(N)} of their mass, and any lower bound for (b) must handle circuits that do so.

Theorem B's circuit already cancels exponentially even where it is not needed. On ladders, log₂(Z⁺/τ) grows about 0.6 bits per rung (§9), although nonnegative poly-size PCs exist there.

Related prior work: the low-discrepancy result of [Chattopadhyay–Datta–Ghosal–Mukhopadhyay, ITCS 2022](https://arxiv.org/abs/2109.06941) for oriented spanning trees on K_n under every near-balanced partition is the natural tool for a *weighted* version of X, where Σ_u ‖f_u‖·‖h_u‖ is bounded. The ℓ² form is 1 ≤ √τ·β_*·Σ_u ‖f_u ⊗ h_u‖₂. It follows from 1 = Σ_u ⟨f_u ⊗ h_u, 1_supp⟩ and ‖1_supp‖_op = √τ·β by Cauchy–Schwarz. It bounds the cancellation as well, not the size.

---

## 7. Answers to the five plan items

1. **Structured first.**
   - Exact: rank ≥ 2^{ν}, which is exponential (Theorem A and Corollary A′).
   - Approximate, ε < ½: exponential (Corollary C).
2. **Scope of the rank bound.**
   - On linked expanders: exponential for *every* balanced split, hence every vtree.
   - On small graphs, adversarial *edge* splits (not vertex-induced) reach low rank by making k small: for example rank 13 on Petersen with k = 4. Lemma B3 shows this cannot persist on expanders.
3. **Exact small-graph computations.**
   - Real rank (float and mod p): §9.
   - ε-rank: achievable by truncated SVD, with the proved lower bound from Theorem C.
   - Square-root rank / PSD route: the √r·β bound and the squared-SVD achievable curves.
   - PSD rank itself: not computed. It is NP-hard in general. Theorem S suffices for PSD quadratic forms over shared structured-decomposable features; arbitrary matrix PSD factorizations and general PSD circuits are different size models and are not covered.
4. **Unstructured signed.**
   - Signed rectangle covers: the MM/v1 decomposition survives with signs, but `lem:overlap` does not.
   - Exact equivalence to syntactically multilinear circuits for ST_G (Proposition M).
   - Raz's two-factor partition argument is too weak (§6).
   - Rank measures are barred by TR26-043.
   - Discrepancy controls only bounded-cancellation circuits (Theorem X).
5. **Positive route.**
   - Theorem B is poly-size on bounded width, and is the BCKN determinant approach.
   - It cannot extend to expanders in structured form (A′).
   - It extends in unstructured form iff ST_G ∈ syntactically-multilinear-VP, which is OPEN.

---

## 8. Open problems and where things stop

- **(O1) Unstructured signed decomposable circuits for UST on expanders [OPEN].** Equivalent to multilinear circuit complexity of ST_G (Proposition M).
  - A *superpolynomial* lower bound is beyond current techniques.
  - An *upper* bound would need an entirely new multilinear algorithm for the determinant/ST, and Corollary X says it would have to cancel e^{Ω(N)}.
  - **Most promising concrete lead:** prove the lower bound for circuits with Z⁺/λ ≤ e^{o(N)}·poly, which is immediate by X, then push into intermediate regimes, e.g. Z⁺/λ ≤ e^{cN} for small c. This could use discrepancy in the CDGM style.
- **(O2) ε ∈ [½, 1) for structured signed circuits whose output is a distribution [OPEN].** Theorem C is sharp at ½ only because Q = 0 is allowed.
  - **Conjecture [CONJ].** With Q ≥ 0 and ΣQ = 1, the structured size is still ≥ (1 − TV)^{O(1)} e^{Ω(N)}.
  - **Obstacle.** Inside one permutation block, a nonnegative rank-r matrix can cover all diagonal cells (Q = J/n²) at the price of mass n². The global mass budget must be traded against coverage across blocks, which C3 does not do.
- **(O3) Constants.** All rates inherit v1's lossy p = D^{−(ℓ−1)}. Nothing here improves them.
- **(O4) Positive unital / quantum-information circuits** (arXiv 2506.01824). Not checked against Theorem S. If their output is tr ρ(x) with ρ(x) = Σ_m d_m(x)d_m(x)† for a structured d, then S applies verbatim [CONJ].

---

## 9. Numerics (all CHECKED; scripts in `toy/r3/signed/`)

### Run 2: `rank_blocks.py` via `run_all.sh`

Seed 1, lmax 6, one Fiedler split plus random balanced edge splits. Raw output: `logs/<graph>.jsonl`.

Every line passed every loud check:
- float rank = mod-p rank (when ≤ 1500 rows);
- 2^{ν} ≤ rank ≤ C(2k−2, k−1);
- Lemma C1/C2 block checks (11.9 million cells over all lines);
- the Theorem C bound ≤ exact rank;
- the Perron sign lemma on 20 random sign patterns;
- the squared-circuit bound BC ≤ √r β on 20 random signed rank-r functions.

Pappus: 3 of 6 splits exceeded `--maxdense 6000`, so no rank was computed for them. The pair-block checks still ran.

Column legend:
- **ranked/all**: splits with a computed rank / all splits run.
- **2^ν lb**: the lower bound 2^ν from a greedy packing of disjoint mixed cycles (lower bound on 2^{ν(A)}).
- **max ‖√P∘S‖/β**: the largest ratio over the random sign patterns.
- **max BC/(√r β)**: the largest ratio over the random signed rank-r functions.
- **pair good / 2p² (1 pair)**: observed probability that the single pair is good, against its lower bound.

| graph | N | M | τ | ranked/all | k range | rank range | 2^ν lb | 1/β² range | max ‖√P∘S‖/β | max BC/(√r β) | pair good / 2p² (1 pair) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K33 | 6 | 9 | 81 | 6/6 | 4–6 | 8–16 | 2 | 2.2–3.3 | 1.000 | 0.695 | – |
| cube | 8 | 12 | 384 | 6/6 | 4–7 | 10–23 | 4 | 2.6–4.5 | 1.000 | 0.592 | 0.0208 / 0.0027 |
| petersen | 10 | 15 | 2000 | 6/6 | 5–8 | 20–75 | 4 | 3.3–4.5 | 0.918 | 0.568 | 0.0050 / 0.0003 |
| prism5 | 10 | 15 | 1805 | 6/6 | 4–10 | 13–91 | 4 | 2.6–5.1 | 0.994 | 0.605 | 0.0266 / 0.0027 |
| prism6 | 12 | 18 | 8100 | 6/6 | 4–11 | 14–120 | 4–8 | 2.4–7.6 | 1.000 | 0.682 | 0.0316 / 0.0027 |
| prism7 | 14 | 21 | 35287 | 6/6 | 4–13 | 14–256 | 4–8 | 2.9–11.0 | 0.995 | 0.569 | 0.0340 / 0.0027 |
| prism8 | 16 | 24 | 150528 | 6/6 | 4–13 | 14–625 | 4–16 | 2.6–13.4 | 1.000 | 0.667 | 0.0359 / 0.0027 |
| mobladder7 | 14 | 21 | 35301 | 6/6 | 4–12 | 14–249 | 4–8 | 2.9–9.5 | 0.996 | 0.500 | 0.0340 / 0.0027 |
| heawood | 14 | 21 | 50421 | 6/6 | 6–12 | 69–380 | 4 | 4.4–8.2 | 0.903 | 0.511 | 0.0016 / 0.000034 |
| rr3_14 | 14 | 21 | 42699 | 6/6 | 5–12 | 34–432 | 4–8 | 3.4–9.1 | 0.940 | 0.439 | 0.0077 / 0.0003 |
| rr3_16 | 16 | 24 | 144254 | 6/6 | 3–14 | 5–846 | 2–16 | 2.1–12.3 | 1.000 | 0.632 | 0.0329 / 0.0003 |
| mobius_kantor | 16 | 24 | 248832 | 6/6 | 6–14 | 107–1208 | 4 | 5.1–11.4 | 0.975 | 0.453 | 0.0019 / 0.000034 |
| pappus | 18 | 27 | 1259712 | 3/6 | 7–17 | 160–2784 | 8 | 5.8–16.1 | 0.731 | 0.386 | 0.0021 / 0.000034 |

Observations:
- The sign lemma is attained (ratio 1.000) only on splits where the top singular vector's sign pattern is hit.
- Observed good-pair probabilities exceed the 2p² bound (with p = the smallest `lem:local` constant among the packed cycles) by 8–110×. The v1 constant is lossy.
- Exact ranks are far below 2^{min(|A|,|B|)} but well above 2^{ν}. The truth lies between the A and B1 bounds.

**ε-curves at Fiedler splits** (`curves` field). Each cell is sq / sg / lb, where:
- sq = TV of the rank-r squared-SVD approximation;
- sg = half-ℓ₁ error of the rank-r truncated SVD of P (signed);
- lb = the Theorem S lower bound 1 − √r β on TV for any square-root-rank-r squared circuit.

| graph | k / rank | r = 1 | r = 4 | r = 16 | r = 64 |
|---|---|---|---|---|---|
| petersen | 5 / 20 | .65 / .42 / .45 | .32 / .28 / 0 | .06 / .05 / 0 | – |
| heawood | 6 / 69 | .72 / .44 / .52 | .44 / .40 / .05 | .20 / .19 / 0 | 0 / 0 / 0 |
| prism8 | 4 / 14 | .61 / .31 / .38 | .24 / .23 / 0 | – | – |

Heawood at a random split (k = 9, rank 254), at r = 1 / 8 / 32 / 128:
- sq: .77 / .45 / .24 / .04;
- sg: .51 / .45 / .26 / .06.

So beyond r ≈ 4, signed truncation is no better than squared. Both are governed by the same spectrum, as Theorem S predicts. The lb column becomes trivial once r ≳ 1/β², which on these small graphs is 3–16.

**Theorem C at small sizes.** With one pair, the bound is 2 at ε = 0 and trivial at ε ≥ 0.1 (δ(1) ≈ 0.96). This is expected, since the rate needs n·p² ≫ 1.

### Adversarial edge splits (`adv_splits.py`, `run_adv.sh`)

The search minimises rank over balanced edge splits. Raw output: `logs/adv-*.jsonl`. The K33 and cube rows come from interactive runs of the same script.

| graph | mode | min rank | k | 2^ν lb | C(2k−2,k−1) |
|---|---|---|---|---|---|
| K33 | exhaustive | 5 | 3 | 2 | 6 |
| cube | exhaustive | 10 | 4 | 4 | 20 |
| petersen | exhaustive (14443 splits) | 13 | 4 | 4 | 20 |
| prism5 | exhaustive | 13 | 4 | 4 | 20 |
| prism6 | descent | 13 | 4 | – | 20 |
| heawood | descent | 52 | 8 | 4 | 3432 |
| rr3_14 | descent | 10 | 5 | 4 | 70 |

The minimisers concentrate A around a few vertices, for example all edges at two adjacent vertices of Petersen, to make k small. This is a finite-size effect (Lemma B3).

### Theorem B circuit (`bw_circuit.py`; re-run log `logs/bw_circuit.jsonl`)

Every run matches det L[V − r] to relative error ~1e−15.

Boolean indicator check: 300 inputs per graph, 0 mismatches on all 8 graphs of the re-run (K33, cube, petersen, ladder5, ladder10, ladder20, rr3_12, prism6), and earlier also on rr3_14.
- Between 157 and 247 of the inputs were trees on each graph except ladder20.
- On ladder20 (M > 30) trees are not enumerated, so all 300 inputs were random non-trees there.

| graph | max boundary w | Laplace terms | table entries |
|---|---|---|---|
| K33 | 3 | 139 | 135 |
| cube | 4 | 411 | 310 |
| petersen | 5 | 2779 | 1091 |
| prism6 | 7 | 8315 | 4364 |
| rr3_12 / 16 / 20 | 6 / 7 / 7 | 4450 / 28101 / 17388 | 2006 / 8872 / 8548 |
| ladder5 / 10 / 20 / 40 | 4 | 496 / 1968 / 5704 / 13330 | 410 / 1139 / 2761 / 6001 |
| prism8 / prism12 | 8 | 208357 / 401697 | 31280 / 50906 |

Cancellation ratio of the same circuit (Theorem X's Z⁺):

| graph | τ | Z⁺ | log₂(Z⁺/τ) |
|---|---|---|---|
| K33 | 81 | 117 | 0.53 |
| cube | 384 | 684 | 0.83 |
| petersen | 2000 | 3724 | 0.90 |
| rr3_12 | 8256 | 23320 | 1.50 |
| prism6 | 8100 | 23172 | 1.52 |
| ladder5 | 209 | 845 | 2.02 |
| ladder10 | 151316 | 6.55e6 | 5.43 |
| ladder20 | 7.93e10 | 3.92e14 | 12.27 |

### Directed variant (`bw_directed.py`, `run_bw.sh`; log `logs/bw_directed.jsonl`)

Out-arborescences at root 0 of the bidirected graph, checked against the weighted determinant (rel. err. ≤ 7e−15). At x = 1 they match brute-force counts: K33 81, cube 384, petersen 2000, ladder5 209.

| graph | arcs | w | size |
|---|---|---|---|
| K33 | 18 | 4 | 1518 |
| cube | 24 | 4 | 1453 |
| petersen | 30 | 5 | 8119 |
| ladder5 | 26 | 4 | 1802 |
| ladder10 | 56 | 5 | 10220 |
| ladder20 | 116 | 6 | 36880 |

---

## 10. Literature and novelty

Searches were run Nov 4 2026 (deep web search, several rounds), and the key pages were read.

**Squared / SOCS / PSD circuits.**
- [Loconte–Mengel–Vergari, "Sum of Squares Circuits", AAAI 2025](https://arxiv.org/abs/2408.11778).
  - Squared circuits can be exponentially weaker than monotone structured PCs (via the square-root-rank lower bound on SUM).
  - SOCS circuits are exponentially more expressive than both.
  - Structured decomposability is assumed. App. B.1 generalises MM Thm 38 to signed circuits.
  - Nothing on UST, and no ε-approximation bound of our kind.
- [Wang–Van den Broeck, AAAI 2025](https://starai.cs.ucla.edu/papers/WangAAAI25.pdf): monotone and squared structured PCs are incomparable (exact, rank-based).
- [arXiv 2506.01824](https://arxiv.org/abs/2506.01824): positive unital circuits; not checked against S (O4).
- [arXiv 2512.17090](https://arxiv.org/abs/2512.17090) (ICLR 2026) and [arXiv 2412.07883](https://arxiv.org/abs/2412.07883): squaring techniques, not lower bounds.

**Rank technique for structured circuits.**
- [de Colnet–Mengel, KR 2021](https://arxiv.org/abs/2110.13014), §5.2: positive (subtraction-allowed) structured ACs need exponential size, via value-matrix rank and induced matchings in expanders. Our L1 adapts the decomposition/rank argument of their Theorem 2 and Lemma 13 to signed outputs at a specified vtree cut.
- Theorem A is the same technique with a new gadget (disjoint mixed cycles) for a different function (UST). Novelty is moderate at best.
- The ITCS 2026 ROABP closure paper uses the same induced-matching rank method.

**Approximation.**
- [arXiv 2506.01281](https://arxiv.org/abs/2506.01281) (NeurIPS 2025): hardness of approximating distributions with PCs, and a decomposable vs deterministic gap under bounded f-divergence.
- [arXiv 2605.12940](https://arxiv.org/abs/2605.12940): separation rank of structured PCs.
- Neither gives an ℓ₁-approximate rank bound against signed structured circuits for an explicit natural distribution.
- I found no prior statement of Theorem C (disjoint induced-permutation blocks + the nuclear-norm threshold + the good-pair domination) or of Theorem S for UST. The ingredients are elementary, so treat the novelty as **"not found in a non-exhaustive search"**.

**Determinant approach.**
- [Bodlaender–Cygan–Kratsch–Nederlof, ICALP 2013 / Inf. Comput. 2015](https://arxiv.org/abs/1211.1505): counting via the matrix-tree theorem over tree decompositions with generalized Laplace expansion. Theorem B is this, written as a structured PC with indicator leaves.
- The rank bound C(2k−2, k−1) is implicit there. **Not new.**

**Pathwidth lower bounds.**
- [Kluk–Nederlof, arXiv 2512.23121 / ICALP 2026](https://arxiv.org/abs/2512.23121): tropical lower bound 2^{Ω(k log log k)} for directed spanning tree at pathwidth k. Corollary K is a two-line combination; KN already frame it as "pure DP vs algebra". **Minor.**

**Monotone spanning-tree bounds.**
- [Chattopadhyay–Datta–Ghosal–Mukhopadhyay, ITCS 2022](https://arxiv.org/abs/2109.06941): 2^{Ω(n)} monotone bound for ST on expanders; ε-sensitive bounds; low discrepancy for oriented spanning trees under every near-balanced partition.
- [ECCC TR20-191](https://eccc.weizmann.ac.il/report/2020/191/): negations give strongly exponential savings, via the same polynomials.
- Hrubeš, [ε-sensitive monotone computations](https://users.math.cas.cz/~hrubes/PDFs/Epsilon.pdf): strong enough ε-sensitive monotone bounds imply general circuit lower bounds. Our Theorem X is a PC-flavoured cancellation measure in the same spirit, but different: it uses the ratio of absolute mass, not an additive perturbation.

**Multilinear circuits.**
- Open whether det has poly-size syntactically multilinear circuits: [Saha, SIGACT column](https://www.csa.iisc.ac.in/~chandan/research/sigact.pdf).
- Best explicit lower bound Ω(n²/log² n): [Alon–Kumar–Volk](https://arxiv.org/abs/1708.02037).
- Rank barrier: [Kush, ECCC TR26-043](https://eccc.weizmann.ac.il/report/2026/043/).
- Proposition M is standard ([Martens–Medabalimi 2014](https://arxiv.org/abs/1411.7717)).

---

## 11. Weakest steps (for a reviewer)

1. **Lemma C2 (ii), disjointness of block rectangles.** The whole ε-rank bound rests on recovering the good-pair set G from a cell. The case analysis uses that a non-good pair, with arcs present, cannot show "exactly one e in X and exactly one f in Y". Check the case where the arcs of q are in t but t contains e_i and f_i of the *same* cycle: that contains C_i, so it cannot happen.
   - Numerically no overlap was found in 11.9M cells with 1 pair per split.
   - Multi-pair runs found none either: `logs/multipair.jsonl`, with lmax 4 and squares as the mixed cycles.
     - ladder9: 2 pairs, ≈ 45k cells per split, 5 random splits, up to g = 2 good pairs in a block.
     - ladder12: 3 pairs, 2.46M cells, one random split, up to g = 3 (256 trees with all three pairs good).
   - All C1/C2 checks passed. The ladder12 job was stopped by its 110 s timeout after two splits; this was not a check failure.
2. **Lemma L1 hypotheses.** It assumes every gate's scope is a vtree-node scope (strong structure after smoothing). Circuits that are structured only in the weak sense (product children ⊆ Var(v_l), Var(v_r)) need the smoothing step. Its cost (Shih et al. 2019) is polynomial, but the exact size accounting was not redone here.
3. **Theorem S's several-vtree version** uses β_* = max over balanced splits. This is fine for UST on linked expanders by v1, but it is weaker than the single-vtree statement if β varies a lot across splits.
4. **Corollary K** depends on the exact form of KN Thm 1.4 (graphs, root convention). I read the statement in the arXiv HTML, not the full proof.

---

## 12. Files

New, all under `toy/r3/signed/`:
- `r3common.py`: graphs, tree enumeration, compressed split matrices, ranks, mixed cycles.
- `rank_blocks.py`: rank, β, ε-curves, Theorem A/C checks, Perron lemma, squared bound.
- `run_all.sh`: driver for `rank_blocks.py`.
- `adv_splits.py` and `run_adv.sh`: minimum rank over balanced edge splits.
- `bw_circuit.py` and `run_bw_circuit.sh`: Theorem B circuit with Boolean and absolute modes.
- `bw_directed.py` and `run_bw.sh`: directed variant.
- logs: `logs/*.jsonl`, `logs/*.err`, `logs/all.done`, `logs/adv.done`, `logs/bw.done`, `logs/bw_circuit.done`, `logs_run_*.txt`.
- `logs/multipair.jsonl` comes from direct `rank_blocks.py --lmax 4` runs on ladder9 and ladder12.

`toy/r2/common.py` was imported and not modified. No other notes were edited. No commits were made.

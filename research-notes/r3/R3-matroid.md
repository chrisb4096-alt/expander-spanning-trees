# R3: matroid generalization (uniform random basis)

Written Nov 4 2026 (CST). Scope: the matroid direction of R3. The paper source `notes/paper-ust/` was read, not modified.

**Verdict: PROVED for the local (short-certificate) version, for every matroid. PROVED for the global (leverage) version for balanced matroids, which include regular matroids and binary matroids with no S8 minor. OPEN for non-balanced matroids in the global version (Lemma 6 is false for S8).**

Reviewed (Sol 6.1, R3-matroid-review.md): holds with fixes, 8/10; fixes applied.

Labels: **[PROVED]** means a full proof is in this note. **[CHECKED]** means exact small-case numerics (§11). **[CONJ]** means conjectured.

---

## 0. Statements in one place

Notation: M is a matroid on E, |E| = n, U(M) is the uniform law of a basis B, and x_e = 1[e ∈ B]. A split is E = X ⊔ Y, balanced if n/3 ≤ |X| ≤ 2n/3. β(U(M); X|Y) is the paper's affinity (`sec-tools`): ‖√μ‖_op, with rows indexed by B∩X and columns by B∩Y.

* **Proposition 1 (duality) [PROVED].** β(U(M); X|Y) = β(U(M*); X|Y) for every split. A PC for U(M) becomes a PC for U(M*) of the same size, smoothness, decomposability, leaf arity and TV by flipping every leaf input.
* **Theorem A (one split, every matroid) [PROVED].** Let K₁,…,K_r be pairwise compatible framed pieces (§1) that meet both X and Y. Then
  β(U(M); X|Y)² ≤ ∏ᵢ (1 − w(Kᵢ)),  with w(K) = ∏_{e∈K} |D_e|⁻¹.
* **Theorem C (PC lower bound, every matroid) [PROVED].** Let 𝒦 be a family of framed pieces of weight ≥ w₀, in which each piece is incompatible with at most Γ others. Suppose every balanced X is crossed by at least m pieces of 𝒦. Then every nonnegative smooth decomposable PC over the n element variables, with leaf arity ≤ 2n/3 and 1 − TV(U(M), Q_C) ≥ δ, has
  |C| ≥ δ² (1 − w₀)^{−⌈m/(Γ+1)⌉} ≥ δ² exp(w₀ m/(Γ+1)).
  With m ≥ c·n for fixed c > 0, w₀ ≥ c₀ > 0 for fixed c₀, and Γ = O(1), this is e^{Ω(n)}. The hypothesis is self-dual.
* **Theorem B (global version, balanced matroids) [PROVED].** Frames are not needed. Each factor 1/|D_e| is replaced by a leverage in a worst-case minor: λ(e; M/F′) for circuit pieces, 1 − λ(e; M\F′) for cocircuit pieces. For **regular** matroids, 1 − λ(e; M\Z) = 1/(1 + R_Z(e)) with an exact Thomson/resistance formula (Lemma 7).
* **Instances.**
  * Graphic: v1 with a **sharper constant**, D^{−ℓ} in place of D^{−2(ℓ−1)} [PROVED].
  * R2 recovered as the cocircuit (vertex-star) case of Theorem B [PROVED reduction].
  * Cographic: free by duality [PROVED].
  * A new non-regular binary family, linked expanders with every edge replaced by a Fano plane (2-sum) [PROVED].
  * Uniform matroids show that high Tutte connectivity alone is not enough [PROVED].
* **Bonus: the whole argument is circuit/cocircuit duality (§9).** v1's cycle bits are circuit pieces of M(G) framed by cuts (vertex stars). R2's vertex-touch bits are cocircuit pieces of M(G). The theorem is invariant under M ↔ M*. The only combinatorial input is orthogonality, |C ∩ D| ≠ 1.

---

## 1. Definitions

* A **cut** of M is a set D such that E∖D is a flat, or equivalently a union of cocircuits.
* A **cycle** is a set D such that E∖D is a flat of M*, or equivalently a union of circuits.
* In a graph, cuts here are unions of bonds (not necessarily a single δ(U)), and cycles are unions of ordinary cycles (not necessarily even subgraphs). Binary cycle/cocycle spaces use symmetric difference and are narrower notions. Vertex stars remain valid cuts.

**Framed piece.** A framed piece is a triple (K, π, (D_e)_{e∈K}) with the following properties:
* π is a cyclic permutation of K that forms a single cycle;
* either (circuit type) K is a circuit and every D_e is a cut, or (cocircuit type) K is a cocircuit and every D_e is a cycle;
* D_e ∩ K = {e, π(e)} for every e ∈ K.

The **weight** is w(K) = ∏_{e∈K} |D_e|⁻¹ and the **territory** is T(K) = K ∪ ⋃_e D_e. Two pieces K, K′ are **compatible** if K′ ∩ T(K) = ∅ and K ∩ T(K′) = ∅. A piece is **mixed** for X|Y if it meets both X and Y.

**Every circuit of size at least two has a frame [PROVED].** Take any cyclic order and set D_e := E ∖ cl(K∖{e, π(e)}).
* cl(·) is a flat, so D_e is a cut.
* K∖{e,π(e)} lies inside the flat.
* e ∉ cl(K∖{e,π(e)}), because K∖π(e) is independent. The same holds for π(e).

So D_e ∩ K = {e, π(e)}. Dually, every cocircuit of size at least two has a frame of cycles. Singleton loops and coloops have no frames, but singletons cannot cross a split. The content of the hypotheses is that the frames are **small** and **avoid the other pieces**.

**Graphs.**
* A cycle C = v₁e₁v₂e₂⋯ of G, with π(e_j) = e_{j+1} and D_{e_j} = δ(v_{j+1}), is a circuit-type piece with weight ∏_{v∈V(C)} deg(v)⁻¹. The star is a cut even when G − v is disconnected.
* T(C) is the set of edges touching V(C), so two cycles are compatible iff they are vertex-disjoint. This is exactly v1's hypothesis.

## 2. Imported tools (distribution-free)

All of these are in `notes/paper-ust/sec-tools.tex`. None uses graph structure. They concern a distribution on {0,1}^E and a split.
* lem:dpi: data processing under maps of each side.
* lem:chain: adaptive chain rule. The kernels may depend on the joint past, and the pieces and their order are fixed in advance (rem:adaptive).
* lem:two-by-two: if K(1,1)=0 then ‖√K‖² = (1+√(1−4K(0,1)K(1,0)))/2 ≤ 1 − K(0,1)K(1,0). The last step is the paper's √(1−x) ≤ 1−x/2.
* cor:reduction: |C| ≥ δ²/max_balanced β². Its balanced splits satisfy |A| > n/3 and |E∖A| ≥ n/3, which are covered by our "n/3 ≤ |X| ≤ 2n/3".

## 3. Duality

**Proposition 1 [PROVED].** B ↦ E∖B is a bijection 𝓑(M) → 𝓑(M*). It sends (B∩X, B∩Y) to (X∖B, Y∖B), a bijection on each side separately. So the matrices √U(M) and √U(M*) agree up to a permutation of rows and a permutation of columns, and β is equal. For PCs, substitute x_e ↦ 1−x_e at the leaves. Scopes are unchanged, so smoothness, decomposability, size and arity are unchanged, and TV is preserved because the substitution is a bijection of {0,1}^E. ∎

## 4. Local conditioning

**Lemma 2 (conditioning is a minor) [PROVED].** Let F ⊆ E and σ ⊆ F with Pr[B∩F = σ] > 0, and put τ = F∖σ. Then {B∖F : B ∈ 𝓑(M), B∩F = σ} = 𝓑(M/σ\τ). Hence, given B∩F = σ, B∖F is a uniform basis of M/σ\τ.

*Proof.* Some basis B₀ satisfies B₀ ⊆ (E∖F) ∪ σ. So r_M((E∖F)∪σ) = r, and r(M/σ\τ) = r − |σ|, with σ independent.
* If B∩F = σ, then B∖F ∪ σ = B is independent, so B∖F is independent in M/σ and has r − |σ| elements. It is therefore a basis of the minor.
* Conversely, for a basis B′ of the minor, B′ ∪ σ is independent in M with r elements, so it is a basis, and it meets F in σ. ∎

**Lemma 3 (exchange bound, every matroid) [PROVED].** If D is a cut of N and e ∈ D, then Pr_{U(N)}[e ∈ B] ≥ 1/|D|. Dually, if C is a cycle of N and e ∈ C, then Pr[e ∉ B] ≥ 1/|C|.

*Proof.* First shrink D to a cocircuit D′ ∋ e with D′ ⊆ D: extend the flat E∖D to a hyperplane H with e ∉ H and take D′ = E∖H. If e is a loop of N, then e lies in no cocircuit, a contradiction, so e is not a loop.

Let B ∌ e. The fundamental circuit C(e,B) ⊆ B+e contains e ∈ D′, and |C ∩ D′| ≠ 1 (orthogonality). So there is f ∈ C ∩ D′ with f ≠ e, f ∈ B, and B − f + e is a basis. Fix one such f for each B.

The map B ↦ B − f + e sends bases without e to bases with e. Each image B′ has at most |D′| − 1 preimages, namely B′ − e + f with f ∈ D′∖e. So #{B ∌ e} ≤ (|D′|−1)·#{B ∋ e}, which gives Pr[e∈B] ≥ 1/|D′| ≥ 1/|D|. The dual statement is the same claim in N*. ∎

In graphs this gives Pr[e ∈ T] ≥ 1/deg(v) and Pr[e ∉ T] ≥ 1/|C|, which is R_eff(e) ≤ (|C|−1)/|C|. We did not find Lemma 3 stated for matroids. It is folklore-level. Douthitt–Oxley's basis-count bound ∏|C*_i| has the same flavour.

**Lemma 4 (certificates survive minors) [PROVED].** Let D be a cut of M, e ∈ D, and σ, τ ⊆ E∖{e} disjoint with D ∩ σ = ∅. Put N = M/σ\τ. Then Pr_{U(N)}[e ∈ B] ≥ 1/|D∖τ| ≥ 1/|D|.

*Proof.* Let H = E∖D, a flat containing σ, and Z = H∖(σ∪τ). For any x in the ground set of N, x ∈ cl_N(Z) iff x ∈ cl_M(Z∪σ). Since Z∪σ ⊆ H and e ∉ H, we get e ∉ cl_N(Z). So (E∖(σ∪τ))∖cl_N(Z) is a cut of N that contains e and lies inside D∖(σ∪τ) = D∖τ. Apply Lemma 3. ∎

**Lemma 5 (frame bound) [PROVED].** Let (K, π, D) be a circuit-type framed piece, S ⊊ K, and F ⊆ E∖K with D_e ∩ F = ∅ for every e ∈ S. Then, for every σ with Pr[B∩F=σ] > 0,
  Pr[S ⊆ B | B∩F = σ] ≥ ∏_{e∈S} |D_e|⁻¹.
For a cocircuit-type piece with cycle frames, the same holds for Pr[S ∩ B = ∅ | B∩F = σ].

*Proof.* Since S ≠ K and π is one cycle, S splits into maximal runs e, π(e), …, π^k(e), and no run wraps around. List the runs one after another, each in π-order, as s₁, …, s_k.
* D_{s_j} ∩ S ⊆ {s_j, π(s_j)}, and π(s_j) is either outside S or later in the list.
* So D_{s_j} avoids F ∪ {s₁,…,s_{j−1}}.

By induction on j, the event {B∩F=σ, s_{<j} ⊆ B} has positive probability. On it, B ∩ (F∪s_{<j}) = σ∪s_{<j}, and by Lemma 2 the remaining law is U(M/(σ∪s_{<j})\τ). Lemma 4, with σ∪s_{<j} in the role of σ, gives a factor ≥ 1/|D_{s_j}|, which is positive, and this keeps the induction going. Multiply the factors.

The cocircuit case is Lemma 5 for M* applied to E∖B (Proposition 1). There B∩F=σ becomes (E∖B)∩F = F∖σ, and S∩B=∅ becomes S ⊆ E∖B. ∎

## 5. Theorem A

**Theorem A [PROVED].** Let K₁,…,K_r be pairwise compatible framed pieces, each mixed for X|Y. Then β(U(M);X|Y)² ≤ ∏ᵢ (1 − w(Kᵢ)).

*Proof.*

*Bits.* For circuit type, set uᵢ = 1[Kᵢ∩X ⊆ B] and vᵢ = 1[Kᵢ∩Y ⊆ B]. For cocircuit type, set uᵢ = 1[Kᵢ∩X∩B = ∅] and vᵢ = 1[Kᵢ∩Y∩B = ∅].
* uᵢ = vᵢ = 1 is impossible: a basis contains no circuit and meets every cocircuit.
* uᵢ depends only on B∩X and vᵢ only on B∩Y, so lem:dpi gives β ≤ ‖√law(u,v)‖.

*Chain.* Put Fᵢ = K₁∪⋯∪K_{i−1}. The past (u_{<i}, v_{<i}) is a function of B∩Fᵢ. By compatibility, every D_e with e ∈ Kᵢ avoids Fᵢ. Lemma 5, with S = Kᵢ∩X ⊊ Kᵢ, gives Pr[uᵢ=1 | B∩Fᵢ=σ] ≥ pᵢ := ∏_{e∈Kᵢ∩X}|D_e|⁻¹ for every σ. Averaging over σ given any positive-probability past π gives K(1,0|π) ≥ pᵢ, and in the same way K(0,1|π) ≥ qᵢ := ∏_{e∈Kᵢ∩Y}|D_e|⁻¹.

*Two-by-two.* lem:two-by-two gives γᵢ² ≤ 1 − pᵢqᵢ = 1 − w(Kᵢ). Then lem:chain finishes the proof. ∎

## 6. Theorem C and necessity

*Proof of Theorem C [PROVED].* Fix a balanced X. The incompatibility graph on the ≥ m mixed pieces has maximum degree ≤ Γ, so greedy selection gives ⌈m/(Γ+1)⌉ pairwise compatible mixed pieces. Theorem A bounds β², and cor:reduction gives the PC bound. The last inequality is −log(1−w₀) ≥ w₀. Self-duality: (M, 𝒦) ↦ (M*, 𝒦) swaps the piece types and the cut/cycle roles of the frames, and keeps w, T, compatibility and crossing counts. ∎

**Intrinsic version [PROVED from C].** Call M (ℓ, d, Δ)-local if the following hold:
* 𝒦 consists of pieces with |K| ≤ ℓ and |D_e| ≤ d;
* every element lies in at most Δ pieces and in at most Δ territories.

Then Γ ≤ ℓ(d+1)Δ + ℓΔ and w₀ ≥ d^{−ℓ}. "Expansion" means that the hypergraph 𝒦 on E has edge expansion: every balanced X is crossed by ≥ c·n pieces. This gives |C| ≥ δ² exp(c n d^{−ℓ}/(ℓ(d+2)Δ + 1)).

**Why some hypothesis of this kind is needed [PROVED].**
1. *Separations.* If M = M₁ ⊕ M₂ with E(M₁) balanced, then U(M) = U(M₁) ⊗ U(M₂) and β = 1. This is the matroid form of prop:pendant: coloops are the extreme case. No circuit or cocircuit crosses a separator, so the hypothesis of Theorem C fails correctly.
2. *Connectivity is not enough.* U_{k,2k} has the maximum connectivity function λ(X) = k for |X| = k, yet β(U; X|Y)² = max_j Pr[|B∩X| = j] = Θ(k^{−1/2}) [PROVED: √μ is a direct sum of all-ones blocks; CHECKED k ≤ 7]. It also has an exact PC of size O(n k²): the count dynamic program on a balanced tree is smooth and decomposable. Its circuits have size k+1, so pieces have weight ~k^{−k} and Theorem C correctly gives nothing.

    So the right "expansion" is crossing by *short framed* pieces, a locality condition. It is not Tutte connectivity, and not edge expansion of some graph on E.

## 7. Instances

**7.1 Graphic: v1 with a sharper constant [PROVED].** Take the hypotheses of thm:general: max degree D, edge expansion h, linked by 𝒞 with parameters (ℓ, κ), and r₀ as in the paper.
* lem:mixed (`sec-geometry`, with m₀ = M/3) gives r₀ vertex-disjoint mixed cycles of length ≤ ℓ.
* These are compatible star-framed pieces (§1). Theorem A gives
  β(UST(G);A|B)² ≤ ∏_k (1 − ∏_{v∈V(C_k)} deg(v)⁻¹) ≤ (1 − D^{−ℓ})^{r₀},  hence |C| ≥ δ² exp(r₀ D^{−ℓ}).

The paper has (1 − D^{−2(ℓ−1)})^{r₀}. The new exponent constant is larger by a factor D^{ℓ−2}. For cubic prisms (D=3, ℓ=4), the constant goes from 3^{−6} to 3^{−4}, which is 9× larger.

The gain comes from bounding the *product* K(0,1)K(1,0) = pq directly. Each edge of the cycle is charged to its forward vertex, so the A-half and the B-half together charge every vertex once. v1 bounds p and q separately by the same ∏_{y≠x₀}, which charges almost every vertex twice. The paper's per-cycle p_k can still be better for very irregular cycles, for example one vertex of huge degree with all others of degree 2, so take the max of both.

v1's lem:contract and lem:local are not needed: conditioning on the previous pieces alone suffices. [CHECKED: chain ≤ this bound on all tested splits of K₄, K₅, K₃,₃, Q₃, Petersen and K₄□K₂.]

**7.2 Cographic [PROVED].** U(M*(G)) is the complement of a UST. By Proposition 1, every graph covered by v1 or R2 gives a cographic matroid with the same β and the same PC bounds. Here cycles of G are cocircuit pieces of M*(G) and stars are their cycle-frames. M*(G) is not graphic when G is nonplanar, so these are formally new matroids, but the result is trivially equivalent to the graphic one.

**7.3 Fano-decorated linked expanders: binary, non-regular [PROVED].** Let G satisfy the hypotheses of thm:general, with M edges. Let M_G be the 2-sum of M(G) with one copy F⁽ᵉ⁾ of F₇ at every edge e, where e is identified with a point p_e of F₇ and then deleted.

*Properties.* M_G has n = 6M elements. It is binary (2-sums of binary matroids are binary), has F₇ minors, and so is non-regular, non-graphic and non-cographic.

*Structure used.* These are standard 2-sum facts (Oxley §7.1). The cocircuit fact is the circuit fact for the dual, since (M⊕₂N)* = M*⊕₂N*.
* (a) Circuits or cocircuits of F⁽ᵉ⁾ avoiding p_e remain circuits or cocircuits of M_G.
* (b) If C is a cycle of G and {p_e, a_e, b_e} is a line of F⁽ᵉ⁾ for each e ∈ C, then ⋃_{e∈C}{a_e, b_e} is a circuit.
* (c) A bond β ⊆ δ(v) of G lifts to the cocircuit ⋃_{f∈β}(D_f − p_f), for any cocircuits D_f ∋ p_f of F⁽ᶠ⁾.

*Local pieces.* Each of the 4 lines of F₇ that avoid p gives a 3-element circuit. Its frame is D_x = F₇ ∖ (the line through p and the third point), a 4-element cocircuit avoiding p. These pieces have weight 4⁻³ and territory inside one gadget.

If the six non-base points of a gadget are split nontrivially, some p-avoiding line is mixed. Two lines of F₇ meet in exactly one point, and that point is not p for two p-avoiding lines. So if all four were monochromatic, any two of them would share a point and have the same colour, and the four lines cover all six points, so the whole gadget would have one colour.

*Lifted cycles.* Take the cyclic order a₁, b₁, a₂, b₂, … on a cycle of G. The frames are:
* D(a_e) = F⁽ᵉ⁾ minus another line through p (size 4);
* D(b_e) = the lifted star at the vertex shared by e and the next edge e′, choosing D_e ∋ b_e, ∌ a_e and D_{e′} ∋ a_{e′}, ∌ b_{e′}. Such cocircuits exist: take the complement of the line through a_e that avoids p.

So |D| ≤ 3D and the weight is ≥ (3D)^{−2ℓ}. Vertex-disjoint cycles of G give compatible lifts, because no territory leaves the edges at V(C).

*Counting.* Let ε := ½·min{hN/(4Dℓ²κM), 1/(2ℓ²κD(1+D/h)), 1/12}. Fix a balanced split of E(M_G) and let c be the number of cut gadgets.
* If c ≥ εM, there are εM compatible local pieces, one per cut gadget.
* If c < εM, colour each edge of G by its gadget's colour, with cut gadgets coloured arbitrarily. Each colour class has ≥ M/3 − εM ≥ M/4 edges. lem:mixed with m₀ = M/4 gives ≥ t/(ℓκ) ≥ 2εM vertex-disjoint mixed cycles. At most c of them contain a cut gadget. Each remaining cycle lifts to a mixed piece, because its X-gadget and its Y-gadget are both whole.

Either way r ≥ εM, and Theorem A plus cor:reduction give |C| ≥ δ² exp(εM·(3D)^{−2ℓ}).

*Checked.* On F₇-decorated K₃ (n = 18, 9216 bases), all local pieces, the lifted triangle and every frame were verified, and Theorem A held on random splits [CHECKED, `f7deco.py`]. The counting step was not run at scale.

The same works with any gadget N where "every nontrivial split of N−p has a mixed framed local piece". Decorating with S8 would give a non-balanced family, but its gadget property was not checked. [CONJ]

**7.4 General matroids.** Theorems A and C use only orthogonality and basis exchange. So **neither regularity, binary representability, strong Rayleigh nor Lorentzian structure is needed for the local version.** The answer to "binary non-regular and general matroids" is: yes, for the local version, with no negative dependence at all. [PROVED; Lemma 3 and the bounds CHECKED on F₇, F₇*, AG(3,2), S8 and U_{2,4}, which are not regular.]

## 8. Global version: leverage replaces effective resistance

Write λ(e; N) := Pr_{U(N)}[e ∈ B] for the leverage. M is **balanced** (Feder–Mihail) if U(N) is negatively correlated for every minor N. This class is closed under duality, because NC for complements is equivalent to NC. It contains:
* regular matroids (Feder–Mihail 1992);
* binary matroids with no S8 minor (Choe–Wagner: these are exactly the binary Rayleigh matroids, and Rayleigh implies balanced).

**Lemma 6 (worst case is "contract everything" / "delete everything") [PROVED for balanced M].** Let F ⊆ E, Pr[B∩F=σ] > 0, and S = {s₁,…,s_k} ⊆ E∖F in any order. Then
  Pr[S ⊆ B | B∩F=σ] ≥ ∏_j λ(s_j; M/F/s_{<j}),   Pr[S∩B=∅ | B∩F=σ] ≥ ∏_j (1 − λ(s_j; M\F\s_{<j})).

*Proof.* Put τ = F∖σ.

*Factorization.* By Lemma 2 and induction, while the factors stay positive, Pr[S⊆B|σ] = ∏_j λ(s_j; N_j) with N_j = M/σ/s_{<j}\τ.

*Comparison.* Pass from N_j to M/F/s_{<j} by turning the elements f ∈ τ from deleted into contracted, one at a time. Each step compares N′\f with N′/f for a minor N′ of M.
* If f is a loop or coloop of N′, the two minors coincide.
* Otherwise λ(e;N′/f) = Pr_{N′}[e|f∈B] ≤ Pr_{N′}[e|f∉B] = λ(e;N′\f) by NC of N′.

So λ(s_j; N_j) ≥ λ(s_j; M/F/s_{<j}).

*Positivity.* If the right-hand product is 0 there is nothing to prove. Otherwise every left factor is positive, which justifies the factorization.

*Second inequality.* Apply the first in M* (Proposition 1, with λ(e;N*) = 1 − λ(e;N)). ∎

Only NC of single pairs in minors is used, not negative association. This matches R2, which uses Rayleigh monotonicity rather than negative association. An NA-based proof also works for strongly Rayleigh measures and arbitrary increasing events, but it is not needed.

**Lemma 7 (regular: leverage = resistance) [PROVED].** Let A be a totally unimodular real matrix of full row rank representing a regular matroid M. By Cauchy–Binet, U(M) is the projection DPP with kernel Aᵀ(AAᵀ)⁻¹A (Lyons 2003), so λ(e; M) = a_eᵀ(AAᵀ)⁻¹a_e. For Z ⊆ E∖e define
  R_Z(e) := min{ ‖x‖² : x ∈ ℝ^{E∖(Z∪e)}, Σ_f x_f a_f = a_e }   (∞ if infeasible).
Then 1 − λ(e; M\Z) = 1/(1 + R_Z(e)).

*Proof.* Let A′ be the columns E∖(Z∪e), restricted to a row basis of the columns E∖Z.
* If e is a coloop of M\Z, both sides are 0.
* Otherwise A′ has full row rank. Sherman–Morrison gives λ = R/(1+R) with R = a_eᵀ(A′A′ᵀ)⁻¹a_e, and least-norm/Thomson gives R = R_Z(e). ∎

Any explicit x gives an upper bound on R, so a lower bound on 1 − λ. For example, m signed circuits through e, each of size ≤ ℓ, that avoid Z and pairwise share only e give R ≤ (ℓ−1)/m. More generally, R ≤ ∑_j (|C_j|−1)/m². This is the "fractional circuit flow" form of Thomson's principle. For graphs, R_Z(e) = R_{G−Z−e}(u,v). Dually, λ(e; M/Z) = 1/(1 + R*_Z(e)), computed in M*; for graphs this is R_eff^{G/Z}(e). [CHECKED to 3·10⁻¹⁵ on R₁₀ and R₁₂ with signed TU representations.]

**Theorem B [PROVED for balanced M].** Let K₁,…,K_r be pairwise disjoint mixed circuits or cocircuits, in a fixed order; no frames are needed. Put Fᵢ = K₁∪⋯∪K_{i−1}, and order each half H of Kᵢ as h₁,…,h_k. Define
* φ(h_j) = λ(h_j; M/(Fᵢ∪h_{<j})) for a circuit piece;
* φ(h_j) = 1 − λ(h_j; M\(Fᵢ∪h_{<j})) for a cocircuit piece.

Then β² ≤ ∏ᵢ (1 − ∏_{e∈Kᵢ} φ(e)). For regular M, the factors are 1/(1+R*) and 1/(1+R) from Lemma 7.

*Proof.* The proof of Theorem A, with Lemma 6 in place of Lemma 5. ∎

**R2 is the cocircuit case [PROVED reduction].** In a 2-connected G, the star S(w) of each w in an independent set W of mixed vertices is a bond, so a cocircuit piece, and these pieces are disjoint. For s = uw ∈ S_A(w), Theorem B gives the factor 1/(1 + R_{G−Fᵢ−s_{≤j}}(u,w)) ≥ 1/(1 + R_u), where R_u = R_{G−(W∖w)−S_A(w)}(u,w), by Rayleigh: the second graph has fewer edges. That is exactly R2's Lemma 2, p_A(w) = ∏_{u∈N_A(w)}(1+R_u)⁻¹. R2's Lemmas 3–4 (expansion of G − W) are the graph-specific input that bounds R_u. They have no general-matroid analogue here; for a regular matroid, "R_Z(e) ≤ ρ for the relevant Z" must be assumed.

**Where the global route fails [CHECKED + reason].**
* *Lemma 6 is false without NC.* On S8, the 8-element binary matroid with 48 bases and a positively correlated pair, exhaustive search finds 4 of 1624 cases with |F| ≤ 2 and |S| ≤ 2 where min_σ Pr[S⊆B|σ] < Pr_{M/F}[S⊆B]; the worst ratio is 0.933. It held on all the balanced test matroids, including F₇, F₇* and AG(3,2).
* *Lorentzian structure does not rescue it.* Huh–Schröter–Wang give Pr[e,f] ≤ 2Pr[e]Pr[f] for all matroids, which costs a constant factor per element flipped in Lemma 6. Since |F| = Θ(n), the loss compounds to e^{−Θ(n)} and the bound is useless.

So for non-balanced matroids, only the local Theorem A/C is available. A global bound would need a "local influence" version of Lemma 6 that is not known for matroids. [OPEN]

## 9. Bonus: the argument is circuit/cocircuit duality

1. **The forbidden cell.** A basis contains no circuit and meets every cocircuit; its complementary cobasis contains no cocircuit and meets every circuit. These are the dual facts read in M and in M*. The cell (1,1) is forbidden for circuit pieces of M exactly when it is forbidden for the corresponding cocircuit pieces of M*.
2. **The conditional bound.** Lemma 3 uses orthogonality once: the fundamental circuit meets D a second time. Frames use it again: |D_e ∩ K| = 2 is the smallest value orthogonality allows. A circuit is certified by cocircuits, and a cocircuit by circuits.
3. **The theorem.** Theorems A and C are invariant under (M, circuit pieces with cut frames) ↔ (M*, cocircuit pieces with cycle frames), and so is the PC lower bound (Proposition 1).
4. **The two graph papers are the two halves.** v1 uses circuits of M(G), certified by cuts (vertex stars). R2 uses cocircuits of M(G) (vertex stars), certified globally by leverage (resistance). Each paper's pieces are the other paper's certificates.

    The local version of R2 is the dual half of Theorem A: stars framed by short cycles through consecutive edges at w. This makes R2's grid sketch (`R2-highgirth.md`, "local resistance bounds") rigorous. A 3-separated W in the grid gives unit squares as frames, of weight 4⁻⁴. This adds nothing beyond the paper's cor:grid.

## 10. Weak points and gaps

* **Reviewer target 1: Lemma 4.** This is the step that turns a cut of M avoiding σ into a cut of the random minor M/σ\τ. The proof is two lines, but it is the step every conditional bound rests on. The claim "x ∈ cl_N(Z) iff x ∈ cl_M(Z∪σ)" for N = M/σ\τ is the standard contraction closure. A failure would show up as an exchange-bound violation; none appeared in about 57k exact cases.
* **Reviewer target 2: §7.3.** It relies on cited 2-sum structure facts (a)–(c) and on lem:mixed at m₀ = M/4 (the lemma allows any m₀ > 0). The structure was verified only on F₇-decorated K₃.
* **Order and adaptivity.** Pieces and their order are fixed before sampling, as rem:adaptive requires. Greedy selection in Theorem C depends only on X.
* **Theorem B for regular matroids has no intrinsic "expansion ⇒ bounded R" lemma.** Only graphs have one (R2 Lemma 4). A spectral or Seymour-decomposition version is [OPEN].
* **No new exponential family beyond graphic/cographic for the global version.** §7.3 is local.

## 11. Numerics (exact enumeration, `toy/r3/matroid/`)

All runs are exact enumeration of bases, with probabilities in floating point. The runner `run_checks.sh` (6 workers) exits nonzero if any worker fails, and each script exits nonzero if a PROVED inequality is violated. All exited 0. Logs are in `toy/r3/matroid/logs/`.

| check | matroids | result |
|---|---|---|
| Prop 1: β(M) = β(M*) | 15 matroids (graphic K₄, K₅, K₃,₃, Q₃, Petersen, K₄□K₂; cographic duals; R₁₀, R₁₂, F₇, F₇*, AG(3,2), S8, U₂,₄) | max diff 6·10⁻¹⁶ |
| Lemma 3: min_e min_{D∋e} Pr[e∈B]·\|D\| ≥ 1 (M and M*) | same | min 1.5 |
| Lemma 5 (random S ⊊ K, random F) | same, M and M* | 0 violations / 56,944 cases |
| Lemma 6, random | same | 0 violations, S8 included (random F missed the bad cases) |
| Lemma 6, exhaustive \|F\|,\|S\| ≤ 2 | 11 matroids | 0 violations except **S8: 4/1624, ratio 0.933** |
| β ≤ ‖√μ_bits‖ ≤ ∏γ_k ≤ (Thm B bound) ≤ (Thm A bound) | 180 (split, piece-family) pairs over 15 matroids | all hold |
| Lemma 7, resistance formula | R₁₀, R₁₂ (TU signings found by brute force) | error ≤ 3.4·10⁻¹⁵ |
| U_{k,2k}: β² = max_j Pr[\|B∩X\|=j] | k = 3…7 | exact |
| §7.3 structure + Thm A | F₇-decorated K₃ | all hold |

Typical sizes: on Petersen with four 3-star pieces, β = 0.43, ∏γ = 0.953, Theorem B bound 0.953, Theorem A bound 0.988. On K₃,₃ with three triangles in the dual, the bits norm is 0.82. The bounds hold everywhere, but the local bound is loose, roughly 1 − d^{−ℓ} versus the true kernel factor.

Library notes: R₁₀ and R₁₂ were confirmed regular by exhibiting TU signings (`tu_check.py`). F₇ and AG(3,2) have no TU signing. "S8" is [I₄|J−I+e₄₄], with 48 bases and a positively correlated pair.

## 12. Literature (searched Nov 4 2026; deep and fast rounds)

Closest work on circuit lower bounds for basis-generating objects:
* arXiv 2109.06941, *Monotone complexity of spanning tree polynomial re-visited*: a 2^{Ω(n)} exact monotone-arithmetic lower bound for the **directed rooted spanning-tree polynomial** on bidirected constant-degree expanders. Identifying the two directed variables gives the undirected graphic basis-generating polynomial, but does not automatically preserve the lower bound. Its discrepancy and epsilon-sensitive results are relevant prior art, not the TV/affinity theorem proved here or a result for general matroids.
* ECCC TR20-191, *Negations provide strongly exponential savings*; ECCC TR18-049 (K_n).
* Fomin–Grigoriev–Koshevoy, arXiv 1307.8425 (subtraction-free complexity).
* **arXiv 2512.19515, *Negations are powerful even in small depth* (Dec 2025):** an explicit matroid whose monotone Boolean function needs superpolynomial monotone circuits, plus exponential monotone vs non-monotone arithmetic separations. **Closest scoop risk.** It is exact computation of a specific matroid's function, not distributional approximation of U(M) under a general hypothesis. The v2 related-work section should cite it and say how this differs.
* None of these give TV/affinity lower bounds for U(M) of general matroids.

Distributional side:
* Mészáros–Peled, arXiv 2610.11554 (Oct 2026), *Product approximations for USTs and determinantal processes*: USTs (and DPPs including Kalai hypertrees) are close to product-like 1-out models when the minimum degree k is large, O(n log k/√k) edits. This is consistent with our setting: high degree makes w₀ = d^{−ℓ} tiny, and our bounds degrade exactly there.
* Lapidoth–Pfister (the β tool); the paper's own PC references.

Negative dependence and matroids:
* Feder–Mihail 1992 (balanced matroids; regular ⇒ balanced).
* Choe–Wagner, *Rayleigh matroids*, arXiv math/0307096 (binary Rayleigh ⇔ no S8 minor; Rayleigh ⊂ balanced).
* Jerrum–Son, arXiv math/0404200 (a non-balanced paving matroid).
* Borcea–Brändén–Liggett (strong Rayleigh ⇒ negative association).
* Lyons 2003 (determinantal measures; regular matroids).
* Huh–Schröter–Wang, arXiv 1806.02675 (factor-2 correlation bound for all matroids).
* Anari–Oveis Gharan–Vinzant, arXiv 1807.00929, and Brändén–Huh (Lorentzian).
* Fu, arXiv 2610.11829 (bases vs dimension in binary matroids; background only).
* Douthitt–Oxley, *Bounding the number of bases of a matroid* (∏|C*_i|; same flavour as Lemma 3).
* Chandrasekaran–Chekuri–Zhu, arXiv 2610.33918 (random deletion for circuits in regular matroids, Oct 2026): unrelated technique, no overlap.

**Novelty status.** We found no prior result giving TV, affinity or PC lower bounds for the uniform basis of any non-graphic matroid family, and no circuit/cocircuit-framed chain argument. Lemma 3 is probably folklore. The search was not exhaustive, and one engine was rate-limited (CAPTCHA) on the last rounds.

## 13. Files

* `notes/approaches/R3-matroid.md` (this note).
* `toy/r3/matroid/matroid.py`: matroid toolkit; imports, and does not modify, `toy/r2/common.py`.
* `toy/r3/matroid/tu_check.py`: TU-signing search.
* `toy/r3/matroid/checks.py` and `toy/r3/matroid/run_checks.sh`: duality, Lemmas 3/5/6, split-level chain.
* `toy/r3/matroid/extra.py`: S8, exhaustive Lemma 6, resistance formula, uniform matroids.
* `toy/r3/matroid/f7deco.py`: §7.3.
* `toy/r3/matroid/logs/`: all outputs.

Reproduce: `cd toy/r3/matroid && ./run_checks.sh && uv run --with numpy --with scipy --with networkx python extra.py && uv run --with numpy --with scipy --with networkx python f7deco.py`.

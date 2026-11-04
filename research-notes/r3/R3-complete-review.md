# Adversarial review of R3: complete-graph spanning trees

Reviewed November 4, 2026 (America/Chicago).

## Verdict: HOLDS-WITH-FIXES

**Theorem K survives.** I found no mathematical failure in the bound on product affinity, including its uniformity over arbitrary balanced edge splits. The conditioning argument gives the bounds for every feasible detailed conditioning, which is enough for every positive-probability past in the Hellinger chain. Proposition U also survives, including smoothness over all edge variables.

**Corollary K needs its leaf-arity hypothesis in the statement.** Literally unrestricted leaves give a one-node counterexample. The note acknowledges this exclusion later, so I regard it as a repairable statement defect rather than a failure of the intended result. The reported numerical size of the upper circuit is also not correctly counted by the supplied counter.

- FATAL: none in Theorem K or the intended bounded-arity circuit result.
- MAJOR: F1, missing leaf-arity restriction in Corollary K.
- MINOR: F2, upper-circuit size counter; F3, unsupported inference from edit distance to TV in the literature table.
- Score for the result as stated: **8/10**. With F1 repaired, the main proof is convincing. I found no prior fixed-TV result that makes this claim redundant, but the literature check was a search, not a novelty certificate.

I did not modify the claimed proof, its programs, or anything under `/home/chris/code/eth-research/notes/paper-ust/`.

## Findings and concrete fixes

### F1. MAJOR: Corollary K does not state the arity restriction

The opening statement quantifies over every smooth decomposable nonnegative PC. Under the model in `/home/chris/code/eth-research/notes/paper-ust/sec-prelim.tex`, leaves may be nonnegative functions on several variables, with arity a. The cited balanced-decomposition lemma in `/home/chris/code/eth-research/notes/paper-ust/sec-tools.tex` requires **a <= 2M/3**. Its circuit-reduction corollary explicitly inherits this hypothesis.

Without any restriction, take one leaf whose function is the UST probability mass function on all M variables. Smoothness, decomposability and determinism are vacuous; its size is one. For n = 3600 and delta = 1, the displayed lower bound would require size at least e/2 > 1. This is a concrete counterexample to the unrestricted wording.

The note's later “What is not covered” paragraph recognizes the exclusion above 2M/3. That repairs the reader's understanding, but it should not be separated from the theorem's quantifier.

**Fix:** state Corollary K for leaves of arity at most 2M/3, with size measured as nodes plus wires, as in the cited lemma. Alternatively state it for univariate-leaf PCs. No change to Theorem K, Proposition U, or the exponential conclusion for the usual PC model is needed.

### F2. MINOR: the upper-circuit counter omits a multiplicity

In `/home/chris/code/eth-research/toy/r3/complete/upper.py`, `children` includes the t choices of u in S', but `cut_lits` does not:

```python
children = sum(comb(s - 2, t - 1) * t for t in range(1, s))
cut_lits = sum(comb(s - 2, t - 1) * (t * (s - t) - 1) for t in range(1, s))
```

The displayed construction has a separate summand for each u. Its zero-literal set also depends on u, because the omitted cut edge is vu. Without a specified sharing construction, the cut-literal wires need the same factor t. This does not count the described circuit's size correctly, even though the proof's asymptotic O(n^4 3^n) bound is correct.

For a concrete implementation with flattened product nodes, shared univariate indicator leaves, and one shared constant for singleton states, the count is

```
2M + 1 + sum_{s=2}^n C(n,s) s
             [1 + sum_{t=1}^{s-1} C(s-2,t-1) t (4 + t(s-t))].
```

At n = 12 this is **95,300,245**, versus the counter's **38,524,500**. The difference is not a lower bound on all implementations: additional sharing could reduce the count. It is a correction for the straightforward implementation being counted.

**Fix:** include t in `cut_lits` and count the shared leaves, or specify an actual shared cut-product construction and count its nodes and wires. Keep Proposition U's asymptotic statement unchanged.

### F3. MINOR: expected edit distance does not establish a TV lower bound

The Mészáros–Peled row says the 1-out approximation has TV approximately one “since” many edges differ. A bound on expected symmetric difference in a particular coupling does not imply that conclusion. Even a large expected difference in one coupling says nothing about the best coupling's equality probability.

The conclusion is nevertheless true here, for a separate support-counting reason. For the usual undirected 1-out graph on K_n, each vertex independently chooses one of its n-1 neighbors. The result is a spanning tree precisely when the functional graph has one component and a unique directed 2-cycle. Each labelled tree arises from n-1 choices of the doubled edge. Consequently

```
P_1-out[spanning tree] = n^(n-2) / (n-1)^(n-1) ~ e/n.
```

Conditional on being a tree, it is uniform. Thus its TV distance from UST is exactly one minus this probability.

**Fix:** replace the edit-distance inference with this support calculation, or just say the coupling theorem does not give a TV approximation.

## Line-by-line mathematical audit

### Setup, forbidden cell, and electrical identities

The variables are the M undirected edge indicators. Both side bits are functions of their own sides of the split; neither uses information from the other side.

For fixed z, membership in W'_z is fixed. If w belongs to W'_z, every w–W edge is absent. A spanning tree on n >= 24 cannot leave w isolated, so at least one w–U edge is present. Hence `(a_w,b_w) = (1,1)` is impossible. This conclusion would fail for a vertex touched by z; the proof correctly discards those vertices.

F1 is the ordinary deletion/contraction law for feasible UST pinnings. Parallel edges must retain their identities and multiplicities. F2 uses resistance in the graph with the candidate edge deleted when computing its absence probability; that is the correct graph and direction. F3 is used only to upper-bound resistance by deleting helpful edges or undoing identifications.

For F4, the symmetry equations give

```
R_Hd(w,u) = 1/d + (d+1)/(dm),  u not in D.
```

The formula's applications have 1 <= d <= m-1. Identifying some original neighbors in a conditioned minor does not invalidate the comparison to this uncontracted graph.

### Lemma 1: averaged rectangle conditioning

The fibers really are rectangles in the full assignment space:

```
R_z = {x_A : x_A restricted to E(W) cap A equals z cap A}
      x {x_B : x_B restricted to E(W) cap B equals z cap B}.
```

They partition the full cube, not just tree assignments. Some cells have zero UST mass; retain them in the partition and give their contribution value zero. This avoids accidentally assuming that the positive-UST cells carry all of q's mass.

For positive q(R_z), restriction and normalization of a product q is again a product. Therefore

```
BC(mu,q) <= sum_z sqrt(mu(z) q(R_z)) beta(mu_z)
         <= sqrt(sum_z mu(z) beta(mu_z)^2)
            sqrt(sum_z q(R_z)).
```

The last factor is one. Zero-q cells contribute zero without needing to define their normalized q_z. The maximization over q is outside the resulting bound. No unjustified exchange of maxima and averages occurs.

**Pass.** This is an average of squared conditional bounds, not a worst-conditioned-cell bound.

### Lemma 2: mixed vertices for every balanced edge split

The L–H counting uses no vertex-induced assumption about A. Every edge between L and H must be counted either in the A-degree of its L endpoint or the B-degree of its H endpoint. Consequently

```
|L||H| < c(n-1)(|L|+|H|),
```

unless L and H are both empty, in which case all vertices are already mixed. Taking the larger of |L| and |H| at least half their sum gives `min(|L|,|H|) < 2cn`.

If H is the small set, bounding the total A-degree by contributions from L, H and Mid gives

```
n(n-1)/3 <= 2|A| < (n-1)(3cn + |Mid|).
```

The other case uses B. The strict inequalities remain valid in the empty-set edge cases; if both extreme sets are empty, handle that case directly. With c = 1/12, this proves `|Mid| > n/12 >= 2 floor(n/24)`.

An A-clique is not a counterexample. Its interior vertices have both clique edges and edges leaving the clique; for a balanced clique size both degrees are linear in n. The outside, nearly pure-B vertices need not be selected.

**Pass.** The lemma is weaker than possible, but covers arbitrary edge colorings.

### Lemma 3: the difficult conditioning step

Fix a positive-probability pair `(z,sigma)` and a vertex w in W'_z.

1. **The current star is undecided.** The sigma stars belong to other W vertices and contain only edges to U. They share no edges with w–U. All w–W edges are already absent by z. No forced edge touches w, so contraction cannot merge w with another vertex.
2. **The clique is retained.** Neither z nor sigma decides a U–U edge. Some U vertices may be identified by forced paths through other W vertices, but this only shorts the clique. Original star edges can become parallel, never self-loops, because w remains separate.
3. **Removing w leaves a connected graph.** The images of U form a connected clique image. Every other isolated W vertex has at least one forced edge to U in any feasible complete sigma. Every nontrivial z-component has its U-edges undecided. These attach it to U. This justifies the use of the positive tree count tau(G'').

For the leaf bound, fixing one particular labelled B-edge as the sole star edge gives exactly tau(G'') trees. These events are disjoint for different B-edges, including when their other endpoints were identified. Thus the numerator is `d_B tau(G'')`, with no missing multiplicity.

Adding all m star edges back one at a time gives

```
tau(G') / tau(G'') = product_{j=2}^m (1 + R_j)
                  <= product_{j=2}^m [j/(j-1) (1+1/m)]
                  = m (1+1/m)^(m-1) < em.
```

The first added edge is a leaf edge and contributes no factor. For the j-th edge, compare with the uncontracted clique plus its j-1 preceding star edges. The actual graph has additional connections and identifications, so its resistance is **smaller**. The direction of Rayleigh monotonicity is correct. This proves the stated leaf bound pointwise in sigma.

For the deletion bound, after deleting any prefix of A-edges, every B-edge remains. Deleting the current A-edge still leaves w attached, provided d_B >= 1. The graph used in the resistance denominator is connected. Comparing to the clique plus only the B-star gives

```
R <= 1/d_B + (d_B+1)/(d_B m) <= 1/d_B + 2/m.
```

Multiplying the absence probabilities and using `log(1+x) <= x` gives the exact finite-m bound and the stated exponential bound. Contraction and retained extra edges make absence more likely, not less likely.

**Pass.** These are not merely unconditional or average-over-sigma bounds.

The optional general-graph version also follows by distributing the current over the existing star edges and averaging flows between neighbors. Convexity bounds the routing energy by rho, and the star energy is 1/(j-1). It presupposes d > 0 and finite neighbor-to-neighbor resistance to give a useful positive bound. This extension is not needed for Theorem K.

### From detailed sigma to every chain past

This is the central quantifier check. Let H be the earlier pairs of side bits in the chosen order for fixed z. They are a function of sigma. For every past h with positive conditional probability,

```
P[a_w=1 | z,H=h]
  = sum_{sigma: H(sigma)=h} P[sigma | z,H=h] P[a_w=1 | z,sigma]
  >= L_A(w).
```

The same calculation gives the B bound. Crucially, L_A(w) and L_B(w) depend on w, W and the split, not on sigma. Therefore the averaging proves the inequality **for each past h**, rather than only after averaging over h. Zero-probability pasts can be omitted.

Sigma can include future stars without causing a problem: it is only a finer partition used to establish a uniform bound. The coarse past is measurable with respect to it. The forbidden cell also remains absent after this averaging.

The set W'_z and its order are fixed before applying the chain rule to mu_z. The proof does not choose later coordinates using future bits. This meets the supplied adaptive-chain hypothesis.

### Lemma 5: internal-edge tail

For a forest F with k edges, the exact containment probability is

```
P[F subset T] = product_i c_i / n^k <= (2/n)^k,
```

because each component size satisfies `c_i <= 2^(c_i-1)`. Cyclic F has probability zero. The union bound over j-subsets of E(W) gives exactly the claimed tail.

With j = ceil(r/4),

```
P[|z| >= j] <= (e r^2/(jn))^j
             <= (4er/n)^j
             <= (e/6)^(r/4).
```

Here `4er/n <= e/6 < 1`; the last exponent comparison has the right direction. For very small r the event can be impossible, which only helps.

On the complement, integer rounding gives `|z| <= ceil(r/4)-1 < r/4`, so `|W'_z| >= r-2|z| > r/2`.

**Pass.** No independence assumption about tree edges is used.

### Constants and Theorem K

The degree conversion is correct, although it should be read as a ratio calculation, not as independently bounding a numerator and denominator in the same direction:

```
d_A/m >= (cn-r)/(n-r) >= (c-epsilon)/(1-epsilon) = 1/23.
```

The second inequality uses r/n <= epsilon and monotonicity of `(c-t)/(1-t)`. The same holds for B.

After swapping the side names if necessary, theta_A is in [1/23,1/2]. The leaf and deletion bounds give

```
p_A p_B >= theta_A exp(-1-4 theta_A).
```

This function increases to theta_A = 1/4 and then decreases; comparing the two endpoints proves its minimum on the stated interval occurs at 1/23.

Independent evaluation of the actual constants gives

```
p = exp(-1-4/23)/23            = 0.013441514564739847
lambda = (1+sqrt(1-4p))/2      = 0.9863727844310988
-log(lambda)/2                = 0.0068604591582345515
(log(6)-1)/4                  = 0.19793986730701374
1/150                         = 0.006666666666666667.
```

Both exponential comparisons have positive slack. The note's decimal rate is conservative. Averaging with Lemma 1 now gives the claimed `2 exp(-floor(n/24)/150)` bound. For small n it is vacuous but valid.

### Circuit reduction and balance

The paper defines balanced splits by `|A|,|B| >= M/3`, in edge variables. The decomposition lemma produces

```
M/3 < |A_u| <= 2M/3,
|E minus A_u| >= M/3.
```

These are among the splits handled by Theorem K. There is no hidden vertex partition and no need to strengthen Theorem K to another balance constant.

I also checked the decomposition mechanism. Binarizing products costs at most the original nodes-plus-wires size. Above scope 2M/3, smooth sum nodes can be expanded through all children; a binary product has a child above M/3. That child is either already medium-sized or can be expanded further. Decomposability keeps its complementary factor disjoint. The leaf-arity ceiling prevents this descent from stopping at an oversized leaf. The paper also handles a root whose scope omits variables.

The overlap inequality allows non-tree mass and different splits for different terms. It does not require deterministic sums. With F1's hypothesis restored, it gives exactly

```
|C| >= (delta^2/2) exp(floor(n/24)/150).
```

### Proposition U: exactness, determinism, and scope accounting

In a tree on S, delete the edge from v into the component of T-v containing s_0. This uniquely determines S', its attachment vertex u, and the remaining tree containing v. Conversely, two recursive trees and exactly that one cut edge form a tree. This proves exactness and disjoint supports of the sum's children.

The factor scopes are precisely

```
E(S'), E(S minus S'), {vu}, E(S',S minus S') minus {vu}.
```

They are pairwise disjoint and their union is E(S). Every zero constraint on a crossing edge is present. Thus each product is decomposable and each sum is smooth. Repeating the argument recursively shows that each edge variable occurs once in a complete parse, despite DAG sharing between alternative parses.

Singleton constants have empty scope. Normalization gives the uniform law because each accepted assignment has weight one and there are n^(n-2) trees.

Finally, `t(1+t(s-t)) = O(s^3)` and

```
sum_s C(n,s) s^4 2^s = O(n^4 3^n),
```

which proves the stated size upper bound independently of the incorrect finite-size counter. Consequently 2^{Theta(n)} is the right scale for fixed positive overlap, or 2^{Theta(sqrt(M))} in the number of edge variables.

## Independent numerical checks

I wrote fresh enumeration, conditioning and determinant code. I did not import the author's sampler or the common tree enumerator. Jobs ran from `/home/chris/code/eth-research/toy/r3/complete-review/`, with one BLAS thread, no subagents, and at most one test process at a time. Temporary scripts and logs are under `/tmp/`; no numerical files were added to the repository.

### Exact-support affinity checks

Trees were generated independently from every Prüfer sequence, using a heap-based decoder. All balanced splits modulo exchanging A and B were checked at n = 4 and n = 5. At larger n I used fresh random splits and clique, cut and incident-edge splits.

“Exact” here means the full discrete distribution, not exact algebraic arithmetic: its square-root matrices were evaluated with double-precision dense SVD or sparse largest-singular-value computation.

| n | Trees | Splits | Range of beta^2 | Detailed conditioning cells |
|---|---:|---:|---:|---:|
| 4 | 16 | 25 | 0.364277 to 0.562500 | 200 |
| 5 | 125 | 336 | 0.271080 to 0.512000 | 5,376 |
| 6 | 1,296 | 27 | 0.167884 to 0.482253 | 864 |
| 7 | 16,807 | 14 | 0.124683 to 0.260308 | 12,684 |
| 8 | 262,144 | 6 | 0.101656 to 0.237305 | not enumerated |

For all 408 split cases, I checked

```
beta^2 <= B1^2 <= B2^2 <= B3^2 <= 1,
```

including the per-z inequalities. B1 is the rectangle-conditioned exact affinity, B2 is the exact side-bit affinity, and B3 uses the proved finite-m local bounds. There were no failures at tolerance 3e-9. Some B3 bounds equal one; these are genuinely uninformative bounds, not evidence of exponential decay.

The small tests use W of size two or three, since floor(n/24) is zero in this range. They test the lemmas and the B1/B3 mechanism, not the asymptotic rate of Theorem K.

### Exhaustive feasible conditioning at small n

For the tested W and splits at n <= 7, I enumerated every positive-probability cell obtained by fixing z and every subset of the other W stars. This includes nonempty z, complete star revelations, partial revelations, and rare cells. I tested both sides of Lemma 3 and the forbidden-cell assertion in **19,124 conditioning cells**. No inequality failed. Some exact probabilities meet the finite-m lower bound to numerical precision, so there is no large slack masking a systematic error in this test.

### Fresh adversarial sigma constructions at larger n

Instead of taking sigma from a random UST, I constructed feasible forced forests directly, deleted the other revealed star edges, contracted the forest, and computed

```
P[no A-star edge | pins] = tau(minor minus A-star) / tau(minor)
```

using Laplacian log-determinants. Parallel-edge multiplicities were retained.

- 240 cases at n = 24, 48, 72, 96 used random, clique, cut and incident-edge splits, with z empty. Star revelations ranged from leaves to high-degree stars that merged many U vertices. W alternated between the theorem's size and the larger stress-test size floor(n/6).
- 200 further cases at n = 48, 96 used partial star revelations and internal W forests. Of these, 197 had nonempty z and 180 left some other star unrevealed. The target w was always isolated in z.

All 440 cases satisfied the stronger finite-m bounds, hence also the displayed exponential bounds. The minimum ratio of an exact conditional probability to its stronger lower bound was greater than 1.018. These tests deliberately include conditioning that a typical-tree sampler rarely visits.

The first large-n test setup incorrectly demanded mixed degree n/4, stronger than the claimed n/12, and rejected an incident-edge split before checking an inequality. I corrected the setup and reran only that stage. An initial scratch constants calculation was also superseded; the final check uses exactly the p printed above. Neither issue was a counterexample to the claim.

### Upper construction

An independent memoized recurrence gave n^(n-2) for n = 2 through 10. Separately, I enumerated **8,422 summands** through n = 7 and checked that their recursive and literal scopes partition E(S) exactly. This directly checks the crossing-edge issue that counting accepted trees alone would miss. The finite-size counts in F2 were computed independently.

Temporary artifacts:

- `/tmp/r3-independent-review.py`
- `/tmp/r3-extra-review.py`
- `/tmp/r3-review-results.jsonl` (small exact-support checks)
- `/tmp/r3-review-large.jsonl` (large direct-pinning checks; use the corrected constants in the next log)
- `/tmp/r3-review-extra.jsonl` (nonempty z, partial sigma, final constants, scope and size checks)

No test jobs or browsers were left running.

## Independent prior-art check

I ran several deep web searches on November 4, 2026, using combinations of uniform spanning trees, complete graphs, labelled trees, fixed total variation, approximation, decomposable PCs/SPNs and multilinear formulas. I read the relevant primary-source statements rather than treating search snippets as evidence.

1. **Martens and Medabalimi, “On the Expressive Efficiency of Sum Product Networks,” Theorem 37.** [Primary source, section 9.2](https://arxiv.org/html/1411.7717#S9.SS2). The theorem explicitly defines “approximated arbitrarily well” as a sequence of bounded-size D&C SPNs whose outputs converge to the target. Its support-based argument is not a fixed-error TV lower bound allowing non-tree mass. It is important prior work, but not the claimed robust conclusion.
2. **Zhang, Holtzen and Van den Broeck, UAI 2020.** [Primary PDF](https://starai.cs.ucla.edu/papers/ZhangUAI20.pdf). The spanning-tree/DPP separation concerns exact representation and invokes the preceding SPN lower bound. I found no fixed-TV theorem there.
3. **Leland and Choi, “On the Hardness of Approximating Distributions with Tractable Probabilistic Models.”** [Primary source](https://arxiv.org/html/2506.01281v2). Its unconditional approximate-representation separation is against deterministic decomposable PCs; its general tractable-marginal result is computational hardness. Neither supplies the claimed unconditional bound against nondeterministic decomposable PCs for UST(K_n).
4. **Chattopadhyay, Datta, Ghosal and Mukhopadhyay, ITCS 2022.** [Primary PDF](https://drops.dagstuhl.de/storage/00lipics/lipics-vol215-itcs2022/LIPIcs.ITCS.2022.39/LIPIcs.ITCS.2022.39.pdf). Exact monotone spanning-tree lower bounds and the epsilon-sensitive polynomials F +/- epsilon ST are relevant context. Their coefficients/representation and approximation notion are different from fixed-TV approximation of the undirected tree law over edge indicators.
5. **Mészáros and Peled, October 2026.** [Primary source](https://arxiv.org/abs/2610.11554). This is an edit-distance coupling to a random 1-out graph, not a small-TV decomposable circuit representation. See F3 for the distinct TV argument.

I found no prior fixed-positive-overlap or constant-TV lower bound for this exact distribution and this nondeterministic decomposable circuit model. The defensible novelty claim is that specific robust result. The exact exponential lower bound is old, and the subset-DP upper construction should be presented as elementary rather than sold as a separate new lower-bound contribution.

## Bottom line

The proposed escape from the bounded-degree loss works: conditioning does not consume the current star, the retained clique controls resistance after arbitrary feasible pinnings, and the uniform detailed bounds survive every chain past. I could not break those steps analytically or numerically.

Put the a <= 2M/3 hypothesis into Corollary K, repair or remove the finite-size circuit count, and correct the edit-distance explanation. With those changes, the stated exponential lower bound and matching exponential upper bound hold in the intended circuit model.

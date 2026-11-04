# Uniform spanning trees of expanders and complete graphs versus nonnegative probabilistic circuits

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23934515.svg)](https://doi.org/10.5281/zenodo.23934515)

Chris Blaschke (Ax0s). Version 2, November 4, 2026. Preprint, not peer reviewed. [Read the paper (PDF).](paper/expander-spanning-trees.pdf)

The uniform spanning tree of a graph is easy in one sense: by the matrix-tree theorem it has a polynomial-size generating circuit, and every marginal is a ratio of determinants. Those circuits use subtraction. This paper shows that on every 2-connected bounded-degree expander, and on the complete graph, nonnegative probabilistic circuits cannot approximate it at all unless they are exponentially large, even at constant total-variation error and without assuming determinism.

## Results

Let $G$ be a 2-connected graph with $N$ vertices, $M$ edges, maximum degree at most $D$ and edge expansion at least $h>0$, and let $P$ be its uniform spanning tree on the $M$ edge indicators.

- **Main theorem (Theorem 1.1).** A smooth decomposable probabilistic circuit with nonnegative weights and leaves of arity at most $2M/3$ whose distribution is within total variation $\varepsilon<1$ of $P$ has size at least $(1-\varepsilon)^2e^{cN}$, with explicit $c=c(D,h)>0$, polynomial in $h$ for fixed $D$. No girth condition: it covers random $d$-regular graphs (with high probability) and the Lubotzky–Phillips–Sarnak Ramanujan graphs. A spectral constant gives $c\ge8\cdot10^{-8}$ for 6-regular Ramanujan graphs.
- **Separation (Corollary 1.2).** Since $M=\Theta(N)$, the optimal size is $2^{\Theta(M)}$ at every fixed error below one, while $P$ has polynomial-size determinant circuits with subtraction. Exact-support versions hold for circuits and monotone arithmetic circuits.
- **Prisms, squares and the grid (Corollary 1.3, Section 8).** A local cycle-bit argument gives much better constants for prisms and squares of expanders (about $3\cdot10^{-7}$ for a cubic example) and $e^{\Omega(L)}$ for the $L\times L$ grid.
- **Complete and dense graphs (Section 5).** For $K_n$, size at least $(1-\varepsilon)^2e^{n/146}/(e^3\sqrt n)$ at every TV error $\varepsilon<1$, and the optimal size is $2^{\Theta(n)}$. Graphs with spectral gap at least $\kappa n$ get $2^{\Omega(n)}$; two cliques joined by a perfect matching show minimum degree $n/2$ is not enough for the method.
- **Matroids (Section 9).** The cycle-bit and vertex-bit arguments are the circuit and cocircuit cases of one argument for uniform bases of matroids: a local form for every matroid, a global form for balanced matroids, and $2^{\Omega(n)}$ bounds for binary matroids that are neither graphic nor cographic.
- **Does subtraction help? (Section 10).** Structured-decomposable circuits with subtraction need $2^{\Omega(N)}$ on every 2-connected bounded-degree expander, without smoothness, at error $\varepsilon<\tfrac12$, and at every $\varepsilon<1$ when the output sums to one. Normalized squares of structured circuits satisfy the same lower-bound expression as nonnegative circuits. For unstructured circuits: support-exact bounds for bounded product-depth, formulas, and sums of structured circuits.
- **Limits (Section 11).** When balanced splits can work at all (blocks), what star bits can and cannot certify, matroid limits, and spanning forests.

## What this paper does not claim

- The information-theoretic tools are known: the adaptive chain rule is a case of an exercise in Polyanskiy and Wu, and the optimized affinity as a matrix norm is due to Lapidoth and Pfister. The circuit reduction is from [our dimer paper](https://github.com/chrisb4096-alt/square-grid-dimers). Several elementary lemmas may be known; the paper says which.
- Exponential exact lower bounds for the spanning-tree polynomial and the uniform spanning tree of the complete graph were known before (Jerrum–Snir, Martens–Medabalimi, Jukna–Seiwert, Chattopadhyay et al. for directed expanders). The literature search was targeted, not exhaustive; the border-complexity literature was not checked.
- The constants are tiny, so the bounds only bite for very large graphs. Some unstructured-signed constants are astronomically small, and asymptotic only.
- General decomposable circuits with subtraction are open; for exact representation the question is equivalent to the syntactically multilinear complexity of the spanning-tree polynomial. Expanders with $h\to0$ (beyond the grid) and the conjectured $1/\mathrm{poly}(D,1/h)$ constant are open.

## Verify

The checks in `verify/` test the lemmas, chains and constants on small and medium graphs, in exact arithmetic where possible. None of the proofs depends on them. From the repository root, with `g++` and `uv` installed:

```sh
bash verify/run_all.sh
```

It takes about four minutes, stops at the first failure and ends with `ALL CHECKS PASSED`. Saved output is in `verify/out/`. Appendix A of the paper says what each script checks.

To compile the paper with Tectonic:

```sh
cd paper
tectonic expander-spanning-trees.tex
```

## Research notes

`research-notes/` holds the working notes and their adversarial reviews from which version 2 was assembled, kept as a dated record. The paper is the reviewed, final form; where they differ, the paper wins.

## Versions

- v2.0 (November 4, 2026): every 2-connected bounded-degree expander, complete and dense graphs, matroids, signed circuits.
- v1.0 (November 3, 2026): expander prisms and squares.

## Cite

See `CITATION.cff`.

## License

Paper text: CC BY 4.0. Code: MIT (see `LICENSE`).

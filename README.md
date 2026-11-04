# Uniform spanning trees of expander prisms and squares versus nonnegative probabilistic circuits

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23934515.svg)](https://doi.org/10.5281/zenodo.23934515)

Chris Blaschke (Ax0s), November 3, 2026. Preprint, not peer reviewed. [Read the paper (PDF).](paper/expander-spanning-trees.pdf)

The uniform spanning tree of a graph is easy in one sense: by the matrix-tree theorem it has a polynomial-size generating circuit, and every marginal is a ratio of determinants. Those circuits use subtraction. This paper shows that on natural sparse graphs, nonnegative probabilistic circuits cannot approximate it at all unless they are exponentially large, even at constant total-variation error and without assuming determinism.

## Results

Let $H$ be a connected graph with at least three vertices, maximum degree at most $\Delta$ and edge expansion at least $h_0>0$. Let $G$ be the prism $H\square K_2$ or the square $H^2$, with $N$ vertices and $M$ edges, and let $P$ be the uniform spanning tree of $G$ on its $M$ edge indicators.

- **Main theorem (Theorem 1.1).** A smooth decomposable probabilistic circuit with nonnegative weights and leaves of arity at most $2M/3$ whose distribution is within total variation $\varepsilon<1$ of $P$ has size at least $(1-\varepsilon)^2e^{cN}$, with explicit $c=c(\Delta,h_0)>0$.
- **Separation (Corollary 1.2).** Since $M=\Theta(N)$, the optimal size is $2^{\Theta(M)}$ at every fixed error below one, while $P$ has polynomial-size determinant circuits with subtraction.
- **Square grid (Corollary 6.4).** The same argument gives $e^{\Omega(L)}$ for the $L\times L$ grid, tight up to the constant in the exponent.
- **Exact support (Corollary 6.6).** Exponential lower bounds for circuits with exactly the spanning trees as support, and for monotone arithmetic circuits for the spanning-tree polynomial of these graphs.
- **Limits (Section 7).** Pendant edges, large girth, complete graphs, weighted trees, and an extension to uniform spanning forests.

## What this paper does not claim

- The information-theoretic tools are known: the adaptive chain rule is a case of an exercise in Polyanskiy and Wu, and the optimized affinity as a matrix norm is due to Lapidoth and Pfister. The circuit reduction is from [our dimer paper](https://github.com/chrisb4096-alt/square-grid-dimers). What we believe is new is the local cycle-bit reduction for spanning trees and the constant-error bound on sparse bounded-degree families.
- Exponential exact lower bounds for the spanning-tree polynomial and the uniform spanning tree of the complete graph were known before (Jerrum–Snir, Martens–Medabalimi, Jukna–Seiwert, Chattopadhyay et al. for directed expanders). The literature search was targeted, not exhaustive.
- The constants are tiny: $c\approx2\cdot10^{-8}$ for cubic prisms, so the bound is only meaningful for very large graphs. High-girth expanders and unbounded degree are open.

## Verify

The checks in `verify/` test the local lemma, the chain rule, the geometry and the constants on small graphs, in exact arithmetic where possible. None of the proofs depends on them. From the repository root, with `g++` and `uv` installed:

```sh
bash verify/run_all.sh
```

It stops at the first failure and ends with `ALL CHECKS PASSED`. Saved output is in `verify/out/`. Appendix A of the paper says what each script checks.

To compile the paper with Tectonic:

```sh
cd paper
tectonic expander-spanning-trees.tex
```

## Cite

See `CITATION.cff`.

## License

Paper text: CC BY 4.0. Code: MIT (see `LICENSE`).

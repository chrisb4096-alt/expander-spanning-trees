# R3 research notes (November 4, 2026)

These are working notes for follow-up results to the paper in this repository, posted here to timestamp them. They are not peer-reviewed. Each note was checked by one adversarial review, saved next to it as `*-review.md`, and that review's fixes were applied. The v2 paper will contain the polished versions.

| Note | Result | Review |
|---|---|---|
| [R3-complete.md](R3-complete.md) | UST(K_n): every balanced edge split has β² ≤ 2e^{−⌊n/24⌋/150}, so nonnegative smooth decomposable PCs (leaf arity ≤ 2M/3) within TV 1−δ need size ≥ (δ²/2)e^{⌊n/24⌋/150}. An exact O(n⁴3ⁿ) circuit makes 2^{Θ(n)} tight. | [holds with fixes](R3-complete-review.md) |
| [R3-matroid.md](R3-matroid.md) | Uniform matroid bases: the argument is circuit/cocircuit duality. Framed-piece bound for every matroid; leverage version for balanced matroids, which include regular matroids. | [holds with fixes](R3-matroid-review.md) |
| [R3-signed.md](R3-signed.md) | Structured (fixed-vtree) signed circuits need exponential size, exactly or within TV < 1/2; so do squares and PSD forms of structured-decomposable circuits. The unstructured signed case is equivalent to syntactically multilinear circuits for the spanning-tree polynomial (open). | [holds with fixes](R3-signed-review.md) |
| [R3-constants.md](R3-constants.md) | Better constants for 2-connected bounded-degree expanders: ln(1/c) = O(D log(D/h)) for h ≤ 1, and near-polynomial in D when λ₂ > 1. | [holds with fixes](R3-constants-review.md) |

File paths in the notes refer to the private working repository. The numerics scripts are not included here.

# /// script
# dependencies = ["numpy"]
# ///
"""Random stress test of the chain inequality (Lemma 3.3) combined with the 2x2 bound (Lemma 3.5).

mu is a law on ({0,1}^2)^r built by the chain rule with arbitrary conditionals kappa_k(.|past) that may
depend on the whole joint past, put zero mass on (1,1), and put at least p on (0,1) and on (1,0).
Each conditional is drawn from the polytope {a+b+c=1, a>=0, b>=p, c>=p}, whose vertices are
(1-2p,p,p), (0,1-p,p) and (0,p,1-p). The extreme modes include zero mass on (0,0) and p up to 1/2.

Checks, for every chain:
  ||sqrt mu||_op <= prod_k max_past ||sqrt kappa_k(.|past)||_op <= ((1+sqrt(1-4p^2))/2)^(r/2) <= (1-p^2)^(r/2).
Singular values are floating point; the tolerance is 1e-12.
Exits nonzero (AssertionError) on any violation.

Copied from toy/r1/proof/check_tensorization.py and extended with the boundary kernels
(zero (0,0) mass, p in {0.45, 0.5}) requested in the referee review.
"""
import sys
import numpy as np

rng = np.random.default_rng(7)
TOL = 1e-12
worst_slack = np.inf
counts = {}
TRIALS = 6000
for trial in range(TRIALS):
    r = int(rng.integers(1, 7))
    p = float(rng.choice([0.0, 0.01, 0.05, 0.2, 1 / 3, 0.45, 0.5]))
    mode = int(rng.integers(4))
    counts[mode] = counts.get(mode, 0) + 1
    verts = np.array([[1 - 2 * p, p, p], [0.0, 1 - p, p], [0.0, p, 1 - p]])
    law = {((), ()): 1.0}
    gam = []
    for k in range(r):
        new = {}
        g = 0.0
        for (u, v), w in law.items():
            if mode == 0:  # random interior point of the polytope
                x = rng.dirichlet([0.3, 0.3, 0.3]) @ verts
            elif mode == 1:  # extreme, past-dependent: alternate between vertices by parity of the past
                x = verts[(sum(u) + 2 * sum(v)) % 3]
            elif mode == 2:  # near-degenerate, copy-like dependence on the past
                x = rng.dirichlet([0.05, 0.05, 0.05]) @ verts
            else:  # boundary: zero mass on (0,0), split between the two cross cells
                s = float(rng.uniform(0, 1))
                x = s * verts[1] + (1 - s) * verts[2]
            a, b, c = x  # (0,0), (0,1), (1,0); (1,1) forbidden
            assert b >= p - 1e-15 and c >= p - 1e-15 and a >= -1e-15
            a = max(a, 0.0)
            M = np.sqrt(np.array([[a, b], [c, 0.0]]))
            g = max(g, np.linalg.norm(M, 2))
            for (uu, vv, q) in [(0, 0, a), (0, 1, b), (1, 0, c)]:
                if q > 0:
                    key = (u + (uu,), v + (vv,))
                    new[key] = new.get(key, 0) + w * q
        law = new
        gam.append(g)
    idx = lambda t: int("".join(map(str, t)) or "0", 2)
    S = np.zeros((1 << r, 1 << r))
    for (u, v), w in law.items():
        S[idx(u), idx(v)] += w
    op = np.linalg.norm(np.sqrt(S), 2)
    bound = float(np.prod(gam))
    sharp = ((1 + np.sqrt(max(0.0, 1 - 4 * p * p))) / 2) ** (r / 2)
    weak = (1 - p * p) ** (r / 2)
    assert op <= bound + TOL, (trial, op, bound)
    assert bound <= sharp + TOL, (trial, bound, sharp)
    assert sharp <= weak + TOL, (trial, sharp, weak)
    worst_slack = min(worst_slack, bound - op)
print(f"{TRIALS} random chains (modes {dict(sorted(counts.items()))}): "
      f"||sqrt mu|| <= prod gamma_k <= gamma(p)^r <= (1-p^2)^(r/2) held; min slack {worst_slack:.3e}")
sys.exit(0)

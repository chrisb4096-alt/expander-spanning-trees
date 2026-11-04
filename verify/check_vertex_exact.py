# /// script
# dependencies = ["numpy", "scipy", "networkx"]
# ///
"""Exact checks of the vertex-touch argument on small 2-connected graphs (ported from toy/r2/exact_small.py).

For P = UST(G), a balanced edge split A|B (|A|, |B| >= M/3) and an independent set W of mixed
vertices, all spanning trees are enumerated and the following are checked exactly (floating point,
tolerance 1e-9):

  beta = ||sqrt P||_op <= ||sqrt mu||_op <= prod_k gamma_k <= prod_k sqrt((1+sqrt(1-4 pA pB))/2)

(Lemma "data processing", Lemma "adaptive chain rule", Theorem "vertex-touch bound for a split"),
where mu is the joint law of the bits a_w = 1[T has no A-edge at w], b_w = 1[T has no B-edge at w]
and pA, pB are the resistance bounds of Lemma "touch probabilities". The lemma itself is checked
for every configuration sigma of T on the stars of W minus {w} (its exact hypothesis), and also for
the coarser conditioning on the stars of the earlier vertices only, which is what the chain rule uses.
Finally pA*pB >= (2+rho_w)^(-deg w) is checked. When G - W is connected, the telescoped bound of
Lemma "telescoped touch bound" is also checked for every sigma: Pr[a_w=1 | sigma] >= (b/d)(1+kappa)^(-a)
and Pr[b_w=1 | sigma] >= (a/d)(1+kappa)^(-b), with a, b the numbers of A- and B-edges at w, d = a+b
and kappa = 1/lambda_2(L_{G-W}); and the chain is checked against the product of these bounds.
Exits with status 1 if any check fails.

Usage: uv run check_vertex_exact.py --graph petersen --cuts 2 --search 1 --steps 100 --seed 1
"""
import argparse, json, sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import svds
from common import *

TOL = 1e-9


def op_sqrtP(T, maskA, tau):
    a = T & np.uint64(maskA)
    b = T & ~np.uint64(maskA)
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    M = sp.csr_matrix((np.full(len(T), 1.0 / np.sqrt(tau)), (ia, ib)), shape=(len(ua), len(ub)))
    if min(M.shape) <= 2:
        return float(np.linalg.norm(M.toarray(), 2))
    return float(svds(M, k=1, return_singular_vectors=False)[0])


def balanced(mask, m):
    k = bin(mask).count("1")
    return 3 * k >= m and 3 * (m - k) >= m


def random_balanced_mask(m, rng):
    while True:
        k = int(rng.integers(m // 3, 2 * m // 3 + 2))
        idx = rng.choice(m, size=min(k, m), replace=False)
        mask = 0
        for i in idx:
            mask |= 1 << int(i)
        if balanced(mask, m):
            return mask


def local_search(T, m, tau, rng, steps):
    """Hill-climb on beta over balanced splits (looks for the hardest splits)."""
    mask = random_balanced_mask(m, rng)
    cur = op_sqrtP(T, mask, tau)
    for _ in range(steps):
        new = mask ^ (1 << int(rng.integers(m)))
        if not balanced(new, m):
            continue
        val = op_sqrtP(T, new, tau)
        if val >= cur:
            mask, cur = new, val
    return mask, cur


def bisection_mask(G, n, E, rng):
    """Fiedler vertex bisection; A = edges inside side 0 plus cut edges until balanced (few mixed vertices)."""
    L = nx.laplacian_matrix(G, nodelist=range(n)).toarray().astype(float)
    w, V = np.linalg.eigh(L)
    f = V[:, 1] + 1e-9 * rng.standard_normal(n)
    side = f < np.median(f)
    mask = 0
    cut = []
    for i, (u, v) in enumerate(E):
        if side[u] and side[v]:
            mask |= 1 << i
        elif side[u] != side[v]:
            cut.append(i)
    rng.shuffle(cut)
    m = len(E)
    for i in cut:
        if balanced(mask, m):
            break
        mask |= 1 << i
    return mask if balanced(mask, m) else None


def run_cut(G, n, E, T, maskA, rng, sep):
    m = len(E); tau = len(T)
    inA = [bool((maskA >> i) & 1) for i in range(m)]
    S = star(n, E)
    Z = mixed_vertices(n, E, inA)
    if not Z:
        return None
    # W: greedy separated subset of Z. Prefer the largest W with G - W connected (larger bounds);
    # the bounds are valid without that, so fall back to the largest W found.
    best, best_conn = None, None
    for t in range(40):
        W = greedy_separated(G, Z, sep, rng)
        keep = [x for x in range(n) if x not in set(W)]
        if best is None or len(W) > len(best):
            best = W
        if nx.is_connected(G.subgraph(keep)) and (best_conn is None or len(W) > len(best_conn)):
            best_conn = W
    W = best_conn if best_conn is not None else best
    beta = op_sqrtP(T, maskA, tau)
    rb = touch_bounds(n, E, W, inA)
    SAm = [sum(1 << i for i in S[w] if inA[i]) for w in W]
    SBm = [sum(1 << i for i in S[w] if not inA[i]) for w in W]
    starm = [sum(1 << i for i in S[w]) for w in W]
    r = len(W)
    U = np.zeros(len(T), dtype=np.int64); Vb = np.zeros(len(T), dtype=np.int64)
    past_mask = 0
    gam, key_past, key_full, theo, coarse_ok, flow = [], [], [], [], [], []
    keepW = [x for x in range(n) if x not in set(W)]
    connW = nx.is_connected(G.subgraph(keepW))
    kappa = 1.0 / algebraic_connectivity(laplacian(n, E, keep_vertices=set(keepW)), keepW) if connW else None
    for k in range(r):
        a = (T & np.uint64(SAm[k])) == 0
        b = (T & np.uint64(SBm[k])) == 0
        if np.any(a & b):
            raise AssertionError("a spanning tree misses every edge at a vertex")
        key = U * (1 << r) + Vb
        loc = a.astype(np.int64) * 2 + b.astype(np.int64)
        keys, inv = np.unique(key, return_inverse=True)
        C = np.zeros((len(keys), 4)); np.add.at(C, (inv, loc), 1.0)
        C /= C.sum(axis=1, keepdims=True)
        g = 0.0
        for row in C:
            g = max(g, np.linalg.norm(np.sqrt(np.array([[row[0], row[1]], [row[2], row[3]]])), 2))
        gam.append(g)
        mins = []
        for cond in (past_mask, sum(starm[j] for j in range(r) if j != k)):
            sig = T & np.uint64(cond)
            sk, sinv = np.unique(sig, return_inverse=True)
            tot = np.bincount(sinv).astype(float)
            hA = np.bincount(sinv, weights=a.astype(float), minlength=len(sk))
            hB = np.bincount(sinv, weights=b.astype(float), minlength=len(sk))
            mins.append((float(np.min(hA / tot)), float(np.min(hB / tot))))
        key_past.append(mins[0]); key_full.append(mins[1])
        pA, pB = rb[k]["A"], rb[k]["B"]
        theo.append((pA, pB))
        coarse_ok.append(pA * pB >= rb[k]["A_coarse"] * rb[k]["B_coarse"] * (1 - TOL))
        if connW:
            na, nb = bin(SAm[k]).count("1"), bin(SBm[k]).count("1"); d = na + nb
            flow.append((nb / d * (1 + kappa) ** (-na), na / d * (1 + kappa) ** (-nb)))
        U = U | (a.astype(np.int64) << k); Vb = Vb | (b.astype(np.int64) << k)
        past_mask |= starm[k]
    Mu = np.zeros((1 << r, 1 << r)); np.add.at(Mu, (U, Vb), 1.0); Mu /= tau
    op_mu = float(np.linalg.norm(np.sqrt(Mu), 2))
    chain = float(np.prod(gam))
    bound = float(np.prod([np.sqrt((1 + np.sqrt(max(0.0, 1 - 4 * pA * pB))) / 2) for pA, pB in theo]))
    def key_holds(km):
        return all(x[0] >= t[0] - TOL and x[1] >= t[1] - TOL for x, t in zip(km, theo))
    ok = bool(beta <= op_mu + TOL and op_mu <= chain + TOL and chain <= bound + TOL)
    ratios = [min(x[0] / t[0], x[1] / t[1]) for x, t in zip(key_full, theo) if t[0] > 0 and t[1] > 0]
    if connW:
        flow_bound = float(np.prod([np.sqrt((1 + np.sqrt(max(0.0, 1 - 4 * x * y))) / 2) for x, y in flow]))
        flow_ok = bool(all(x[0] >= t[0] - TOL and x[1] >= t[1] - TOL for x, t in zip(key_full, flow))
                       and chain <= flow_bound + TOL)
        flow_ratio = round(min(min(x[0] / t[0], x[1] / t[1]) for x, t in zip(key_full, flow)), 4)
    else:
        flow_bound, flow_ok, flow_ratio = None, True, None
    return dict(Z=len(Z), W=len(W), G_minus_W_connected=best_conn is not None,
                beta=round(beta, 6), op_mu=round(op_mu, 6), chain=round(chain, 6), bound=round(bound, 6),
                chain_ok=ok, key_full_ok=key_holds(key_full), key_past_ok=key_holds(key_past),
                coarse_ok=bool(all(coarse_ok)), flow_ok=flow_ok, flow_bound=None if flow_bound is None else round(flow_bound, 6),
                flow_min_ratio=flow_ratio,
                key_min_ratio=round(min(ratios), 4) if ratios else None,
                pA=[round(t[0], 4) for t in theo], pB=[round(t[1], 4) for t in theo])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default="petersen")
    ap.add_argument("--cuts", type=int, default=2)
    ap.add_argument("--search", type=int, default=1)
    ap.add_argument("--steps", type=int, default=100)
    ap.add_argument("--sep", type=int, default=2, help="pairwise distance in W (2 = independent)")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    G = named_graph(args.graph, seed=args.seed)
    n = G.number_of_nodes(); E = edge_list(G); m = len(E)
    T = enum_trees(n, E)
    tau = len(T)
    assert tau == round(np.exp(np.linalg.slogdet(laplacian(n, E)[1:, 1:])[1])), "tree count != matrix-tree theorem"
    print(f"# {args.graph}: N={n} M={m} tau={tau} girth={nx.girth(G)} kappa={nx.node_connectivity(G)}")
    cuts = [("random", random_balanced_mask(m, rng)) for _ in range(args.cuts)]
    bm = bisection_mask(G, n, E, rng)
    if bm is not None:
        cuts.append(("bisection", bm))
    for _ in range(args.search):
        mask, _v = local_search(T, m, tau, rng, args.steps)
        cuts.append(("search", mask))
    failed = 0
    for kind, mask in cuts:
        assert balanced(mask, m)
        res = run_cut(G, n, E, T, mask, rng, args.sep)
        if res is None:
            print(json.dumps(dict(graph=args.graph, kind=kind, note="no mixed vertex")))
            failed += 1          # impossible for a balanced split of a connected graph
            continue
        good = res["chain_ok"] and res["key_full_ok"] and res["key_past_ok"] and res["coarse_ok"] and res["flow_ok"]
        failed += not good
        res = dict(graph=args.graph, kind=kind, N=n, M=m, status="PASS" if good else "FAIL", **res)
        print(json.dumps(res))
    print(f"# {args.graph}: {len(cuts) - failed}/{len(cuts)} splits pass")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

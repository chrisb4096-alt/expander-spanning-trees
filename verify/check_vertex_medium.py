# /// script
# dependencies = ["numpy", "scipy", "networkx"]
# ///
"""Medium-size checks of the vertex-touch argument (ported from toy/r2/key_lemma.py).

For a 2-connected graph G (cages up to girth 12, and 2-connected random cubic graphs) and three
balanced splits per graph (random half, Fiedler bisection, BFS ball):
  * |Z| = number of mixed vertices, W = greedy subset of Z with pairwise distance >= sep,
  * lambda_2(L_{G-W}) and rho = max_w max_{x,y in N(w)} R_{G-W}(x,y),
  * the bounds pA(w), pB(w) of Lemma "touch probabilities" and the coarse bound (2+rho_w)^(-deg w),
  * the exponent sum_w -1/2 log((1+sqrt(1-4 pA pB))/2) of Theorem "vertex-touch bound for a split",
    so that beta <= exp(-exponent) for that split (floating point; a bound for this split only),
  * spot checks of Lemma "touch probabilities": sample T ~ UST(G) by Wilson's algorithm, pick
    w in W, condition on T restricted to the stars of W minus {w}, and compute
    Pr[T has no A-edge at w | sigma] exactly as a ratio of spanning-tree counts (determinants of
    contracted/deleted Laplacians). It must be >= pA(w), and likewise for B.
Exits with status 1 on any violation, or if a split has no mixed vertex.

The rr<d>_<n> graphs are resampled until 2-connected, so they are not evidence for statements
about random regular graphs in general.

Usage: uv run check_vertex_medium.py --graph rr3_200 --sep 3 --samples 20 --seed 2
"""
import argparse, json, sys
import numpy as np
from common import *


def wilson(n, adj, rng, root=0):
    in_tree = np.zeros(n, bool); in_tree[root] = True
    nxt = -np.ones(n, int); nxt_e = -np.ones(n, int)
    for s in rng.permutation(n):
        u = s
        while not in_tree[u]:
            j = rng.integers(len(adj[u]))
            v, e = adj[u][j]
            nxt[u] = v; nxt_e[u] = e
            u = v
        u = s
        while not in_tree[u]:
            in_tree[u] = True
            u = nxt[u]
    return set(int(nxt_e[u]) for u in range(n) if u != root)


def logtau_minor(n, E, contract, delete, extra_delete=()):
    """log tau of (G - delete - extra_delete) / contract, keeping parallel edges, dropping loops."""
    par = list(range(n))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i in contract:
        a, b = f(E[i][0]), f(E[i][1])
        assert a != b, "contracted set is not a forest"
        par[a] = b
    roots = sorted({f(x) for x in range(n)})
    idx = {r: k for k, r in enumerate(roots)}
    k = len(roots)
    L = np.zeros((k, k))
    dl = set(delete) | set(extra_delete) | set(contract)
    for i, (u, v) in enumerate(E):
        if i in dl:
            continue
        a, b = idx[f(u)], idx[f(v)]
        if a == b:
            continue
        L[a, a] += 1; L[b, b] += 1; L[a, b] -= 1; L[b, a] -= 1
    s, ld = np.linalg.slogdet(L[1:, 1:])
    return ld if s > 0 else -np.inf


def make_split(G, n, E, kind, rng):
    m = len(E)
    if kind == "random":
        idx = set(rng.choice(m, size=m // 2, replace=False).tolist())
        return [i in idx for i in range(m)]
    if kind == "bisection":   # Fiedler vertex bisection: few mixed vertices
        L = nx.laplacian_matrix(G, nodelist=range(n)).toarray().astype(float)
        w, V = np.linalg.eigh(L)
        side = V[:, 1] < np.median(V[:, 1])
        inA = [bool(side[u] and side[v]) for u, v in E]
        cut = [i for i, (u, v) in enumerate(E) if side[u] != side[v]]
        rng.shuffle(cut)
        for i in cut:
            if 3 * sum(inA) >= m:
                break
            inA[i] = True
        return inA
    if kind == "bfs":   # the m//2 edges closest to a random vertex
        src = int(rng.integers(n))
        dist = nx.single_source_shortest_path_length(G, src)
        es = sorted(range(m), key=lambda i: (min(dist[E[i][0]], dist[E[i][1]]), max(dist[E[i][0]], dist[E[i][1]])))
        A = set(es[: m // 2])
        return [i in A for i in range(m)]
    raise ValueError(kind)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default="rr3_200")
    ap.add_argument("--splits", default="random,bisection,bfs")
    ap.add_argument("--sep", type=int, default=3)
    ap.add_argument("--samples", type=int, default=20)
    ap.add_argument("--seed", type=int, default=2)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    G = named_graph(args.graph, seed=args.seed)
    n = G.number_of_nodes(); E = edge_list(G); m = len(E)
    S = star(n, E)
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(E):
        adj[u].append((v, i)); adj[v].append((u, i))
    lam_G = algebraic_connectivity(laplacian(n, E), range(n))
    print(f"# {args.graph}: N={n} M={m} girth={nx.girth(G)} lambda2(L_G)={lam_G:.4f} D={max(d for _, d in G.degree())}")
    failed = 0
    for kind in args.splits.split(","):
        inA = make_split(G, n, E, kind, rng)
        assert 3 * sum(inA) >= m and 3 * (m - sum(inA)) >= m, "split not balanced"
        Z = mixed_vertices(n, E, inA)
        if not Z:
            print(json.dumps(dict(graph=args.graph, split=kind, note="no mixed vertex", status="FAIL")))
            failed += 1
            continue
        best, best_conn = None, None
        for t in range(10):
            W = greedy_separated(G, Z, args.sep, rng)
            keep = [x for x in range(n) if x not in set(W)]
            if best is None or len(W) > len(best):
                best = W
            if nx.is_connected(G.subgraph(keep)) and (best_conn is None or len(W) > len(best_conn)):
                best_conn = W
        W = best_conn if best_conn is not None else best
        keep = [x for x in range(n) if x not in set(W)]
        lamH = algebraic_connectivity(laplacian(n, E, keep_vertices=set(keep)), keep)
        rb = touch_bounds(n, E, W, inA)
        rho = max(x["rho"] for x in rb)
        site = [-0.5 * np.log((1 + np.sqrt(max(0.0, 1 - 4 * x["A"] * x["B"]))) / 2) for x in rb]
        site_coarse = [-0.5 * np.log((1 + np.sqrt(max(0.0, 1 - 4 * x["A_coarse"] * x["B_coarse"]))) / 2) for x in rb]
        coarse_viol = sum(x["A"] * x["B"] < x["A_coarse"] * x["B_coarse"] * (1 - 1e-9) for x in rb)
        worst_ratio = np.inf; viol = 0; checked = 0
        Wpos = {w: k for k, w in enumerate(W)}
        for s in range(args.samples):
            T = wilson(n, adj, rng)
            assert len(T) == n - 1
            w = W[int(rng.integers(len(W)))]
            past = set(i for x in W if x != w for i in S[x])
            contract = [i for i in past if i in T]
            delete = [i for i in past if i not in T]
            lt = logtau_minor(n, E, contract, delete)
            assert np.isfinite(lt), "observed configuration has probability 0"
            for side in ("A", "B"):
                rem = [i for i in S[w] if inA[i] == (side == "A")]
                pr = float(np.exp(logtau_minor(n, E, contract, delete, extra_delete=rem) - lt))
                bound = rb[Wpos[w]][side]
                checked += 1
                if bound > 0:
                    worst_ratio = min(worst_ratio, pr / bound)
                viol += pr < bound * (1 - 1e-9)
        good = viol == 0 and coarse_viol == 0
        failed += not good
        out = dict(graph=args.graph, split=kind, status="PASS" if good else "FAIL", N=n, sep=args.sep,
                   Z=len(Z), Z_frac=round(len(Z) / n, 4), W=len(W), W_frac=round(len(W) / n, 4),
                   G_minus_W_connected=best_conn is not None,
                   lambda2_H=round(float(lamH), 5), rho=round(float(rho), 3),
                   pA_min=round(float(min(x["A"] for x in rb)), 6), pB_min=round(float(min(x["B"] for x in rb)), 6),
                   exponent_per_N=float(f"{sum(site) / n:.4g}"), coarse_exponent_per_N=float(f"{sum(site_coarse) / n:.4g}"),
                   key_checks=checked, key_violations=int(viol), coarse_violations=int(coarse_viol),
                   key_min_ratio=round(float(worst_ratio), 4))
        print(json.dumps(out), flush=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()

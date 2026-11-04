# /// script
# dependencies = ["numpy", "networkx", "scipy"]
# ///
"""Checks of the path-packing theorem for sums of structured circuits (Section "Does subtraction help?",
Lemmas lem:connector, lem:blob, lem:bond, lem:sep-delete, Theorem thm:packing and Theorem thm:signed-sums-linear).
Ported from toy/r5/approx-path-packing/check_pp.py.

[A] Lemma lem:bond (bond through two edges): on random small 2-connected (multi)graphs, for every
    pair of distinct edges e, f, the two-disjoint-paths construction returns a vertex set Y with
    G[Y] and G - Y connected and e, f both crossing.  Cross-checked against brute force over all
    vertex subsets for N <= 9.
[B] The bond gadget in context: for short connector paths U in expanders and pocketed graphs,
    B(U) = V - giant(G - U), |B(U)| <= |U| (1 + D/h_est), and the bond gadget Y is a valid
    co-connected gadget with e, f on its boundary.
[C] The packing construction of Theorem thm:packing (short connectors -> blob -> bond gadget, greedy with
    blocked balls) on random cubic, random 4-regular, pocketed cubic, truncated cubic and prism graphs, for
    far P, Q (delta = 1/3, 1/9), vtree-refined P, Q (r = 2, 3) and adversarial pocket P, Q, with gadgets
    pairwise at distance >= d for d = 2, 4, 6 (the proof uses d = sigma + 1, far larger).  Reports
    gadgets/N, sizes, absorption rate, connectivity and lambda_2 of G - W.  Numerical only; sizes from
    the environment variable PP_SIZES.
[D] Exact checks on small graphs (spanning trees enumerated): for the constructed gadget family
    (including absorbed, non-path blobs), the pendant restriction gives prod XOR (rank 2^m),
    pair blocks are cubes with disjoint permutation rectangles, and the conditional pendant bound
    of Lemma lem:pendant-prob holds with kappa' = 2 kappa (and kappa' = kappa when G/U is simple).
[E] Explicit constants of Theorem thm:packing and Theorem thm:signed-sums-linear (floating point, log10).
[F] Lemmas lem:connector and lem:blob, and Lemma lem:sep-delete with one gadget, against exact edge
    expansion h on small graphs; imported Cheeger-Mohar bound lambda_2 >= mu(D, h).

Exits nonzero if any check fails.
"""
import itertools
import math
import os
import random
import sys
import time
from collections import defaultdict, deque

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "4")

import networkx as nx
import numpy as np
from common import enum_trees

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []
GAM = ("gamma",)


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL:", msg, flush=True)
    return cond


# ===================================================================== utilities
def lam2(G):
    if G.number_of_nodes() < 2 or not nx.is_connected(G):
        return 0.0
    if G.number_of_nodes() <= 400:
        L = nx.laplacian_matrix(G).toarray().astype(float)
        return float(np.linalg.eigvalsh(L)[1])
    return float(nx.algebraic_connectivity(G, method="tracemin_lu", seed=1))


def ekey(u, v):
    return (u, v) if u <= v else (v, u)


def connected_regular(d, n, seed):
    while True:
        H = nx.random_regular_graph(d, n, seed=seed)
        if nx.is_connected(H) and nx.is_biconnected(H):
            return H
        seed += 7919


def pocket_graph(H):
    """Each edge uv of H becomes the path u-a-b-v plus a degree-2 vertex c adjacent to a and b.
    Max degree max(3, Delta(H)); 2-connected if H is; {a, b} cuts off the pocket {c}."""
    G = nx.Graph()
    G.add_nodes_from(H.nodes())
    nxt = max(H.nodes()) + 1
    units = []
    for u, v in H.edges():
        a, b, c = nxt, nxt + 1, nxt + 2
        nxt += 3
        G.add_edges_from([(u, a), (a, b), (b, v), (a, c), (c, b)])
        units.append((u, a, b, v, c))
    return G, units


def truncate(H):
    """Replace every vertex by a triangle (cubic H -> cubic graph with many triangles)."""
    G = nx.Graph()
    port = {}
    nxt = 0
    for v in H.nodes():
        nb = list(H.neighbors(v))
        ids = list(range(nxt, nxt + len(nb)))
        nxt += len(nb)
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                G.add_edge(ids[i], ids[j])
        for w, x in zip(nb, ids):
            port[(v, w)] = x
    for u, v in H.edges():
        G.add_edge(port[(u, v)], port[(v, u)])
    return G


def prism(H):
    n = H.number_of_nodes()
    G = nx.Graph()
    for u, v in H.edges():
        G.add_edge(u, v)
        G.add_edge(u + n, v + n)
    for x in range(n):
        G.add_edge(x, x + n)
    return G


# --- P, Q sources (copied from toy/r5/approx-restrictions/check_r5a.py and toy/r4/.../r4s_lib.py)
def geometric_vtree(G, E, root, rng):
    dist = nx.single_source_shortest_path_length(G, root)
    order = sorted(range(len(E)), key=lambda i: (min(dist[E[i][0]], dist[E[i][1]]), rng.random()))

    def build(lo, hi):
        if hi - lo == 1:
            return order[lo]
        mid = (lo + hi) // 2
        return (build(lo, mid), build(mid, hi))
    return build(0, len(order))


def vtree_nodes(T):
    out = []

    def rec(x):
        s = rec(x[0]) | rec(x[1]) if isinstance(x, tuple) else frozenset([x])
        out.append(s)
        return s
    rec(T)
    return out


def balanced_node(scopes, U):
    U = set(U)
    for S in scopes:
        if len(U) / 3 <= len(S & U) <= 2 * len(U) / 3:
            return S
    raise AssertionError("no balanced node")


def refine(vtrees_scopes, E0):
    S = balanced_node(vtrees_scopes[0], E0)
    P, Q = sorted(set(E0) & S), sorted(set(E0) - S)
    n = min(len(P), len(Q))
    P, Q = P[:n], Q[:n]
    for scopes in vtrees_scopes[1:]:
        S = balanced_node(scopes, set(P) | set(Q))
        a = [x for x in P if x in S]; b = [x for x in P if x not in S]
        c = [x for x in Q if x in S]; d = [x for x in Q if x not in S]
        P, Q = (a, d) if min(len(a), len(d)) >= min(len(b), len(c)) else (b, c)
        n = min(len(P), len(Q))
        P, Q = P[:n], Q[:n]
    return P, Q


def far_sets(G, E, delta, rng):
    root = rng.choice(list(G.nodes()))
    dist = nx.single_source_shortest_path_length(G, root)
    key = sorted(range(len(E)), key=lambda i: (min(dist[E[i][0]], dist[E[i][1]]), rng.random()))
    k = int(delta * len(E))
    return [E[i] for i in key[:k]], [E[i] for i in key[-k:]]


# ===================================================================== Lemma CB
def two_paths_bond(GB, e, f):
    """GB simple connected graph (2-connected), e=(a,a2), f=(b,b2) its edges (as endpoint pairs,
    possibly equal pairs = parallel edges).  Returns (Y0, Z0): disjoint connected vertex sets with
    a in Y0, a2 in Z0 and f's endpoints split between them (Menger: two disjoint paths)."""
    a, a2 = e
    b, b2 = f
    if {a, a2} == {b, b2}:
        return {a}, {a2}
    common = {a, a2} & {b, b2}
    if common:
        x = common.pop()
        y = ({a, a2} - {x}).pop()
        z = ({b, b2} - {x}).pop()
        path = nx.shortest_path(GB.subgraph(set(GB) - {x}), y, z)
        return ({x}, set(path)) if x == a else (set(path), {x})
    aux = GB.copy()
    aux.add_edges_from([("s*", a), ("s*", a2), ("t*", b), ("t*", b2)])
    paths = list(nx.node_disjoint_paths(aux, "s*", "t*"))
    assert len(paths) == 2, "2-connectivity violated"
    pa = next(p for p in paths if p[1] == a)
    pa2 = next(p for p in paths if p[1] == a2)
    return set(pa[1:-1]), set(pa2[1:-1])


def extend_partition(GB, Y0, Z0, prefer=None):
    """Extend disjoint connected Y0, Z0 to a partition of V(GB) into two connected sides.
    If prefer (a vertex) is given, it is attached first and its side is grown maximally."""
    sides = [set(Y0), set(Z0)]
    assigned = lambda: sides[0] | sides[1]
    if prefer is not None and prefer in GB and prefer not in assigned():
        par = {prefer: None}
        dq = deque([prefer])
        hit = None
        while dq and hit is None:
            x = dq.popleft()
            for y in GB.neighbors(x):
                if y in assigned():
                    hit = (x, y)
                    break
                if y not in par:
                    par[y] = x
                    dq.append(y)
        x, y = hit
        s = 0 if y in sides[0] else 1
        while x is not None:
            sides[s].add(x)
            x = par[x]
    order = [0, 1]
    if prefer is not None and prefer in sides[1]:
        order = [1, 0]
    for s in order:
        dq = deque(sides[s])
        while dq:
            x = dq.popleft()
            for y in GB.neighbors(x):
                if y not in sides[0] and y not in sides[1]:
                    sides[s].add(y)
                    dq.append(y)
    assert len(sides[0]) + len(sides[1]) == GB.number_of_nodes()
    return sides


def giant_complement(G, U):
    U = set(U)
    rest = G.subgraph([v for v in G if v not in U])
    comps = list(nx.connected_components(rest))
    giant = max(comps, key=len) if comps else set()
    return set(G) - giant, giant


def bond_gadget(G, U, e, f):
    """Lemma CB: U connected, e, f distinct edges each with an endpoint in U.  Returns (Y, B):
    Y subset of B = V - giant(G-U), G[Y] and G-Y connected, e and f in the boundary of Y."""
    B, giant = giant_complement(G, U)
    m = lambda x: x if x in B else GAM
    GB = nx.Graph()
    GB.add_nodes_from(B)
    for u in B:
        for v in G.neighbors(u):
            GB.add_edge(u, m(v))
    ee = (m(e[0]), m(e[1]))
    ff = (m(f[0]), m(f[1]))
    assert ee[0] != ee[1] and ff[0] != ff[1]
    if GB.number_of_nodes() == 2:  # B = {u}
        return set(B), B
    Y0, Z0 = two_paths_bond(GB, ee, ff)
    sides = extend_partition(GB, Y0, Z0, prefer=GAM if giant else None)
    Y = sides[1] if GAM in sides[0] else sides[0]
    return set(Y), B


def is_cobond_side(G, Y):
    Y = set(Y)
    if not Y or len(Y) == G.number_of_nodes():
        return False
    return nx.is_connected(G.subgraph(Y)) and nx.is_connected(G.subgraph(set(G) - Y))


def crosses(e, Y):
    return (e[0] in Y) != (e[1] in Y)


def check_A(rng):
    print("[A] Lemma lem:bond: every two edges of a 2-connected graph lie on a common bond (constructive)", flush=True)
    ngraph = npairs = 0
    t0 = time.time()
    for trial in range(400):
        n = rng.randint(3, 11)
        p = rng.uniform(0.25, 0.7)
        G = nx.gnp_random_graph(n, p, seed=rng.randrange(10 ** 9))
        if not nx.is_biconnected(G):
            continue
        ngraph += 1
        E = list(G.edges())
        brute = None
        if n <= 9:
            brute = set()
            nodes = list(G)
            for mask in range(1, 2 ** (n - 1)):  # vertex 0 never in Y: each bond once
                Y = {nodes[i + 1] for i in range(n - 1) if mask >> i & 1}
                if is_cobond_side(G, Y):
                    cut = frozenset(ekey(*x) for x in E if crosses(x, Y))
                    brute.add(cut)
        for e, f in itertools.combinations(E, 2):
            npairs += 1
            Y0, Z0 = two_paths_bond(G, e, f)
            sides = extend_partition(G, Y0, Z0)
            Y = sides[0]
            ok = is_cobond_side(G, Y) and crosses(e, Y) and crosses(f, Y)
            check(ok, f"CB construction failed on n={n}, e={e}, f={f}")
            if brute is not None:
                check(any(ekey(*e) in c and ekey(*f) in c for c in brute), f"CB brute force: no common bond {e} {f}")
    # multigraph case: parallel edges represented by the same endpoint pair
    for trial in range(100):
        G = nx.cycle_graph(rng.randint(3, 8))
        e = rng.choice(list(G.edges()))
        Y0, Z0 = two_paths_bond(G, e, e)
        Y = extend_partition(G, Y0, Z0)[0]
        check(is_cobond_side(G, Y) and crosses(e, Y), "CB parallel-edge case")
    print(f"  {ngraph} random 2-connected graphs (n <= 11), {npairs} edge pairs: construction valid; "
          f"brute force agrees for n <= 9  ({time.time() - t0:.1f}s)", flush=True)


# ===================================================================== PP construction
def incidence(P):
    inc = defaultdict(list)
    for e in P:
        inc[e[0]].append(e)
        inc[e[1]].append(e)
    return inc


def ball(G, S, r):
    seen = set(S)
    fr = list(S)
    for _ in range(r):
        nx_ = []
        for x in fr:
            for y in G.neighbors(x):
                if y not in seen:
                    seen.add(y)
                    nx_.append(y)
        fr = nx_
    return seen


def connector(G, starts, ends, blocked):
    """Shortest path in G - blocked from a vertex of `starts` to a vertex of `ends`."""
    starts = [s for s in starts if s not in blocked]
    hit = next((s for s in starts if s in ends), None)
    par = {s: None for s in starts}
    dq = deque(starts)
    while dq and hit is None:
        x = dq.popleft()
        for y in G.neighbors(x):
            if y in blocked or y in par:
                continue
            par[y] = x
            if y in ends:
                hit = y
                break
            dq.append(y)
    if hit is None:
        return None
    path = [hit]
    while par[path[-1]] is not None:
        path.append(par[path[-1]])
    return path[::-1]


def build_family(G, P, Q, rng, mind=2, block=None):
    """Greedy robust connectors (blocked = block-balls of previous gadgets, default block = mind),
    bond gadgets via Lemma CB.  Discards a connector whose gadget comes within distance mind - 1
    of an earlier gadget, so gadgets are pairwise at distance >= mind (the proof uses
    mind = sigma + 1 and block = sigma + r0, both far larger than here)."""
    block = mind if block is None else block
    Pinc, Qinc = incidence(P), incidence(Q)
    VQ = set(Qinc)
    blocked, dead, near = set(), set(), set()
    fam = []
    while True:
        starts = [x for x in Pinc if x not in dead]
        path = connector(G, starts, VQ - blocked, blocked)
        if path is None:
            break
        x0, xk = path[0], path[-1]
        e = rng.choice(Pinc[x0])
        f = rng.choice([g for g in Qinc[xk] if g != e])
        Y, B = bond_gadget(G, path, e, f)
        if Y & near:
            dead.add(x0)
            continue
        fam.append(dict(Y=Y, e=e, f=f, U=path, B=B))
        near |= ball(G, Y, mind - 1)
        blocked |= ball(G, Y, block)
    return fam


def family_stats(G, fam, P, Q, mind=2):
    Pset, Qset = set(map(lambda x: ekey(*x), P)), set(map(lambda x: ekey(*x), Q))
    W = set()
    ok = True
    for g in fam:
        Y = g["Y"]
        ok &= nx.is_connected(G.subgraph(Y)) and nx.is_connected(G.subgraph(set(G) - Y))
        ok &= ekey(*g["e"]) in Pset and ekey(*g["f"]) in Qset and g["e"] != g["f"]
        ok &= crosses(g["e"], Y) and crosses(g["f"], Y)
        ok &= Y <= g["B"]
        ok &= not (W & ball(G, Y, mind - 1))
        W |= Y
    return ok, W


def lp_bound(G, fam, lam):
    D = max(d for _, d in G.degree())
    W = set().union(*[g["Y"] for g in fam])
    out = []
    for g in fam:
        Y = g["Y"]
        dY = sum(1 for y in Y for z in G.neighbors(y) if z not in Y)
        simple = all(sum(1 for z in G.neighbors(x) if z in Y) <= 1 for x in ball(G, Y, 1) - Y)
        kap = (1 if simple else 2) / lam
        out.append(D ** (-(len(Y) - 1)) / dY * (1 + kap) ** (-(dY - 1)))
    return out


def pq_cases(G, rng, units=None):
    E = [ekey(*e) for e in G.edges()]
    M = len(E)
    cases = []
    for delta in (1 / 3, 1 / 9):
        P, Q = far_sets(G, E, delta, rng)
        cases.append((f"far 1/{round(1 / delta)}", P, Q))
    for r in (2, 3):
        vts = [geometric_vtree(G, E, rng.choice(list(G.nodes())), rng) for _ in range(r)]
        Pi, Qi = refine([vtree_nodes(T) for T in vts], list(range(M)))
        cases.append((f"vtree r={r}", [E[i] for i in Pi], [E[i] for i in Qi]))
    if units is not None:
        # adversarial: P = pocket edges c-a on the near half, Q = pocket edges c-b on the far half
        root = rng.choice(list(G.nodes()))
        dist = nx.single_source_shortest_path_length(G, root)
        us = sorted(units, key=lambda t: dist[t[4]])
        h = len(us) // 2
        cases.append(("pocket-adv", [ekey(t[4], t[1]) for t in us[:h]], [ekey(t[4], t[2]) for t in us[h:]]))
    return cases


def check_BC(rng, sizes, minds=(2, 4, 6)):
    print("\n[B,C] Theorem thm:packing construction: robust connectors -> bond gadgets (Lemma lem:bond), gadgets pairwise at distance >= d", flush=True)
    print("  per row: cand/N at d=2; abs = fraction of connectors U with B(U) != U (absorption needed); |U|, |Y| max;"
          " then for each d: gadgets/N, G-W connected?, lambda2(G-W)/lambda2(G), min Lemma lem:pendant-prob bound (kappa' = 2 kappa)", flush=True)
    rows = []
    for N0 in sizes:
        ev = lambda x: max(10, int(round(x / 2)) * 2)
        Gp, units = pocket_graph(connected_regular(3, ev(N0 / 5.5), rng.randrange(10 ** 6)))
        graphs = [("cubic", connected_regular(3, N0, rng.randrange(10 ** 6)), None),
                  ("4-reg", connected_regular(4, N0, rng.randrange(10 ** 6)), None),
                  ("pocket", Gp, units),
                  ("trunc", nx.convert_node_labels_to_integers(truncate(connected_regular(3, ev(N0 / 3), rng.randrange(10 ** 6)))), None),
                  ("prism", prism(connected_regular(3, ev(N0 / 2), rng.randrange(10 ** 6))), None)]
        for gname, G, un in graphs:
            check(nx.is_biconnected(G), f"{gname} N={G.number_of_nodes()} is 2-connected")
            N = G.number_of_nodes()
            l2G = lam2(G)
            D = max(d for _, d in G.degree())
            for label, P, Q in pq_cases(G, rng, un):
                out = []
                for mind in minds:
                    fam = build_family(G, P, Q, rng, mind=mind)
                    ok, W = family_stats(G, fam, P, Q, mind)
                    check(ok, f"{gname} N={N} {label} d={mind}: gadgets valid (connected, co-connected, marks cross, separated)")
                    hl = l2G / 2  # Cheeger: h(G) >= lambda2/2, so D/h <= 2D/lambda2
                    check(all(len(g["B"]) <= len(g["U"]) * (1 + D / hl) + 1e-9 for g in fam),
                          f"{gname} N={N} {label}: |B(U)| <= |U|(1+D/h)")
                    H = G.subgraph([v for v in G if v not in W])
                    conn = nx.is_connected(H)
                    l2 = lam2(H) if conn else 0.0
                    pmin = min(lp_bound(G, fam, l2)) if (conn and fam) else float("nan")
                    out.append((mind, len(fam) / N, conn, l2 / l2G, pmin, fam))
                fam2 = out[0][5]
                absorbed = np.mean([len(g["B"]) > len(g["U"]) for g in fam2]) if fam2 else 0.0
                Umax = max((len(g["U"]) for g in fam2), default=0)
                Ymax = max((len(g["Y"]) for g in fam2), default=0)
                Ymean = np.mean([len(g["Y"]) for g in fam2]) if fam2 else 0.0
                rows.append((gname, N, label, absorbed, out))
                print(f"  {gname:6s} N={N:5d} {label:10s} abs={absorbed:.2f} |U|max={Umax:2d} |Y|max={Ymax:2d} mean={Ymean:5.2f} | "
                      + " | ".join(f"d={m}: {c:.4f} {'conn' if cn else 'DISC'} {r:.2f} {pm:.0e}" for m, c, cn, r, pm, _ in out),
                      flush=True)
    # summary per graph type and d
    print("  summary (all sizes and P,Q sources): fraction of instances with G-W connected, median lambda2 ratio", flush=True)
    for gname in ("cubic", "4-reg", "pocket", "trunc", "prism"):
        parts = []
        for k, mind in enumerate(minds):
            rs = [o[k] for (gn, _, _, _, o) in rows if gn == gname]
            parts.append(f"d={mind}: conn {np.mean([x[2] for x in rs]):.2f}, l2 ratio {np.median([x[3] for x in rs]):.2f}, gad/N {np.median([x[1] for x in rs]):.4f}")
        print(f"    {gname:6s} " + " | ".join(parts), flush=True)
    return rows


# ===================================================================== D. exact checks
def exact_check(name, G, fam, rng):
    n = G.number_of_nodes()
    E = [ekey(*e) for e in G.edges()]
    idx = {e: i for i, e in enumerate(E)}
    W = set().union(*[g["Y"] for g in fam])
    H = G.subgraph([v for v in G if v not in W])
    if not nx.is_connected(H) or len(fam) < 2:
        return None
    lam = lam2(H)
    D = max(d for _, d in G.degree())
    info = []
    for g in fam:
        Y = g["Y"]
        bd = sorted(idx[ekey(y, z)] for y in Y for z in G.neighbors(y) if z not in Y)
        inner = sorted(idx[ekey(*x)] for x in G.subgraph(Y).edges())
        simple = all(sum(1 for z in G.neighbors(x) if z in Y) <= 1 for x in ball(G, Y, 1) - Y)
        dY = len(bd)
        info.append(dict(Y=Y, bd=bd, bdmask=sum(1 << i for i in bd), allmask=sum(1 << i for i in bd + inner),
                         e=idx[ekey(*g["e"])], f=idx[ekey(*g["f"])], simple=simple,
                         p2=D ** (-(len(Y) - 1)) / dY * (1 + 2 / lam) ** (-(dY - 1)),
                         p1=D ** (-(len(Y) - 1)) / dY * (1 + 1 / lam) ** (-(dY - 1))))
    m = len(fam)
    # (1) pendant restriction lemma: rank 2^m, value = prod XOR
    ones = [idx[ekey(*x)] for x in nx.minimum_spanning_tree(H).edges()]
    for g in fam:
        ones += [idx[ekey(*x)] for x in nx.minimum_spanning_tree(G.subgraph(g["Y"])).edges()]
    free_e = [d["e"] for d in info]
    free_f = [d["f"] for d in info]

    def is_tree(chosen):
        if len(chosen) != n - 1:
            return False
        par = list(range(n))

        def fd(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for i in chosen:
            a, b = fd(E[i][0]), fd(E[i][1])
            if a == b:
                return False
            par[a] = b
        return True
    Mat = np.zeros((2 ** m, 2 ** m), dtype=int)
    ok_xor = True
    for xe in range(2 ** m):
        for xf in range(2 ** m):
            ch = ones + [free_e[i] for i in range(m) if xe >> i & 1] + [free_f[i] for i in range(m) if xf >> i & 1]
            v = int(is_tree(ch))
            Mat[xe, xf] = v
            ok_xor &= v == int(all((xe >> i & 1) ^ (xf >> i & 1) for i in range(m)))
    check(ok_xor, f"{name}: pendant restriction = prod XOR")
    check(np.linalg.matrix_rank(Mat) == 2 ** m, f"{name}: restricted rank 2^m")
    # (2) blocks
    trees = [int(t) for t in enum_trees(n, E)]
    treeset = set(trees)
    pairs = [(2 * i, 2 * i + 1) for i in range(m // 2)]

    def pat(t, k):
        s = t & info[k]["bdmask"]
        return "e" if s == 1 << info[k]["e"] else ("f" if s == 1 << info[k]["f"] else None)

    def good(t, q):
        return (pat(t, q[0]), pat(t, q[1])) in (("e", "f"), ("f", "e"))
    special = {k: (1 << info[k]["e"]) | (1 << info[k]["f"]) for k in range(m)}
    blocks = defaultdict(list)
    for t in trees:
        Gq = tuple(qi for qi, q in enumerate(pairs) if good(t, q))
        F = t
        for qi in Gq:
            F &= ~(special[pairs[qi][0]] | special[pairs[qi][1]])
        blocks[(Gq, F)].append(t)
    owner = {}
    ok_cube = ok_perm = ok_disj = True
    for key, ts in blocks.items():
        Gq, F = key
        g = len(Gq)
        ok_cube &= len(ts) == 2 ** g
        rows = [sum(1 << info[pairs[qi][b]]["e"] for qi, b in zip(Gq, bits)) for bits in itertools.product((0, 1), repeat=g)]
        cols = [sum(1 << info[pairs[qi][b]]["f"] for qi, b in zip(Gq, bits)) for bits in itertools.product((0, 1), repeat=g)]
        Mt = np.zeros((len(rows), len(cols)), dtype=int)
        for a, rs in enumerate(rows):
            for b, cs in enumerate(cols):
                cell = F | rs | cs
                Mt[a, b] = cell in treeset
                if owner.get(cell, key) != key:
                    ok_disj = False
                owner[cell] = key
        ok_perm &= bool((Mt.sum(0) == 1).all() and (Mt.sum(1) == 1).all())
    check(ok_cube, f"{name}: blocks are cubes")
    check(ok_perm, f"{name}: block rectangles are permutations")
    check(ok_disj, f"{name}: block rectangles disjoint")
    # (3) Lemma LP, conditional on everything inside and on the boundary of the other gadgets
    worst2 = worst1 = np.inf
    for k in range(m):
        others = [info[j]["allmask"] for j in range(m) if j != k]
        groups = defaultdict(lambda: defaultdict(int))
        for t in trees:
            groups[tuple(t & mm for mm in others)][t & info[k]["bdmask"]] += 1
        for cnt in groups.values():
            tot = sum(cnt.values())
            for b in info[k]["bd"]:
                pr = cnt.get(1 << b, 0) / tot
                worst2 = min(worst2, pr / info[k]["p2"])
                if info[k]["simple"]:
                    worst1 = min(worst1, pr / info[k]["p1"])
    check(worst2 >= 1 - 1e-9, f"{name}: LP bound kappa'=2kappa (worst ratio {worst2:.3f})")
    if worst1 < np.inf:
        check(worst1 >= 1 - 1e-9, f"{name}: LP bound kappa'=kappa, simple gadgets (worst ratio {worst1:.3f})")
    print(f"  {name}: N={n} M={len(E)} tau={len(trees)} m={m} sizes={[len(d['Y']) for d in info]} "
          f"simple={[d['simple'] for d in info]} lambda2(G-W)={lam:.3f} blocks={len(blocks)} "
          f"LP ratio (2k)={worst2:.2f} (k, simple)={worst1:.2f}", flush=True)
    return True


def check_D(rng):
    print("\n[D] exact checks on small graphs: restriction rank, blocks, Lemma lem:pendant-prob with absorbed blobs", flush=True)
    done = absorbed_seen = 0
    tries = 0
    graphs = []
    # pocket graphs on small bases (absorption happens when a connector crosses a unit a-b edge)
    for base in [nx.complete_graph(4), nx.complete_bipartite_graph(2, 3), nx.cycle_graph(4)]:
        Gp, units = pocket_graph(nx.convert_node_labels_to_integers(base))
        graphs.append(("pocket(" + str(base.number_of_nodes()) + "," + str(base.number_of_edges()) + ")", nx.convert_node_labels_to_integers(Gp)))
    graphs += [("petersen", nx.petersen_graph()), ("trunc(K4)", truncate(nx.complete_graph(4))),
               ("rr3_14", connected_regular(3, 14, 5)), ("rr3_16", connected_regular(3, 16, 9))]
    for gname, G in graphs:
        G = nx.convert_node_labels_to_integers(G)
        E = [ekey(*e) for e in G.edges()]
        if len(E) > 62:
            continue
        L = nx.laplacian_matrix(G).toarray().astype(float)
        tau = round(math.exp(np.linalg.slogdet(L[1:, 1:])[1]))
        if tau > 3_000_000:
            print(f"  {gname}: tau={tau} too large, skipped")
            continue
        cnt = 0
        for trial in range(60):
            tries += 1
            # random P, Q: edges near/far from a root, smallish delta so connectors are paths
            delta = rng.choice([0.2, 0.3, 0.4])
            P, Q = far_sets(G, E, delta, rng)
            fam = build_family(G, P, Q, rng, mind=2, block=1)
            if len(fam) < 2:
                continue
            fam = fam[: min(len(fam), 4)]
            fam = fam[: len(fam) - len(fam) % 2]
            ok, W = family_stats(G, fam, P, Q)
            check(ok, f"{gname}: family valid")
            ab = any(len(g["B"]) > len(g["U"]) or len(g["Y"]) != len(g["U"]) for g in fam)
            if cnt >= 2 and not ab:
                continue
            r = exact_check(f"{gname} trial {trial}{' [absorbed/blob]' if ab else ''}", G, fam, rng)
            if r:
                done += 1
                cnt += 1
                absorbed_seen += ab
            if cnt >= 4:
                break
    check(done >= 10, f"[D] enough exact instances ({done})")
    print(f"  exact instances: {done}, with absorbed or non-path gadgets: {absorbed_seen}", flush=True)


# ===================================================================== E. constants
def mu(Delta, i):
    return Delta - math.sqrt(Delta ** 2 - i ** 2)


def constants(D, h, delta):
    """log10 of the explicit constants of Theorems thm:packing and thm:signed-sums-linear, deterministic route."""
    beta = 2 * delta / D
    tstar = math.floor(math.log(1 / beta) / math.log(1 + h / (2 * D))) + 1
    ell = 2 * tstar + 1
    L = math.floor(ell * (1 + D / h))
    Dp = D * L
    hp = h / (Dp + 1)
    lam0 = mu(D, hp)
    sigma = 2 * math.ceil((Dp + 1) / h)
    r0 = math.floor(D * ell / h)
    s = sigma + r0
    log10_c = math.log10(h * beta / (8 * L)) - (s + 1) * math.log10(D)
    log10_p = -(L - 1) * math.log10(D) - math.log10(Dp) - (Dp - 1) * math.log10(1 + 2 / lam0)
    return dict(ell=ell, L=L, s=s, sigma=sigma, log10_lam0=math.log10(lam0), log10_c=log10_c, log10_p=log10_p,
                log10_g=2 * log10_p + log10_c - math.log10(2))


def check_E():
    print("\n[E] explicit constants (log10): gadgets >= c N; per-gadget pendant bound p; approximate exponent g ~ p^2 c N / 2;"
          " exact exponent = c N", flush=True)
    for D, h in [(3, 0.1), (3, 0.2), (4, 0.4), (6, 0.76)]:
        for r in (1, 2, 3, 5):
            C = constants(D, h, 3.0 ** (-r))
            print(f"  D={D} h={h:<4} r={r}: connector <= {C['ell']} vertices, gadget <= {C['L']}, sigma={C['sigma']}, blocked radius {C['s']}, "
                  f"log10 lambda0={C['log10_lam0']:.3g}, log10 c={C['log10_c']:.4g}, log10 p={C['log10_p']:.4g}, "
                  f"log10(g/N)={C['log10_g']:.4g}", flush=True)
            check(math.isfinite(C["log10_g"]), "constants finite")
    # monotone in r (log of the constant is O(r log r) for fixed D, h)
    for key in ("log10_c", "log10_g"):
        xs = [-constants(3, 0.2, 3.0 ** (-r))[key] for r in range(1, 9)]
        print(f"  D=3 h=0.2: -{key} for r=1..8: " + ", ".join(f"{x:.4g}" for x in xs)
              + "; increments " + ", ".join(f"{b - a:.4g}" for a, b in zip(xs, xs[1:])), flush=True)


# ===================================================================== F. elementary lemmas with exact h
def exact_h(G):
    n = G.number_of_nodes()
    masks = np.arange(1, 2 ** n, dtype=np.int64)
    size = np.zeros(len(masks), dtype=np.int64)
    for v in range(n):
        size += (masks >> v) & 1
    cut = np.zeros(len(masks), dtype=np.int64)
    for u, v in G.edges():
        cut += ((masks >> u) ^ (masks >> v)) & 1
    ok = size <= n // 2
    return float((cut[ok] / size[ok]).min())


def check_F(rng):
    print("\n[F] Lemmas lem:connector, lem:blob, lem:sep-delete (one gadget) against exact h, N <= 18", flush=True)
    worst_rc = worst_gc = 0.0
    n_rc = n_gc = 0
    for trial in range(40):
        n = rng.choice([12, 14, 16, 18])
        d = rng.choice([3, 3, 4])
        G = nx.convert_node_labels_to_integers(connected_regular(d, n, rng.randrange(10 ** 6)))
        if rng.random() < 0.3:
            G = nx.convert_node_labels_to_integers(pocket_graph(nx.cycle_graph(4))[0])
        n = G.number_of_nodes()
        if n > 22:
            continue
        h = exact_h(G)
        D = max(dd for _, dd in G.degree())
        dist = dict(nx.all_pairs_shortest_path_length(G))
        for _ in range(30):
            beta = rng.choice([0.15, 0.25, 0.4])
            k = max(1, math.ceil(beta * n))
            X = set(rng.sample(range(n), k)); Y = set(rng.sample(range(n), k))
            zmax = int(h * beta * n / (4 * D))
            Z = set(rng.sample(range(n), zmax))
            tstar = math.floor(math.log(1 / beta) / math.log(1 + h / (2 * D))) + 1
            path = connector(G, list(X - Z), Y - Z, Z)
            n_rc += 1
            check(path is not None and len(path) <= 2 * tstar + 1, f"RC: connector with <= 2t*+1 vertices (n={n})")
            if path is not None:
                worst_rc = max(worst_rc, len(path) / (2 * tstar + 1))
        for _ in range(30):
            u = rng.randint(1, 3)
            start = rng.randrange(n)
            U = {start}
            while len(U) < u:
                U.add(rng.choice([y for x in U for y in G.neighbors(x) if y not in U]))
            if n <= len(U) * (1 + D / h):
                continue
            B, giant = giant_complement(G, U)
            n_gc += 1
            check(len(giant) > n / 2 and len(B) <= len(U) * (1 + D / h) + 1e-9, f"GC: |B(U)| <= |U|(1+D/h) (n={n})")
            check(all(min(dist[x][y] for y in U) <= math.floor(D * len(U) / h) for x in B), "GC: B(U) within D|U|/h of U")
            worst_gc = max(worst_gc, len(B) / (len(U) * (1 + D / h)))
    print(f"  lem:connector: {n_rc} instances, max (vertices / bound) = {worst_rc:.2f};  lem:blob: {n_gc} instances, max |B|/bound = {worst_gc:.2f}", flush=True)
    # Lemma SD with one gadget (cases j = 0, 1 of the proof): h(G - Y) >= h(G) / (|dY| + 1) for co-connected Y, N > 2|Y|
    worst_sd = np.inf
    n_sd = 0
    for trial in range(60):
        n = rng.choice([10, 12, 14, 16])
        base = connected_regular(rng.choice([3, 3, 4]), n, rng.randrange(10 ** 6))
        if rng.random() < 0.3:
            base = pocket_graph(nx.cycle_graph(rng.choice([3, 4])))[0]
        G = nx.convert_node_labels_to_integers(base)
        n = G.number_of_nodes()
        hG = exact_h(G)
        E = [ekey(*e) for e in G.edges()]
        for _ in range(6):
            x = rng.randrange(n)
            y = rng.randrange(n)
            U = nx.shortest_path(G, x, y)[: rng.randint(1, 4)]
            e = rng.choice([ekey(U[0], z) for z in G.neighbors(U[0])])
            f = rng.choice([g for g in (ekey(U[-1], z) for z in G.neighbors(U[-1])) if g != e])
            Y, B = bond_gadget(G, U, e, f)
            if not (n > 2 * len(Y)):
                continue
            check(is_cobond_side(G, Y) and crosses(e, Y) and crosses(f, Y), "SD: CB gadget valid")
            dY = sum(1 for a in Y for b in G.neighbors(a) if b not in Y)
            H = nx.convert_node_labels_to_integers(G.subgraph([v for v in G if v not in Y]))
            hH = exact_h(H)
            n_sd += 1
            check(hH >= hG / (dY + 1) - 1e-12, f"SD: h(G-Y) >= h(G)/(|dY|+1) (n={n}, |Y|={len(Y)})")
            Dm = max(d for _, d in G.degree())
            check(lam2(H) >= mu(Dm, hH) - 1e-9, "imported Cheeger-Mohar: lambda2(H) >= mu(D, h(H))")
            worst_sd = min(worst_sd, hH / (hG / (dY + 1)))
    print(f"  lem:sep-delete (one gadget): {n_sd} instances, min h(G-Y) / (h(G)/(|dY|+1)) = {worst_sd:.2f}; Cheeger-Mohar lambda2 >= mu(D,h) held", flush=True)


def main():
    rng = random.Random(20261104)
    t0 = time.time()
    check_A(rng)
    sizes = [int(x) for x in os.environ.get("PP_SIZES", "200,800,2400").split(",")]
    check_BC(rng, sizes)
    check_D(rng)
    check_E()
    check_F(rng)
    print(f"\ntotal {time.time() - t0:.0f}s", flush=True)
    if FAIL:
        print(f"\n{len(FAIL)} CHECK(S) FAILED")
        sys.exit(1)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()

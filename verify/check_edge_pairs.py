# /// script
# dependencies = ["numpy", "scipy", "networkx"]
# ///
"""Checks for the 2x2 lemma, the hub graphs, the edge-pair chain and Conjecture conj:pairs.

  P1  Lemma lem:two-by-two-general on random kernels on {0,1}^2: the norm formula, Delta <= 1/2, the rare-cell
      bound (a) and the covariance identity and bounds (b)
  P2  Lemmas lem:edge-pair-law and lem:edge-pair-res for every pair of edges at a common vertex of seven small
      graphs (one with parallel edges), against the exact law from all spanning trees
  P3  Theorem thm:edge-pair on small graphs, exactly: for every past sigma of positive probability the conditional
      law of the two bits equals the transfer-current law in K_sigma (Lemma lem:markov), Delta >= B_sigma, and
      beta <= ||sqrt(law of the pair bits)|| <= prod_k (1 - Delta_k^2)^(1/2)
  P4  Hub graphs: Lemma lem:subdivide on random subdivisions (exact h), Proposition prop:hub-expansion on the two
      smallest hub graphs (exact h), Proposition prop:hub (Pr[a_w = 1] <= L^-a), and Proposition prop:hub-pairs on
      G_{L,13}(2,2) with random cubic cores and L from the hypothesis, for sampled pasts
  P5  The example after Conjecture conj:pairs (a triangle attached at each vertex): beta = 1, Y = 0, h, degree
  P6  Numerical evidence for Conjecture conj:pairs on random regular graphs, D in {3,4}, N in {200,1000}

Exits 1 on any failure.
"""
import itertools
import math
import sys

import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from common import enum_trees

FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL:", msg, flush=True)


# ------------------------------------------------------------------ 2x2 kernels
def delta_of(K):
    K = np.clip(K, 0, None)
    return abs(math.sqrt(K[0, 0] * K[1, 1]) - math.sqrt(K[0, 1] * K[1, 0]))


def p1_two_by_two(rng):
    worst = {"formula": 0.0, "upper": -1.0, "rare": -1.0, "cov_id": 0.0, "cov_lb": -1.0}
    cnt = 0
    for t in range(20000):
        alpha = rng.choice([0.1, 0.5, 1.0, 5.0])
        K = rng.dirichlet(np.full(4, alpha))
        if t % 5 == 0:
            K[rng.integers(4)] = 0.0
            K /= K.sum()
        K = K.reshape(2, 2)
        D = delta_of(K)
        nrm = np.linalg.norm(np.sqrt(K), 2) ** 2
        worst["formula"] = max(worst["formula"], abs(nrm - (1 + math.sqrt(max(0.0, 1 - 4 * D * D))) / 2))
        worst["upper"] = max(worst["upper"], nrm - (1 - D * D), D - 0.5)
        p, pp, dl = K[1, 0] * rng.uniform(0.5, 1), K[0, 1] * rng.uniform(0.5, 1), K[1, 1] * rng.uniform(1, 2)
        worst["rare"] = max(worst["rare"], max(0.0, math.sqrt(p * pp) - math.sqrt(dl)) - D)
        px, py, pxy = K[1].sum(), K[:, 1].sum(), K[1, 1]
        cov = pxy - px * py
        worst["cov_id"] = max(worst["cov_id"], abs(K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0] - cov))
        sx = math.sqrt((K[0, 0] + K[0, 1]) * (K[1, 0] + K[1, 1]))
        sy = math.sqrt((K[0, 0] + K[1, 0]) * (K[0, 1] + K[1, 1]))
        cov_det = K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]  # equal to cov (checked above), without cancellation
        lb = abs(cov_det) / min(sx, sy) if min(sx, sy) > 0 else 0.0
        ssum = math.sqrt(K[0, 0] * K[1, 1]) + math.sqrt(K[0, 1] * K[1, 0])
        Ds = abs(K[0, 0] * K[1, 1] - K[0, 1] * K[1, 0]) / ssum if ssum > 0 else 0.0  # = Delta, without cancellation
        worst["cov_lb"] = max(worst["cov_lb"], (lb - Ds) / max(Ds, 1e-300) if lb > 0 else 0.0, 2 * abs(cov_det) - lb)
        cnt += 1
    print(f"P1 {cnt} random kernels: norm formula error {worst['formula']:.1e}; largest excess of "
          f"||sqrt K||^2 - (1 - Delta^2) and Delta - 1/2: {worst['upper']:.1e}; of the rare-cell bound over Delta: "
          f"{worst['rare']:.1e}; covariance identity error {worst['cov_id']:.1e}; largest excess of "
          f"|Cov|/min(sigma) over Delta and of 2|Cov| over it (relative): {worst['cov_lb']:.1e}")
    check(worst["formula"] < 1e-12 and worst["upper"] < 1e-12 and worst["rare"] < 1e-12
          and worst["cov_id"] < 1e-15 and worst["cov_lb"] < 1e-12, "P1 two-by-two-general")


# ------------------------------------------------------------------ electrical helpers
def laplacian(n, E):
    L = np.zeros((n, n))
    for u, v in E:
        if u != v:
            L[u, u] += 1
            L[v, v] += 1
            L[u, v] -= 1
            L[v, u] -= 1
    return L


def resistance(n, E, x, y):
    """Effective resistance between x and y in the multigraph (range(n), E); inf if disconnected; 0 if x = y."""
    if x == y:
        return 0.0
    G = nx.MultiGraph()
    G.add_nodes_from(range(n))
    G.add_edges_from((u, v) for u, v in E if u != v)
    comp = nx.node_connected_component(G, x)
    if y not in comp:
        return math.inf
    comp = sorted(comp)
    idx = {v: i for i, v in enumerate(comp)}
    Es = [(idx[u], idx[v]) for u, v in E if u in idx and v in idx]
    k = len(comp)
    M = laplacian(k, Es) + np.ones((k, k)) / k
    b = np.zeros(k)
    b[idx[x]], b[idx[y]] = 1.0, -1.0
    return float(b @ np.linalg.solve(M, b))


def pair_quantities(n, E, e, f):
    """Y, Pr[e in T], Pr[f in T] for edges e, f with a common endpoint, in the multigraph (range(n), E)."""
    (a, b), (c, d) = E[e], E[f]
    w = ({a, b} & {c, d}).pop() if a != b and c != d else None
    u = a if b == w else b
    v = c if d == w else d
    Ruw, Rvw, Ruv = resistance(n, E, u, w), resistance(n, E, v, w), resistance(n, E, u, v)
    return 0.5 * (Ruw + Rvw - Ruv), Ruw, Rvw, (u, v, w)


def law_from(pe, pf, Y):
    K11 = pe * pf - Y * Y
    return np.array([[1 - pe - pf + K11, pf - K11], [pe - K11, K11]])


def exact_pair_law(T, e, f):
    x = ((T >> np.uint64(e)) & np.uint64(1)).astype(bool)
    y = ((T >> np.uint64(f)) & np.uint64(1)).astype(bool)
    N = len(T)
    return np.array([[np.sum(~x & ~y), np.sum(~x & y)], [np.sum(x & ~y), np.sum(x & y)]]) / N


# ------------------------------------------------------------------ P2
def small_graphs(rng):
    out = {"K4": nx.complete_graph(4), "K33": nx.complete_bipartite_graph(3, 3), "Petersen": nx.petersen_graph(),
           "prism C5": nx.circular_ladder_graph(5), "cube": nx.hypercube_graph(3), "wheel W6": nx.wheel_graph(7)}
    gs = {}
    for name, G in out.items():
        G = nx.convert_node_labels_to_integers(G)
        gs[name] = (G.number_of_nodes(), [tuple(sorted(e)) for e in G.edges()])
    n, E = gs["K4"]
    gs["K4 + parallel edge"] = (n, E + [E[0], E[3]])
    return gs


def p2_edge_pair_law(rng):
    worst = {"cov": 0.0, "Yneg": 0.0, "D2": -1.0, "Dq": -1.0, "res": -1.0, "law": 0.0}
    cnt = 0
    for name, (n, E) in small_graphs(rng).items():
        T = enum_trees(n, E)
        for e, f in itertools.combinations(range(len(E)), 2):
            ends_e, ends_f = set(E[e]), set(E[f])
            if not ends_e & ends_f:
                continue
            for w in ends_e & ends_f:
                u = (ends_e - {w}).pop() if len(ends_e) == 2 else w
                v = (ends_f - {w}).pop() if len(ends_f) == 2 else w
                if u == w or v == w:
                    continue
                Ruw, Rvw, Ruv = resistance(n, E, u, w), resistance(n, E, v, w), resistance(n, E, u, v)
                Y = 0.5 * (Ruw + Rvw - Ruv)
                K = exact_pair_law(T, e, f)
                pe, pf = K[1].sum(), K[:, 1].sum()
                cov = K[1, 1] - pe * pf
                worst["cov"] = max(worst["cov"], abs(cov + Y * Y))
                worst["law"] = max(worst["law"], abs(pe - Ruw), abs(pf - Rvw))
                worst["Yneg"] = min(worst["Yneg"], Y)
                D = delta_of(K)
                worst["D2"] = max(worst["D2"], 2 * Y * Y - D)
                qe, qf = 1 - pe, 1 - pf
                if qe * qf > 0:
                    worst["Dq"] = max(worst["Dq"], Y * Y / (2 * math.sqrt(qe * qf)) - D)
                Em = [E[i] for i in range(len(E)) if i != e]
                Ef = [E[i] for i in range(len(E)) if i != f]
                rho = min(resistance(n, Em, u, w), resistance(n, Ef, v, w))
                Ew = [x for x in E if w not in x]
                Rp = resistance(n, Ew, u, v)
                a = 1.0 if math.isinf(rho) else rho / (1 + rho)
                b = 1.0 if math.isinf(Rp) else Rp / (2 + Rp)
                worst["res"] = max(worst["res"], (a - b) - Y)
                cnt += 1
    print(f"P2 {cnt} pairs of edges at a common vertex: |Cov + Y^2| <= {worst['cov']:.1e}, "
          f"Kirchhoff error {worst['law']:.1e}, smallest Y {worst['Yneg']:.3f}; largest excess of 2Y^2 over Delta "
          f"{worst['D2']:.1e}, of Y^2/(2 sqrt(q_e q_f)) over Delta {worst['Dq']:.1e}, of the bound of "
          f"lem:edge-pair-res over Y {worst['res']:.1e}")
    check(worst["cov"] < 1e-12 and worst["law"] < 1e-12 and worst["Yneg"] > -1e-12 and worst["D2"] < 1e-12
          and worst["Dq"] < 1e-12 and worst["res"] < 1e-12, "P2 edge-pair law and resistance criterion")


# ------------------------------------------------------------------ P3
def contract(n, E, incl, excl):
    """K_sigma: contract incl, delete excl. Returns (m, Es, cls) with Es indexed like E (None for removed edges);
    loops are kept as (c, c)."""
    uf = nx.utils.UnionFind(range(n))
    for i in incl:
        uf.union(*E[i])
    roots = sorted({uf[x] for x in range(n)})
    idx = {r: k for k, r in enumerate(roots)}
    cls = [idx[uf[x]] for x in range(n)]
    Es = [None if (i in incl or i in excl) else (cls[E[i][0]], cls[E[i][1]]) for i in range(len(E))]
    return len(roots), Es, cls


def opnorm_sqrt_joint(rows, cols):
    _, ri = np.unique(rows, return_inverse=True)
    _, ci = np.unique(cols, return_inverse=True)
    M = np.zeros((ri.max() + 1, ci.max() + 1))
    np.add.at(M, (ri, ci), 1.0 / len(rows))
    return float(np.linalg.norm(np.sqrt(M), 2))


def p3_chain(rng):
    graphs = {"Petersen": nx.petersen_graph(), "prism C6": nx.circular_ladder_graph(6),
              "cube": nx.hypercube_graph(3), "K33": nx.complete_bipartite_graph(3, 3),
              "Moebius ladder 10": nx.circulant_graph(10, [1, 5])}
    tot_pasts = 0
    for name, G in graphs.items():
        G = nx.convert_node_labels_to_integers(G)
        n = G.number_of_nodes()
        E = [tuple(sorted(e)) for e in G.edges()]
        M = len(E)
        T = enum_trees(n, E)
        for trial in range(3):
            inA = np.zeros(M, bool)
            inA[rng.choice(M, M // 2, replace=False)] = True
            used, pairs = set(), []
            for w in rng.permutation(n):
                w = int(w)
                ea = [i for i in range(M) if w in E[i] and inA[i] and i not in used]
                eb = [i for i in range(M) if w in E[i] and not inA[i] and i not in used]
                if ea and eb:
                    e, f = int(rng.choice(ea)), int(rng.choice(eb))
                    pairs.append((w, e, f))
                    used |= {e, f}
                if len(pairs) == 4:
                    break
            bitsA = np.zeros(len(T), np.int64)
            bitsB = np.zeros(len(T), np.int64)
            prod = 1.0
            for k, (w, e, f) in enumerate(pairs):
                past = [x for (_, a, b) in pairs[:k] for x in (a, b)]
                pat = np.zeros(len(T), np.int64)
                for j, x in enumerate(past):
                    pat |= (((T >> np.uint64(x)) & np.uint64(1)).astype(np.int64)) << j
                Dk = math.inf
                for s in np.unique(pat):
                    sub = T[pat == s]
                    incl = [x for j, x in enumerate(past) if s >> j & 1]
                    excl = [x for j, x in enumerate(past) if not s >> j & 1]
                    m, Es, cls = contract(n, E, incl, excl)
                    Kx = exact_pair_law(sub, e, f)
                    pe, pf = Kx[1].sum(), Kx[:, 1].sum()
                    loops = Es[e][0] == Es[e][1] or Es[f][0] == Es[f][1]
                    if loops or min(pe, pf, 1 - pe, 1 - pf) < 1e-12:
                        B = 0.0
                    else:
                        Eg = [x for x in Es if x is not None]
                        wc = cls[w]
                        u = Es[e][0] if Es[e][1] == wc else Es[e][1]
                        v = Es[f][0] if Es[f][1] == wc else Es[f][1]
                        Ruw, Rvw, Ruv = resistance(m, Eg, u, wc), resistance(m, Eg, v, wc), resistance(m, Eg, u, v)
                        Y = 0.5 * (Ruw + Rvw - Ruv)
                        Kt = law_from(Ruw, Rvw, Y)
                        check(np.abs(Kt - Kx).max() < 1e-10, f"P3 {name}: K_sigma law differs from the conditional law")
                        B = Y * Y * max(2.0, 1 / (2 * math.sqrt((1 - pe) * (1 - pf))))
                    check(delta_of(Kx) >= B - 1e-12, f"P3 {name}: Delta < B_sigma")
                    Dk = min(Dk, B)
                    tot_pasts += 1
                prod *= 1 - Dk * Dk
                bitsA |= (((T >> np.uint64(e)) & np.uint64(1)).astype(np.int64)) << k
                bitsB |= (((T >> np.uint64(f)) & np.uint64(1)).astype(np.int64)) << k
            Am = np.uint64(sum(1 << i for i in range(M) if inA[i]))
            beta = opnorm_sqrt_joint(T & Am, T & ~Am)
            nb = opnorm_sqrt_joint(bitsA, bitsB)
            ok = beta <= nb + 1e-12 and nb ** 2 <= prod + 1e-12
            print(f"P3 {name:18s} split {trial}: {len(pairs)} pairs, beta = {beta:.4f} <= {nb:.4f} <= "
                  f"{math.sqrt(prod):.6f}  {'OK' if ok else 'FAIL'}")
            check(ok, f"P3 {name} chain")
    print(f"P3 {tot_pasts} pasts checked against K_sigma")


# ------------------------------------------------------------------ hub graphs
def edge_expansion(n, E):
    """Exact h = min over 0 < |S| <= n/2 of e(S, S^c)/|S| (n <= 24)."""
    assert n <= 24
    eu = np.array([u for u, v in E], np.int64)
    ev = np.array([v for u, v in E], np.int64)
    best = math.inf
    chunk = 1 << 20
    for start in range(1, 1 << n, chunk):
        mk = np.arange(start, min(1 << n, start + chunk), dtype=np.int64)
        size = np.zeros(len(mk), np.int64)
        for i in range(n):
            size += (mk >> i) & 1
        ok = 2 * size <= n
        mk, size = mk[ok], size[ok]
        cut = np.zeros(len(mk), np.int64)
        for a, b in zip(eu, ev):
            cut += ((mk >> a) & 1) ^ ((mk >> b) & 1)
        best = min(best, float((cut / size).min()))
    return best


def cubic_core(k, seed):
    if k == 4:
        return nx.complete_graph(4)
    if k == 6:
        return nx.complete_bipartite_graph(3, 3)
    s = seed
    while True:
        X = nx.random_regular_graph(3, k, seed=s)
        if nx.is_connected(X) and nx.node_connectivity(X) >= 2:
            return X
        s += 1000


def hub_graph(a, b, n, L, XA, XB):
    G = nx.Graph()
    side = {}
    for u, v in XA.edges():
        G.add_edge(("A", u), ("A", v))
        side[frozenset((("A", u), ("A", v)))] = "A"
    for u, v in XB.edges():
        G.add_edge(("B", u), ("B", v))
        side[frozenset((("B", u), ("B", v)))] = "B"
    nxt = {"A": 0, "B": 0}
    first = {}
    for g in range(n):
        hub = ("H", g)
        for s, cnt in (("A", a), ("B", b)):
            for p in range(cnt):
                seq = [hub] + [("P", g, s, p, t) for t in range(1, L)] + [(s, nxt[s])]
                nxt[s] += 1
                for x, y in zip(seq, seq[1:]):
                    G.add_edge(x, y)
                    side[frozenset((x, y))] = s
                first[(g, s, p)] = (seq[0], seq[1])
    assert nxt["A"] == XA.number_of_nodes() and nxt["B"] == XB.number_of_nodes()
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    E = [tuple(sorted((idx[u], idx[v]))) for u, v in G.edges()]
    inA = np.array([side[frozenset((nodes[u], nodes[v]))] == "A" for u, v in E])
    eid = {frozenset(e): i for i, e in enumerate(E)}
    firsts = {key: eid[frozenset((idx[x], idx[y]))] for key, (x, y) in first.items()}
    return G, idx, E, inA, firsts


def logtau(n, E):
    s, ld = np.linalg.slogdet(laplacian(n, E)[1:, 1:])
    return ld if s > 0 else -math.inf


def p4_hub(rng):
    # Lemma lem:subdivide
    worst = math.inf
    cnt = 0
    while cnt < 40:
        k = int(rng.integers(5, 8))
        H = nx.gnm_random_graph(k, int(rng.integers(k, 2 * k)), seed=int(rng.integers(1 << 30)))
        if not nx.is_connected(H):
            continue
        EH = [tuple(sorted(e)) for e in H.edges()]
        F = [i for i in range(len(EH)) if rng.random() < 0.5]
        Lmax = int(rng.integers(2, 4))
        ell = {i: int(rng.integers(1, Lmax + 1)) for i in F}
        nG = k + sum(l - 1 for l in ell.values())
        if nG > 18:
            continue
        EG, nxt = [], k
        for i, (u, v) in enumerate(EH):
            l = ell.get(i, 1)
            seq = [u] + list(range(nxt, nxt + l - 1)) + [v]
            nxt += l - 1
            EG += [tuple(sorted(p)) for p in zip(seq, seq[1:])]
        hH = edge_expansion(k, EH)
        nu = max(sum(1 for i in F if EH[i][0] in S and EH[i][1] in S) / len(S)
                 for r in range(1, k + 1) for S in map(set, itertools.combinations(range(k), r)))
        hG = edge_expansion(nG, EG)
        bound = hH / (1 + (max(ell.values(), default=1) - 1) * (nu + hH))
        worst = min(worst, hG / bound)
        cnt += 1
    print(f"P4 Lemma lem:subdivide on {cnt} random subdivisions: smallest h(G)/bound = {worst:.3f}")
    check(worst >= 1 - 1e-12, "P4 lem:subdivide")

    # Propositions prop:hub-expansion and prop:hub on the smallest hub graphs
    for a, b in ((2, 2), (2, 3)):
        n, L = 2, 2
        XA, XB = cubic_core(a * n, 0), cubic_core(b * n, 1)
        h0 = min(edge_expansion(a * n, list(XA.edges())), edge_expansion(b * n, list(XB.edges())))
        G, idx, E, inA, firsts = hub_graph(a, b, n, L, XA, XB)
        N = G.number_of_nodes()
        h = edge_expansion(N, E)
        h1, m = min(h0, 1), min(a, b)
        lo, hi = 2 * h1 * m / (5 * (a + b + 1) * L), 2 / (L - 1)
        Dmax = max(d for _, d in G.degree())
        print(f"P4 G_(L=2,n=2)({a},{b}): N = {N} (formula {n * (1 + (a + b) * L)}), max degree {Dmax}, "
              f"2-connected {nx.node_connectivity(G) >= 2}, h = {h:.4f} in [{lo:.4f}, {hi:.4f}]")
        check(N == n * (1 + (a + b) * L) and Dmax == max(a + b, 4) and nx.node_connectivity(G) >= 2
              and lo <= h + 1e-12 and h <= hi + 1e-12, f"P4 hub expansion ({a},{b})")
    for L in (2, 3, 5, 8):
        a = b = 2
        n = 2
        G, idx, E, inA, firsts = hub_graph(a, b, n, L, cubic_core(4, 0), cubic_core(4, 1))
        N = G.number_of_nodes()
        lt = logtau(N, E)
        for s, cnt in (("A", a), ("B", b)):
            out = {firsts[(0, s, p)] for p in range(cnt)}
            P = math.exp(logtau(N, [E[i] for i in range(len(E)) if i not in out]) - lt)
            print(f"P4 prop:hub, L={L}: Pr[all {s}-edges at a hub are out] = {P:.5f} <= L^-{cnt} = {L ** -cnt:.5f}")
            check(P <= L ** -cnt + 1e-12, f"P4 prop:hub L={L}")

    # Proposition prop:hub-pairs
    a = b = 2
    k, n = 2 * a - 1, 13
    XA, XB = cubic_core(a * n, 3), cubic_core(b * n, 4)
    lamA = sorted(nx.laplacian_spectrum(XA))[1]
    lamB = sorted(nx.laplacian_spectrum(XB))[1]
    L = max(3, math.ceil(4 * k * (1 / lamA + 1 / lamB)))
    G, idx, E, inA, firsts = hub_graph(a, b, n, L, XA, XB)
    N = G.number_of_nodes()
    pair = {g: (firsts[(g, "A", 0)], firsts[(g, "B", 0)]) for g in range(n)}
    Rmid = 2 / lamA + 2 / lamB + L / (n - 1) * (1 / (a - 1) + 1 / (b - 1))
    gen = (2 * L / k - Rmid) / ((L + L / k) * (2 * L + Rmid))
    ymin, qmax, Bmin, npast = math.inf, 0.0, math.inf, 0
    for w in range(4):
        e, f = pair[w]
        others = [g for g in range(n) if g != w]
        allp = [x for g in others for x in pair[g]]
        pasts = [([], []), ([], allp), (allp, [])]
        for _ in range(6):
            inc = [x for x in allp if rng.random() < 0.5]
            pasts.append((inc, [x for x in allp if x not in inc]))
        for inc, exc in pasts:
            Y, pe, pf = sparse_pair(N, E, e, f, idx[("H", w)], inc, exc)
            qe, qf = 1 - pe, 1 - pf
            ymin = min(ymin, Y)
            qmax = max(qmax, qe, qf)
            Bmin = min(Bmin, Y * Y * max(2, 1 / (2 * math.sqrt(qe * qf))))
            npast += 1
    proved = 1 / (405 * (a + b) ** 5 * L ** 3)
    print(f"P4 prop:hub-pairs on G_(L,13)(2,2): lambda_A = {lamA:.4f}, lambda_B = {lamB:.4f}, L = {L}, N = {N}; "
          f"over {npast} pasts, min Y = {ymin:.3e} >= 1/(3(a+b)L) = {1 / (3 * (a + b) * L):.3e} and >= the "
          f"bound before simplification {gen:.3e}; max q = {qmax:.3e} < 1/L = {1 / L:.3e}; "
          f"n B_min^2 / N = {n * Bmin ** 2 / N:.3e} >= 1/(405 (a+b)^5 L^3) = {proved:.3e}")
    check(ymin >= 1 / (3 * (a + b) * L) and ymin >= gen - 1e-12 and qmax < 1 / L and Rmid <= L / k
          and n * Bmin ** 2 / N >= proved, "P4 prop:hub-pairs")
    print(f"P4 the sufficient condition L >= 48k/h0^2 for a = b = 2, h0 = 0.1: L >= {48 * 3 / 0.01:.0f}")


def sparse_pair(n, E, e, f, w, incl, excl):
    """Y, Pr[e in T], Pr[f in T] in K_sigma, by sparse LU of the Laplacian grounded at w."""
    m, Es, cls = contract(n, E, set(incl), set(excl))
    W = cls[w]
    U = Es[e][0] if Es[e][1] == W else Es[e][1]
    V = Es[f][0] if Es[f][1] == W else Es[f][1]
    rows, cols, vals = [], [], []
    for x in Es:
        if x is None or x[0] == x[1]:
            continue
        p, q = x
        rows += [p, q, p, q]
        cols += [p, q, q, p]
        vals += [1.0, 1.0, -1.0, -1.0]
    Lg = sp.csc_matrix((vals, (rows, cols)), shape=(m, m))
    keep = np.array([i for i in range(m) if i != W])
    pos = -np.ones(m, int)
    pos[keep] = np.arange(len(keep))
    lu = spla.splu(Lg[keep][:, keep].tocsc())
    bu = np.zeros(m - 1)
    bu[pos[U]] = 1.0
    bv = np.zeros(m - 1)
    bv[pos[V]] = 1.0
    phi, psi = lu.solve(bu), lu.solve(bv)
    Y, Y2 = phi[pos[V]], psi[pos[U]]
    check(abs(Y - Y2) < 1e-8 * max(1.0, abs(Y)), "transfer current not symmetric")
    return Y, phi[pos[U]], psi[pos[V]]


# ------------------------------------------------------------------ P5
def p5_triangles():
    for name, H in (("K4", nx.complete_graph(4)), ("K33", nx.complete_bipartite_graph(3, 3))):
        H = nx.convert_node_labels_to_integers(H)
        k = H.number_of_nodes()
        E = [tuple(sorted(e)) for e in H.edges()]
        nA = len(E)
        for v in range(k):
            x, y = k + 2 * v, k + 2 * v + 1
            E += [(v, x), (v, y), (x, y)]
        n, M = 3 * k, len(E)
        T = enum_trees(n, E)
        Am = np.uint64(sum(1 << i for i in range(nA)))
        beta = opnorm_sqrt_joint(T & Am, T & ~Am)
        covs = [abs(exact_pair_law(T, i, j)[1, 1] - exact_pair_law(T, i, j)[1].sum() * exact_pair_law(T, i, j)[:, 1].sum())
                for i in range(nA) for j in range(nA, M) if set(E[i]) & set(E[j])]
        hH = edge_expansion(k, E[:nA])
        hG = edge_expansion(n, E)
        Dmax = max(sum(1 for e in E if v in e) for v in range(n))
        print(f"P5 triangles attached to {name}: N = {n}, M = {M}, |A| = {nA}, {len(T)} trees, beta = {beta:.12f}, "
              f"max |Cov| over adjacent A-B pairs = {max(covs):.1e}, h(G) = {hG:.4f} >= min(h(H)/3, 1) = "
              f"{min(hH / 3, 1):.4f}, max degree {Dmax}")
        check(abs(beta - 1) < 1e-9 and max(covs) < 1e-12 and hG >= min(hH / 3, 1) - 1e-12 and Dmax == 5
              and 3 * nA >= M and 3 * (M - nA) <= 2 * M, f"P5 {name}")


# ------------------------------------------------------------------ P6
def wilson(n, adj, rng):
    in_tree = np.zeros(n, bool)
    in_tree[0] = True
    nxt = np.full(n, -1)
    for s in rng.permutation(n):
        u = int(s)
        while not in_tree[u]:
            nxt[u] = adj[u][int(rng.integers(len(adj[u])))]
            u = nxt[u]
        u = int(s)
        while not in_tree[u]:
            in_tree[u] = True
            u = nxt[u]
    return {tuple(sorted((v, int(nxt[v])))) for v in range(1, n)}


def p6_regular(rng):
    worst_ratio = math.inf
    for d in (3, 4):
        for nv in (200, 1000):
            G = nx.random_regular_graph(d, nv, seed=int(rng.integers(1 << 30)))
            while not nx.is_connected(G):
                G = nx.random_regular_graph(d, nv, seed=int(rng.integers(1 << 30)))
            E = [tuple(sorted(e)) for e in G.edges()]
            m = len(E)
            adj = [list(G.neighbors(v)) for v in range(nv)]
            at = [[] for _ in range(nv)]
            for i, (u, v) in enumerate(E):
                at[u].append(i)
                at[v].append(i)
            for split in ("random", "bfs"):
                if split == "random":
                    A = set(rng.permutation(m)[: m // 2].tolist())
                else:
                    dist = nx.single_source_shortest_path_length(G, 0)
                    order = sorted(range(m), key=lambda i: (min(dist[E[i][0]], dist[E[i][1]]), rng.random()))
                    A = set(order[: m // 2])
                Z = [v for v in range(nv) if any(i in A for i in at[v]) and any(i not in A for i in at[v])]
                W, blocked = [], set()
                for z in rng.permutation(Z):
                    z = int(z)
                    if z in blocked:
                        continue
                    W.append(z)
                    blocked |= set(nx.single_source_shortest_path_length(G, z, cutoff=3))
                    if len(W) >= 30:
                        break
                pairs = []
                for w in W:
                    best = None
                    for i in at[w]:
                        for j in at[w]:
                            if i in A and j not in A:
                                Y = sparse_pair(nv, E, i, j, w, [], [])[0]
                                if best is None or Y > best[0]:
                                    best = (Y, i, j)
                    pairs.append((w, best[1], best[2]))
                trees = [wilson(nv, adj, rng) for _ in range(3)]
                ymin = math.inf
                for k, (w, i, j) in enumerate(pairs):
                    past = [x for (_, a, b) in pairs[:k] for x in (a, b)]
                    pats = [([], past), (past, [])]
                    for _ in range(6):
                        inc = [x for x in past if rng.random() < 0.5]
                        pats.append((inc, [x for x in past if x not in inc]))
                    for T in trees:
                        inc = [x for x in past if E[x] in T]
                        pats.append((inc, [x for x in past if x not in inc]))
                    for incl, excl in pats:
                        if not feasible(nv, E, incl, excl):
                            continue
                        ymin = min(ymin, sparse_pair(nv, E, i, j, w, incl, excl)[0])
                tree = (d - 2) / (d * (d - 1))
                worst_ratio = min(worst_ratio, ymin / tree)
                print(f"P6 D={d} N={nv} {split:6s} split: {len(pairs)} spread-out mixed vertices, "
                      f"min Y_sigma = {ymin:.4f}, ratio to (D-2)/(D(D-1)) = {ymin / tree:.3f}")
    print(f"P6 smallest ratio {worst_ratio:.3f} (numerical evidence only)")
    check(worst_ratio >= 0.81, "P6 smallest ratio quoted in the paper (0.81)")


def feasible(n, E, incl, excl):
    uf = nx.utils.UnionFind(range(n))
    for i in incl:
        if uf[E[i][0]] == uf[E[i][1]]:
            return False
        uf.union(*E[i])
    H = nx.MultiGraph()
    H.add_nodes_from({uf[x] for x in range(n)})
    ex = set(excl) | set(incl)
    H.add_edges_from((uf[E[i][0]], uf[E[i][1]]) for i in range(len(E)) if i not in ex)
    return nx.is_connected(H)


def main():
    rng = np.random.default_rng(20261104)
    for step in (p1_two_by_two, p2_edge_pair_law, p3_chain, p4_hub):
        step(rng)
    p5_triangles()
    p6_regular(rng)
    if FAILS:
        print("FAILED:", FAILS)
        sys.exit(1)
    print("check_edge_pairs: all checks passed")


if __name__ == "__main__":
    main()

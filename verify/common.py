"""Shared helpers for the vertex-touch checks (section "Vertex-touch bits and 2-connected expanders").

Ported from toy/r2/common.py. Graphs, spanning-tree enumeration (C++ helper enum_trees.cpp),
effective resistances, and the bounds pA(w), pB(w) of Lemma "touch probabilities".

Effective resistances are computed as (e_x - e_y)^T (L + J/k)^{-1} (e_x - e_y) on a connected
vertex set of size k, which equals (e_x - e_y)^T L^+ (e_x - e_y). We do not use np.linalg.pinv:
on residual graphs of Tutte-Coxeter it returned resistances off by more than 1 (see the review
notes/approaches/R2-review.md, Sec. 6).
"""
import os, subprocess, tempfile
import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))

LCF = {
    "mobius_kantor": (16, [5, -5], 8),
    "mcgee": (24, [12, 7, -7], 8),
    "tutte_coxeter": (30, [-13, -9, 7, -7, 9, 13], 5),                                         # girth 8
    "foster": (90, [17, -9, 37, -37, 9, -17], 15),                                              # girth 10
    "tutte12": (126, [17, 27, -13, -59, -35, 35, -11, 13, -53, 53, -27, 21, 57, 11, -21, -57, 59, -17], 7),  # girth 12
}


def named_graph(name, seed=0):
    """Simple 2-connected test graphs. rr<d>_<n> is nx.random_regular_graph(d, n), resampled
    (seed += 1000) until it is 2-connected; the rejection means these graphs are not evidence
    for any probabilistic statement about random regular graphs."""
    if name == "petersen":
        G = nx.petersen_graph()
    elif name == "heawood":
        G = nx.heawood_graph()
    elif name == "pappus":
        G = nx.pappus_graph()
    elif name in LCF:
        n, s, r = LCF[name]
        G = nx.LCF_graph(n, s, r)
    elif name == "K33":
        G = nx.complete_bipartite_graph(3, 3)
    elif name == "cube":
        G = nx.hypercube_graph(3)
    elif name.startswith("rr"):
        d, n = name[2:].split("_")
        G = nx.random_regular_graph(int(d), int(n), seed=seed)
        while not nx.is_connected(G) or nx.node_connectivity(G) < 2:
            seed += 1000
            G = nx.random_regular_graph(int(d), int(n), seed=seed)
    else:
        raise ValueError(name)
    G = nx.convert_node_labels_to_integers(G)
    assert nx.is_connected(G) and nx.node_connectivity(G) >= 2, f"{name} is not 2-connected"
    return G


def edge_list(G):
    return [tuple(sorted(e)) for e in G.edges()]


def enum_trees(n, E):
    """All spanning trees, as uint64 edge bitmasks (needs len(E) <= 64)."""
    assert len(E) <= 64
    exe = os.path.join(HERE, "enum_trees")
    if not os.path.exists(exe):
        subprocess.run(["g++", "-O2", "-o", exe, os.path.join(HERE, "enum_trees.cpp")], check=True)
    inp = f"{n} {len(E)}\n" + "".join(f"{u} {v}\n" for u, v in E)
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
        path = f.name
    try:
        subprocess.run([exe, path], input=inp.encode(), check=True, capture_output=True)
        T = np.fromfile(path, dtype=np.uint64)
    finally:
        os.unlink(path)
    return T


# ---------------------------------------------------------------- electrical
def laplacian(n, E, keep_vertices=None, drop_edges=()):
    drop = set(drop_edges)
    L = np.zeros((n, n))
    for i, (u, v) in enumerate(E):
        if i in drop:
            continue
        if keep_vertices is not None and (u not in keep_vertices or v not in keep_vertices):
            continue
        L[u, u] += 1; L[v, v] += 1; L[u, v] -= 1; L[v, u] -= 1
    return L


def resistance_matrix(L, verts):
    """Effective resistances among `verts` in the graph with Laplacian L restricted to verts.
    Returns (R, connected); R is None if the graph on verts is disconnected."""
    verts = list(verts)
    Ls = L[np.ix_(verts, verts)]
    w = np.linalg.eigvalsh(Ls)
    if len(verts) < 2 or w[1] < 1e-9:
        return None, False
    k = len(verts)
    P = np.linalg.inv(Ls + np.ones((k, k)) / k)
    d = np.diag(P)
    R = d[:, None] + d[None, :] - 2 * P
    assert np.all(np.isfinite(R))
    return R, True


def algebraic_connectivity(L, verts):
    Ls = L[np.ix_(list(verts), list(verts))]
    return float(np.linalg.eigvalsh(Ls)[1])


# ---------------------------------------------------------------- splits / mixed vertices
def star(n, E):
    S = [[] for _ in range(n)]
    for i, (u, v) in enumerate(E):
        S[u].append(i); S[v].append(i)
    return S


def mixed_vertices(n, E, inA):
    S = star(n, E)
    return [v for v in range(n) if any(inA[i] for i in S[v]) and any(not inA[i] for i in S[v])]


def greedy_separated(G, cand, sep, rng=None):
    """Greedy subset of cand with pairwise graph distance >= sep."""
    cand = list(cand)
    if rng is not None:
        rng.shuffle(cand)
    chosen, blocked = [], set()
    for v in cand:
        if v in blocked:
            continue
        chosen.append(v)
        blocked |= set(nx.single_source_shortest_path_length(G, v, cutoff=sep - 1))
    return chosen


def touch_bounds(n, E, W, inA):
    """For each w in W (an independent set), the bounds of Lemma "touch probabilities":
       pA(w) = prod_{u in N_A(w)} 1/(1+R_u),  R_u = R_{G - (W minus w) - S_A(w)}(u, w),
    computed on the component of w (R = inf, factor 0, if u is not connected to w), the same
    for B, and the coarse bound (2+rho_w)^{-|S_A(w)|}, rho_w = max_{x,y in N(w)} R_{G-W}(x,y)
    (rho_w = inf if some pair is disconnected in G - W)."""
    S = star(n, E)
    Wset = set(W)
    keepH = [x for x in range(n) if x not in Wset]
    LH = laplacian(n, E, keep_vertices=set(keepH))
    Hg = nx.Graph(); Hg.add_nodes_from(keepH); Hg.add_edges_from(e for e in E if e[0] not in Wset and e[1] not in Wset)
    comp_of = {}
    for ci, comp in enumerate(nx.connected_components(Hg)):
        for x in comp: comp_of[x] = ci
    resH = {}
    out = []
    for w in W:
        nb = [E[i][0] if E[i][1] == w else E[i][1] for i in S[w]]
        rho = 0.0
        for x in nb:
            for y in nb:
                if x == y: continue
                if comp_of[x] != comp_of[y]:
                    rho = np.inf; continue
                ci = comp_of[x]
                if ci not in resH:
                    verts = sorted(v for v in keepH if comp_of[v] == ci)
                    R, ok = resistance_matrix(LH, verts)
                    resH[ci] = (R, {v: k for k, v in enumerate(verts)})
                R, idx = resH[ci]
                rho = max(rho, float(R[idx[x], idx[y]]))
        res = {"w": w, "rho": rho}
        keep = set(keepH) | {w}
        for side in ("A", "B"):
            rem = [i for i in S[w] if inA[i] == (side == "A")]
            Kg = nx.Graph(); Kg.add_nodes_from(keep)
            Kg.add_edges_from(E[i] for i in range(len(E)) if i not in set(rem) and E[i][0] in keep and E[i][1] in keep)
            compw = sorted(nx.node_connected_component(Kg, w))
            p = 1.0
            if len(compw) >= 2:
                L = laplacian(n, E, keep_vertices=set(compw), drop_edges=rem)
                R, ok = resistance_matrix(L, compw)
                assert ok
                idx = {x: k for k, x in enumerate(compw)}
            for i in rem:
                u = E[i][0] if E[i][1] == w else E[i][1]
                p *= (1.0 / (1.0 + R[idx[u], idx[w]])) if u in set(compw) else 0.0
            res[side] = p
            res[side + "_coarse"] = (2 + rho) ** (-len(rem)) if np.isfinite(rho) else 0.0
        out.append(res)
    return out

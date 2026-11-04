# /// script
# dependencies = ["numpy", "networkx"]
# ///
"""Exhaustive checks of the block results and of the electrical lemmas on small graphs.

1. Proposition "the block criterion": for every split A|B of several small graphs, beta = 1 exactly
   when every block lies on one side; a balanced split with beta = 1 exists exactly when no block
   has more than 2M/3 edges. K4 plus a pendant edge: max beta over balanced splits is 0.75.
2. Theorem "block dichotomy" (b), lower half: max over balanced splits of beta^2 >= 2^-(M1-2M/3)-1.
3. Lemma "blocks inherit expansion": h(G1) >= h(G) for every block with at least two edges
   (exact h by enumerating all vertex subsets).
4. One-vertex deletions, the case that uses 2-connectivity. These are the deletion and resistance
   lemmas of the previous version of the paper, which the current version replaces (the current
   lemmas are checked in check_deletion.py): h(G - w) >= h(G)/(2D) for every vertex w, and
   R_{G-w} <= 16 D^3 / h(G)^2, together with h(H)^2 <= 2 D lambda_2(L_H) for H = G and H = G - w,
   which follows from Lemma "Cheeger-Mohar lower bound".
5. Lemma "deleting one edge": Pr[e not in T] = 1/(1 + R_{G-e}(u,w)) for every edge, by enumeration.
Exits with status 1 on any failure.
"""
import sys
import numpy as np
import networkx as nx
from common import enum_trees, edge_list, laplacian, resistance_matrix, algebraic_connectivity, LCF

FAIL = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        FAIL.append(msg)


def beta_of(T, maskA, tau):
    a = T & np.uint64(maskA)
    b = T & ~np.uint64(maskA)
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    M = np.zeros((len(ua), len(ub)))
    np.add.at(M, (ia, ib), 1.0 / np.sqrt(tau))
    return float(np.linalg.norm(M, 2))


def edge_expansion(n, E):
    """Exact h = min |boundary S|/|S| over 0 < |S| <= n/2, by enumerating subsets (n <= 22)."""
    assert n <= 22
    best = np.inf
    masks = np.arange(1, 1 << n, dtype=np.uint32)
    size = np.zeros(len(masks), dtype=np.int32)
    for v in range(n):
        size += ((masks >> np.uint32(v)) & np.uint32(1)).astype(np.int32)
    keep = (2 * size <= n)
    masks, size = masks[keep], size[keep]
    bd = np.zeros(len(masks), dtype=np.int32)
    for u, v in E:
        bd += (((masks >> np.uint32(u)) ^ (masks >> np.uint32(v))) & np.uint32(1)).astype(np.int32)
    return float(np.min(bd / size))


def blocks_of(G, E):
    idx = {e: i for i, e in enumerate(E)}
    return [sorted(idx[tuple(sorted(e))] for e in comp) for comp in nx.biconnected_component_edges(G)]


def graph(edges):
    G = nx.Graph(); G.add_edges_from(edges)
    return nx.convert_node_labels_to_integers(G, ordering="sorted")


K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
SMALL = {
    "K4+leaf": graph(K4 + [(3, 4)]),
    "bowtie": graph([(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (2, 4)]),
    "K4+path2": graph(K4 + [(3, 4), (4, 5)]),
    "C4+triangle": graph([(0, 1), (1, 2), (2, 3), (0, 3), (3, 4), (4, 5), (3, 5)]),
    "prism_K3": nx.circular_ladder_graph(3),
    "two_K4": graph(K4 + [(3, 4), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)]),
}

# ---------------------------------------------------------------- 1, 2: block criterion
for name, G in SMALL.items():
    n = G.number_of_nodes(); E = edge_list(G); m = len(E)
    T = enum_trees(n, E); tau = len(T)
    blocks = blocks_of(G, E)
    M1 = max(len(b) for b in blocks)
    crit_ok, best_bal, argbest = True, 0.0, []
    for mask in range(1, (1 << m) - 1):
        beta = beta_of(T, mask, tau)
        whole = all(len({(mask >> i) & 1 for i in b}) == 1 for b in blocks)
        if (abs(beta - 1) < 1e-9) != whole:
            crit_ok = False
        k = bin(mask).count("1")
        if 3 * k >= m and 3 * (m - k) >= m:
            if beta > best_bal + 1e-12:
                best_bal, argbest = beta, [mask]
            elif abs(beta - best_bal) <= 1e-12:
                argbest.append(mask)
    check(crit_ok, f"{name}: M={m}, blocks {sorted(len(b) for b in blocks)}: beta=1 iff every block on one side (all {2**m - 2} splits)")
    check((abs(best_bal - 1) < 1e-9) == (3 * M1 <= 2 * m),
          f"{name}: max balanced beta = {best_bal:.6f}; balanced beta=1 exists iff M1 <= 2M/3 (M1={M1})")
    if 3 * M1 > 2 * m:
        lb = 2.0 ** (-(M1 - 2 * m / 3) - 1)
        check(best_bal ** 2 >= lb - 1e-12, f"{name}: max balanced beta^2 = {best_bal**2:.6f} >= 2^-(M1-2M/3)-1 = {lb:.6f}")
    if name == "K4+leaf":
        check(abs(best_bal - 0.75) < 1e-12, "K4+leaf: max balanced beta = 0.75")
        tri_mask = sum(1 << E.index(e) for e in [(0, 1), (0, 2), (1, 2)])
        check(abs(beta_of(T, tri_mask, tau) - 0.75) < 1e-12,
              f"K4+leaf: A = a triangle of K4 attains it ({len(argbest)} of the balanced splits attain 0.75)")

# ---------------------------------------------------------------- 3: blocks inherit expansion
pet = nx.petersen_graph()
pend = nx.Graph(pet); pend.add_edges_from((v, v + 10) for v in range(10))
EXP = dict(SMALL)
EXP["petersen+pendants"] = pend
EXP["petersen+K4"] = graph(list(pet.edges()) + [(0, 10), (0, 11), (0, 12), (10, 11), (10, 12), (11, 12)])
for name, G in EXP.items():
    n = G.number_of_nodes(); E = edge_list(G)
    hG = edge_expansion(n, E)
    ok = True
    for b in blocks_of(G, E):
        if len(b) < 2:
            continue
        H = nx.convert_node_labels_to_integers(G.edge_subgraph([E[i] for i in b]).copy())
        hB = edge_expansion(H.number_of_nodes(), edge_list(H))
        ok &= hB >= hG - 1e-12
    check(ok, f"{name}: h(G) = {hG:.4f} <= h(block) for every block with >= 2 edges")

# ---------------------------------------------------------------- 4: deleting one vertex, resistances, Mohar
def cayley_like(name):
    if name in LCF:
        n, s, r = LCF[name]
        return nx.convert_node_labels_to_integers(nx.LCF_graph(n, s, r))
    if name.startswith("rr"):
        d, n = map(int, name[2:].split("_")); seed = 1
        G = nx.random_regular_graph(d, n, seed=seed)
        while not nx.is_connected(G) or nx.node_connectivity(G) < 2:
            seed += 1000; G = nx.random_regular_graph(d, n, seed=seed)
        return G
    return {"petersen": nx.petersen_graph, "heawood": nx.heawood_graph, "pappus": nx.pappus_graph,
            "K33": lambda: nx.complete_bipartite_graph(3, 3), "cube": lambda: nx.hypercube_graph(3),
            "K5": lambda: nx.complete_graph(5), "prism_K4": lambda: nx.cartesian_product(nx.complete_graph(4), nx.path_graph(2)),
            "C10_square": lambda: nx.power(nx.cycle_graph(10), 2)}[name]()

for name in ["K33", "cube", "petersen", "heawood", "mobius_kantor", "pappus", "K5", "prism_K4", "C10_square", "rr3_14", "rr3_16", "rr4_12"]:
    G = nx.convert_node_labels_to_integers(cayley_like(name))
    assert nx.node_connectivity(G) >= 2
    n = G.number_of_nodes(); E = edge_list(G); D = max(d for _, d in G.degree())
    hG = edge_expansion(n, E)
    lam = algebraic_connectivity(laplacian(n, E), range(n))
    ok_mohar = hG ** 2 <= 2 * D * lam + 1e-9
    min_ratio, max_R, ok = np.inf, 0.0, True
    for w in range(n):
        keep = [x for x in range(n) if x != w]
        H = nx.convert_node_labels_to_integers(G.subgraph(keep).copy())
        EH = edge_list(H)
        hH = edge_expansion(n - 1, EH)
        min_ratio = min(min_ratio, hH / (hG / (2 * D)))
        LH = laplacian(n - 1, EH)
        lamH = algebraic_connectivity(LH, range(n - 1))
        R, conn = resistance_matrix(LH, range(n - 1))
        ok &= conn and hH >= hG / (2 * D) - 1e-12 and hH ** 2 <= 2 * D * lamH + 1e-9
        max_R = max(max_R, float(R.max()))
    ok &= max_R <= 16 * D ** 3 / hG ** 2
    check(ok and ok_mohar, f"{name}: N={n} D={D} h={hG:.4f}: min_w h(G-w)/(h/2D) = {min_ratio:.2f} >= 1, "
          f"max R_(G-w) = {max_R:.3f} <= 16D^3/h^2 = {16 * D**3 / hG**2:.0f}, Mohar holds")

# ---------------------------------------------------------------- 5: deleting one edge
for name in ["K4+leaf", "bowtie", "prism_K3", "two_K4"]:
    G = SMALL[name]; n = G.number_of_nodes(); E = edge_list(G)
    T = enum_trees(n, E); tau = len(T)
    err = 0.0
    for i, (u, v) in enumerate(E):
        p_out = float(np.count_nonzero((T >> np.uint64(i)) & np.uint64(1) == 0)) / tau
        L = laplacian(n, E, drop_edges=[i])
        R, conn = resistance_matrix(L, range(n))
        pred = 1.0 / (1.0 + R[u, v]) if conn else 0.0
        err = max(err, abs(p_out - pred))
    check(err < 1e-12, f"{name}: Pr[e not in T] = 1/(1+R_(G-e)) for every edge (max error {err:.1e})")

print(f"{'ALL PASS' if not FAIL else str(len(FAIL)) + ' FAILURES'}")
sys.exit(1 if FAIL else 0)

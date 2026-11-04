# /// script
# dependencies = ["numpy", "networkx"]
# ///
"""Checks of the deletion lemmas behind Theorem "2-connected expanders" (section "Expansion").

1. Lemma "Cheeger-Mohar lower bound": lambda_2(L_H) >= mu(Delta, i) = Delta - sqrt(Delta^2 - i^2)
   >= i^2/(2 Delta) for every connected graph H in a list (K2, K3, K4, paths, cycles, the named
   graphs and all their one-vertex deletions), with i = h(H) computed exactly. Also prints the
   K3 example where the upper form i <= sqrt(lambda_2 (2 Delta - lambda_2)) fails (lambda_2 > Delta).
2. Lemma "counting connected sets": for every connected S containing vertex 0, with |S| = k and
   outer vertex boundary |N(S)| = j, the number of such S is at most binom(k-1+j, j) (exhaustive).
3. Lemma "deleting a 3-separated set": h(G-W) >= h(G) - 1 and lambda_2(G-W) >= lambda_2(G) - 1 for
   random greedy 3-separated sets W (exact h by enumerating subsets, N <= 22), and the spectral part
   on random 6- and 8-regular graphs with N = 60.
4. Lemma "random sparse deletion", deterministic part: for random R inside a greedy independent set,
   R* = R minus every w in a dangerous pair (w, S) satisfies h(G - R*) >= h' = h/(D+1), where h is
   the exact edge expansion of G (all connected S with |S| <= N/2 enumerated, N <= 16). On graphs this
   small every dangerous threshold m(b) is 2, so R* has at most one vertex and the check exercises only
   the one-neighbour case of the proof (the case that uses 2-connectivity).
Exits with status 1 on any failure.
"""
import math, sys
import numpy as np
import networkx as nx
from common import named_graph, edge_list, laplacian, algebraic_connectivity

FAIL = []
TOL = 1e-9


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        FAIL.append(msg)


def edge_expansion(n, E):
    """Exact h = min |boundary S|/|S| over 0 < |S| <= n/2 (n <= 22); inf if n < 2."""
    assert n <= 22
    if n < 2:
        return math.inf
    masks = np.arange(1, 1 << n, dtype=np.uint32)
    size = np.zeros(len(masks), dtype=np.int32)
    for v in range(n):
        size += ((masks >> np.uint32(v)) & np.uint32(1)).astype(np.int32)
    keep = 2 * size <= n
    masks, size = masks[keep], size[keep]
    bd = np.zeros(len(masks), dtype=np.int32)
    for u, v in E:
        bd += (((masks >> np.uint32(u)) ^ (masks >> np.uint32(v))) & np.uint32(1)).astype(np.int32)
    return float(np.min(bd / size))


def relabel(H):
    return nx.convert_node_labels_to_integers(H)


def lam2(G):
    n = G.number_of_nodes()
    return algebraic_connectivity(laplacian(n, edge_list(G)), range(n))


def mu(Delta, i):
    return Delta - math.sqrt(max(0.0, Delta * Delta - i * i))


def nbmasks(G):
    n = G.number_of_nodes()
    nb = [0] * n
    for u, v in G.edges():
        nb[u] |= 1 << v; nb[v] |= 1 << u
    return nb


def is_conn(S, root, nb):
    seen = frontier = 1 << root
    while frontier:
        new = 0; f = frontier
        while f:
            x = (f & -f).bit_length() - 1; f &= f - 1
            new |= nb[x]
        new &= S & ~seen
        seen |= new; frontier = new
    return seen == S


def outer(S, nb):
    out = 0; f = S
    while f:
        x = (f & -f).bit_length() - 1; f &= f - 1
        out |= nb[x]
    return out & ~S


def greedy_sep(G, cand, sep, rng):
    cand = list(cand); rng.shuffle(cand)
    chosen, blocked = [], set()
    for v in cand:
        if v in blocked:
            continue
        chosen.append(v)
        blocked |= set(nx.single_source_shortest_path_length(G, v, cutoff=sep - 1))
    return chosen


NAMES = ["K33", "cube", "petersen", "heawood", "mobius_kantor", "pappus", "rr3_14", "rr3_16", "rr4_12"]
rng = np.random.default_rng(3)

# ---------------------------------------------------------------- 1: Cheeger-Mohar lower bound
graphs = {"K2": nx.complete_graph(2), "K3": nx.complete_graph(3), "K4": nx.complete_graph(4), "K5": nx.complete_graph(5),
          "P5": nx.path_graph(5), "C7": nx.cycle_graph(7), "K4-e": nx.complete_graph(4)}
graphs["K4-e"].remove_edge(0, 1)
for name in NAMES:
    graphs[name] = named_graph(name, seed=1)
worst, count = math.inf, 0
for name, G in list(graphs.items()):
    family = [G] + ([relabel(G.subgraph([x for x in G if x != v]).copy()) for v in G] if G.number_of_nodes() > 3 else [])
    ok = True
    for H in family:
        if H.number_of_nodes() < 2 or not nx.is_connected(H):
            continue
        n = H.number_of_nodes(); Delta = max(d for _, d in H.degree())
        i = edge_expansion(n, edge_list(H)); l2 = lam2(H)
        ok &= i <= Delta + TOL and l2 >= mu(Delta, i) - TOL and mu(Delta, i) >= i * i / (2 * Delta) - TOL
        worst = min(worst, l2 / mu(Delta, i)); count += 1
    check(ok, f"{name} and its one-vertex deletions: lambda_2 >= Delta - sqrt(Delta^2 - h^2) >= h^2/(2 Delta)")
print(f"      ({count} graphs, min lambda_2/mu = {worst:.4f})")
K3 = graphs["K3"]
check(abs(lam2(K3) - 3) < TOL and edge_expansion(3, edge_list(K3)) == 2.0 and 2.0 > math.sqrt(3 * (4 - 3)),
      "K3: lambda_2 = 3 > Delta = 2 and h = 2 > sqrt(lambda_2 (2 Delta - lambda_2)) = sqrt(3): the upper form needs lambda_2 <= Delta")

# ---------------------------------------------------------------- 2: counting connected sets
for name in ["K33", "cube", "petersen", "heawood", "rr4_12", "mobius_kantor", "rr3_16"]:
    G = graphs[name]; n = G.number_of_nodes(); nb = nbmasks(G)
    cnt = {}
    for rest in range(1 << (n - 1)):
        S = 1 | (rest << 1)
        if not is_conn(S, 0, nb):
            continue
        key = (bin(S).count("1"), bin(outer(S, nb)).count("1"))
        cnt[key] = cnt.get(key, 0) + 1
    ratio = max(c / math.comb(k - 1 + j, j) for (k, j), c in cnt.items())
    check(ratio <= 1, f"{name}: N={n}, {sum(cnt.values())} connected sets containing vertex 0, max count/binom(k-1+j, j) = {ratio:.4f}")

# ---------------------------------------------------------------- 3: deleting a 3-separated set
for name in NAMES:
    G = graphs[name]; n = G.number_of_nodes(); E = edge_list(G)
    hG, lG = edge_expansion(n, E), lam2(G)
    worst_h, worst_l, ok = math.inf, math.inf, True
    for t in range(6):
        W = greedy_sep(G, range(n), 3, rng)
        H = relabel(G.subgraph([x for x in range(n) if x not in set(W)]).copy())
        nH = H.number_of_nodes()
        hH = edge_expansion(nH, edge_list(H)); lH = lam2(H)
        ok &= hH >= hG - 1 - TOL and lH >= lG - 1 - TOL
        worst_h, worst_l = min(worst_h, hH - (hG - 1)), min(worst_l, lH - (lG - 1))
    check(ok, f"{name}: h={hG:.3f}, lambda_2={lG:.3f}; 6 greedy 3-separated W: min h(G-W)-(h-1) = {worst_h:.3f}, "
          f"min lambda_2(G-W)-(lambda_2-1) = {worst_l:.3f}")
for d, n in [(6, 60), (8, 60)]:
    G = nx.random_regular_graph(d, n, seed=5); assert nx.is_connected(G)
    lG = lam2(G); ok = True; worst = math.inf
    for t in range(10):
        W = greedy_sep(G, range(n), 3, rng)
        H = relabel(G.subgraph([x for x in range(n) if x not in set(W)]).copy())
        lH = lam2(H); ok &= lH >= lG - 1 - TOL; worst = min(worst, lG - lH)
    check(ok, f"random {d}-regular N={n}: lambda_2 = {lG:.3f}, 10 greedy 3-separated W: largest drop of lambda_2 = {worst:.3f} <= 1")

# ---------------------------------------------------------------- 4: random sparse deletion, deterministic part
for name in ["K33", "cube", "petersen", "heawood", "rr4_12", "rr3_14", "mobius_kantor", "rr3_16"]:
    G = graphs[name]; n = G.number_of_nodes(); E = edge_list(G); nb = nbmasks(G)
    D = max(d for _, d in G.degree()); h = edge_expansion(n, E); hp = h / (D + 1); theta = (1 - hp / h) / D
    conn = []
    for S in range(1, 1 << n):
        k = bin(S).count("1")
        if 2 * k > n or not is_conn(S, (S & -S).bit_length() - 1, nb):
            continue
        b = sum(((S >> u) ^ (S >> v)) & 1 for u, v in E)
        conn.append((S, outer(S, nb), max(2, math.floor(theta * b) + 1)))
    ok, sizes = True, []
    for q in (0.3, 0.6, 1.0):
        for t in range(10):
            W0 = greedy_sep(G, range(n), 2, rng)
            R = [w for w in W0 if rng.random() < q]
            Rm = sum(1 << w for w in R)
            bad = set()
            for S, NS, m in conn:
                if bin(NS & Rm).count("1") >= m:
                    bad |= {w for w in R if (NS >> w) & 1}
            Rs = [w for w in R if w not in bad]
            H = relabel(G.subgraph([x for x in range(n) if x not in set(Rs)]).copy())
            ok &= nx.is_connected(H) and edge_expansion(H.number_of_nodes(), edge_list(H)) >= hp - TOL
            sizes.append(len(Rs))
    check(ok, f"{name}: N={n}, h={h:.3f}, h'={hp:.3f}, {len(conn)} connected sets with |S|<=N/2; 30 random R: "
          f"h(G-R*) >= h' (|R*| from {min(sizes)} to {max(sizes)})")

print("ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES")
sys.exit(1 if FAIL else 0)

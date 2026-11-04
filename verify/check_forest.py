# /// script
# dependencies = []
# ///
"""Exact checks of the local completion bounds (Lemma 4.2 for trees, Proposition 7.4 for forests).

For a small graph G, a cycle C, and every outside configuration F' (the part of the random subgraph
outside I(C), the edges with an endpoint on C), check

  trees:   Pr[T  ∩ E(C) = E(C) - e | T  \\ I(C) = F'] >= 1 / min_{x0} prod_{y != x0} deg(y)
  forests: Pr[Fo ∩ E(C) = E(C) - e | Fo \\ I(C) = F'] >= 1 / min_{x0} prod_{y != x0} (deg(y) + 1)

where T is a uniform spanning tree and Fo a uniform spanning forest (uniform over all acyclic edge
subsets). All counts are exact integers; the comparison is exact (Fractions). It also confirms the
referee's counterexample: on the triangle the tree bound 1/4 fails for forests (actual value 1/7).
Exits nonzero on any violation. Standard library only.
"""
import itertools, sys
from collections import Counter, defaultdict
from fractions import Fraction


def acyclic_subsets(n, E):
    """All acyclic edge subsets (as bitmasks), and the spanning trees among them."""
    forests, trees = [], []
    m = len(E)
    for mask in range(1 << m):
        par = list(range(n))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        ok = True
        cnt = 0
        for i in range(m):
            if (mask >> i) & 1:
                a, b = find(E[i][0]), find(E[i][1])
                if a == b:
                    ok = False
                    break
                par[a] = b
                cnt += 1
        if ok:
            forests.append(mask)
            if cnt == n - 1:
                trees.append(mask)
    return forests, trees


def short_cycles(n, E, lmax=4):
    eid = {}
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(E):
        eid[(u, v)] = eid[(v, u)] = i
        adj[u].append(v)
        adj[v].append(u)
    out = {}
    def dfs(start, path):
        for w in adj[path[-1]]:
            if w == start and len(path) >= 3:
                es = sum(1 << eid[(path[i], path[(i + 1) % len(path)])] for i in range(len(path)))
                out[es] = frozenset(path)
            elif w > start and w not in path and len(path) < lmax:
                dfs(start, path + [w])
    for s in range(n):
        dfs(s, [s])
    return out


def prism(Hedges, nH):
    return 2 * nH, [(u, v) for u, v in Hedges] + [(u + nH, v + nH) for u, v in Hedges] + [(v, v + nH) for v in range(nH)]


def square_graph(Hedges, nH):
    adj = [set() for _ in range(nH)]
    for u, v in Hedges:
        adj[u].add(v); adj[v].add(u)
    E = set()
    for v in range(nH):
        for u in adj[v]:
            E.add(tuple(sorted((u, v))))
        for u, w in itertools.combinations(sorted(adj[v]), 2):
            E.add((u, w))
    return nH, sorted(E)


cyc = lambda k: [(i, (i + 1) % k) for i in range(k)]
K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
GRAPHS = {
    "C3": (3, cyc(3)),
    "K4": (4, K4),
    "prism_C3": prism(cyc(3), 3),
    "prism_C4 (cube)": prism(cyc(4), 4),
    "prism_C5": prism(cyc(5), 5),
    "prism_K4": prism(K4, 4),
    "C7_square": square_graph(cyc(7), 7),
    "two_triangles_bridge": (6, cyc(3) + [(3, 4), (4, 5), (3, 5), (2, 3)]),
}


def check(name, n, E):
    deg = [0] * n
    for u, v in E:
        deg[u] += 1; deg[v] += 1
    forests, trees = acyclic_subsets(n, E)
    cycles = short_cycles(n, E)
    res = {}
    for label, family, shift in (("trees", trees, 0), ("forests", forests, 1)):
        margin = None
        conds = 0
        for c, vs in cycles.items():
            inc = sum(1 << i for i, (a, b) in enumerate(E) if a in vs or b in vs)
            den = min(_prod(deg[y] + shift for y in vs if y != x0) for x0 in vs)
            groups = defaultdict(list)
            for s in family:
                groups[s & ~inc].append(s & c)
            for g in groups.values():
                pat = Counter(g)
                for i in range(len(E)):
                    if (c >> i) & 1:
                        hits = pat[c ^ (1 << i)]
                        conds += 1
                        ratio = Fraction(hits * den, len(g))
                        if ratio < 1:
                            print(f"VIOLATION {name} {label}: cycle {bin(c)} edge {i} ratio {ratio}", file=sys.stderr)
                            sys.exit(1)
                        margin = ratio if margin is None else min(margin, ratio)
        res[label] = (len(family), conds, margin)
    print(f"{name}: N={n} M={len(E)} trees={res['trees'][0]} forests={res['forests'][0]} short cycles={len(cycles)}; "
          f"tree bound: {res['trees'][1]} conditional inequalities, min ratio {res['trees'][2]}; "
          f"forest bound: {res['forests'][1]} conditional inequalities, min ratio {res['forests'][2]}")
    return deg, forests, cycles


def _prod(xs):
    out = 1
    for x in xs:
        out *= x
    return out


for name, (n, E) in GRAPHS.items():
    check(name, n, E)

# The referee's counterexample: uniform forests of the triangle with the unmodified tree bound.
forests, _ = acyclic_subsets(3, cyc(3))
full = (1 << 3) - 1
p_path = Fraction(sum(1 for f in forests if f == full ^ 1), len(forests))
tree_bound = Fraction(1, 2 * 2)
assert len(forests) == 7 and p_path == Fraction(1, 7) and p_path < tree_bound
print(f"triangle: {len(forests)} forests; Pr[F ∩ C = C - e] = {p_path} < {tree_bound} (tree bound), "
      f">= {Fraction(1, 9)} (forest bound (D+1)^-(l-1))")

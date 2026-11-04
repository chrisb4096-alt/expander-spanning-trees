# /// script
# dependencies = ["numpy", "scipy"]
# ///
"""Exact small-graph checks of the local lemma and Theorem "affinity bound for a split" (thm:cut) of the paper.

For a graph G, its uniform spanning tree P, and a balanced edge cut A|B, check the chain

  ||sqrt P||_op  <=  ||sqrt mu||_op  <=  prod_k gamma_k  <=  prod_k sqrt((1+sqrt(1-4 q_k))/2)  <=  prod_k sqrt(1 - q_k)

where mu is the law of the cycle statistics (1[A∩C_i ⊆ T], 1[B∩C_i ⊆ T])_i for a family of
vertex-disjoint mixed short cycles C_i, gamma_k is the largest operator norm of the conditional
2x2 matrices sqrt(kappa_k(.|past)), q_k = max(p_k^2, w_k), p_k = 1/min_{x0} prod_{y in V(C_k), y != x0} deg(y)
is the bound of the local lemma, and w_k = 1/prod_{y in V(C_k)} deg(y) is the star-framed weight.
For every past it also checks kappa(1,0) >= max(p_k, w_k^A) and kappa(0,1) >= max(p_k, w_k^B), where
w_k^A is the product of 1/deg(head(e)) over the A-edges of C_k in a fixed cyclic orientation (Lemma "framed cycles").

Also checks the local lower bound directly (min over outside configurations F).
Usage: uv run check_chain.py [--graph NAME] [--cuts K] [--search S] [--seed s]
"""
import argparse, itertools, os, subprocess, sys, tempfile
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import svds

HERE = os.path.dirname(os.path.abspath(__file__))


def prism(Hedges, nH):
    E = [(u, v) for u, v in Hedges] + [(u + nH, v + nH) for u, v in Hedges] + [(v, v + nH) for v in range(nH)]
    return 2 * nH, E


def cycle_edges(k):
    return [(i, (i + 1) % k) for i in range(k)]


def grid(h, w):
    idx = lambda i, j: i * w + j
    E = []
    for i in range(h):
        for j in range(w):
            if j + 1 < w: E.append((idx(i, j), idx(i, j + 1)))
            if i + 1 < h: E.append((idx(i, j), idx(i + 1, j)))
    return h * w, E


def square_graph(Hedges, nH):
    adj = [set() for _ in range(nH)]
    for u, v in Hedges: adj[u].add(v); adj[v].add(u)
    E = set()
    for v in range(nH):
        for u in adj[v]: E.add(tuple(sorted((u, v))))
        for u, w in itertools.combinations(sorted(adj[v]), 2): E.add((u, w))
    return nH, sorted(E)


K4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
K33 = [(i, 3 + j) for i in range(3) for j in range(3)]
GRAPHS = {
    "prism_C5": lambda: prism(cycle_edges(5), 5),
    "prism_C6": lambda: prism(cycle_edges(6), 6),
    "prism_K4": lambda: prism(K4, 4),
    "prism_K33": lambda: prism(K33, 6),
    "grid_4x4": lambda: grid(4, 4),
    "grid_3x5": lambda: grid(3, 5),
    "C10_square": lambda: square_graph(cycle_edges(10), 10),
    "C9_square": lambda: square_graph(cycle_edges(9), 9),
}


def enum_trees(n, E):
    exe = os.path.join(HERE, "enum_trees")
    if not os.path.exists(exe):
        subprocess.run(["g++", "-O2", "-o", exe, os.path.join(HERE, "enum_trees.cpp")], check=True)
    with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
        path = f.name
    inp = f"{n} {len(E)}\n" + "".join(f"{u} {v}\n" for u, v in E)
    subprocess.run([exe, path], input=inp.encode(), check=True, capture_output=True)
    T = np.fromfile(path, dtype=np.uint64)
    os.unlink(path)
    return T


def short_cycles(n, E, lmax=4):
    """All simple cycles of length 3..lmax, as frozensets of edge indices."""
    eid = {}
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(E):
        eid[(u, v)] = eid[(v, u)] = i
        adj[u].append(v); adj[v].append(u)
    cycles = set()
    def dfs(start, path):
        last = path[-1]
        for w in adj[last]:
            if w == start and len(path) >= 3:
                es = frozenset(eid[(path[i], path[(i + 1) % len(path)])] for i in range(len(path)))
                cycles.add(es)
            elif w > start and w not in path and len(path) < lmax:
                dfs(start, path + [w])
    for s in range(n):
        dfs(s, [s])
    return sorted(cycles, key=lambda c: (len(c), sorted(c)))


def verts_of(cyc, E):
    return sorted({x for i in cyc for x in E[i]})


def op_sqrtP(T, maskA, tau):
    a = T & np.uint64(maskA)
    b = T & ~np.uint64(maskA)
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    M = sp.csr_matrix((np.full(len(T), 1.0 / np.sqrt(tau)), (ia, ib)), shape=(len(ua), len(ub)))
    if min(M.shape) <= 2:
        return np.linalg.norm(M.toarray(), 2)
    return svds(M, k=1, return_singular_vectors=False)[0]


def stats_bits(T, cyc, maskA):
    cm = 0
    for i in cyc: cm |= 1 << i
    cA = np.uint64(cm & maskA); cB = np.uint64(cm & ~maskA & ((1 << 64) - 1))
    u = (T & cA) == cA
    v = (T & cB) == cB
    return u, v


def cycle_heads(cyc, E):
    """Cyclic orientation of a cycle given by edge indices: list of (edge index, head vertex)."""
    cyc = sorted(cyc)
    edges = {i: E[i] for i in cyc}
    first = cyc[0]
    a, b = edges[first]
    out = [(first, b)]
    cur, used = b, {first}
    while len(used) < len(cyc):
        nxt = next(i for i in cyc if i not in used and cur in edges[i])
        x, y = edges[nxt]
        cur = y if x == cur else x
        out.append((nxt, cur)); used.add(nxt)
    assert cur == a, "not a cycle"
    return out


def chain_check(n, E, T, maskA, cycles_sel, deg):
    tau = len(T)
    r = len(cycles_sel)
    U = np.zeros(len(T), dtype=np.int64); Vv = np.zeros(len(T), dtype=np.int64)
    gammas, frame_ok, qtheo, pq_emp = [], True, [], []
    for k, cyc in enumerate(cycles_sel):
        u, v = stats_bits(T, cyc, maskA)
        assert not np.any(u & v), "forbidden (1,1) pattern occurred"
        # conditional 2x2 laws given past (U,Vv)
        key = U * (1 << r) + Vv
        loc = u.astype(np.int64) * 2 + v.astype(np.int64)
        keys, inv = np.unique(key, return_inverse=True)
        C = np.zeros((len(keys), 4))
        np.add.at(C, (inv, loc), 1.0)
        C /= C.sum(axis=1, keepdims=True)
        vs = verts_of(cyc, E)
        ds = sorted(deg[x] for x in vs)
        pk = 1.0 / float(np.prod(ds[:-1]))
        wA = wB = 1.0
        for i, head in cycle_heads(cyc, E):
            if (maskA >> i) & 1:
                wA /= deg[head]
            else:
                wB /= deg[head]
        q = max(pk * pk, wA * wB)
        g = 0.0; pq = 1.0
        for row in C:
            Mx = np.sqrt(np.array([[row[0], row[1]], [row[2], row[3]]]))
            g = max(g, np.linalg.norm(Mx, 2))
            pq = min(pq, row[1] * row[2])
            # row[2] = Pr[(u,v)=(1,0) | past], row[1] = Pr[(0,1) | past]
            if row[2] < max(pk, wA) - 1e-12 or row[1] < max(pk, wB) - 1e-12:
                frame_ok = False
        gammas.append(g); qtheo.append(q); pq_emp.append(pq)
        U = U | (u.astype(np.int64) << k); Vv = Vv | (v.astype(np.int64) << k)
    # full statistic matrix
    Mu = np.zeros((1 << r, 1 << r))
    np.add.at(Mu, (U, Vv), 1.0)
    Mu /= tau
    op_mu = np.linalg.norm(np.sqrt(Mu), 2)
    return op_mu, gammas, qtheo, pq_emp, frame_ok


def local_lemma_check(n, E, T, cyc, deg):
    """min over F = T \\ I(C) and e in C of Pr[T ∩ C = C - e | F], vs proved bound."""
    vs = set(verts_of(cyc, E))
    Imask = 0
    for i, (x, y) in enumerate(E):
        if x in vs or y in vs: Imask |= 1 << i
    cm = 0
    for i in cyc: cm |= 1 << i
    F = T & ~np.uint64(Imask)
    TC = T & np.uint64(cm)
    keys, inv = np.unique(F, return_inverse=True)
    tot = np.bincount(inv)
    worst = 1.0
    for e in cyc:
        target = np.uint64(cm & ~(1 << e))
        hit = np.bincount(inv, weights=(TC == target).astype(float), minlength=len(keys))
        worst = min(worst, float(np.min(hit / tot)))
    ds = sorted(deg[x] for x in vs)
    return worst, 1.0 / float(np.prod(ds[:-1]))


def greedy_disjoint_mixed(cycles, E, maskA, rng=None):
    mixed = [c for c in cycles if 0 < sum((maskA >> i) & 1 for i in c) < len(c)]
    order = list(range(len(mixed)))
    best = []
    for trial in range(30 if rng is not None else 1):
        if rng is not None: rng.shuffle(order)
        used = set(); sel = []
        for j in order:
            vs = set(verts_of(mixed[j], E))
            if vs & used: continue
            sel.append(mixed[j]); used |= vs
        if len(sel) > len(best): best = sel
    return best


def random_balanced_mask(m, rng):
    k = rng.integers(m // 3 + 1, 2 * m // 3 + 1)
    idx = rng.choice(m, size=k, replace=False)
    mask = 0
    for i in idx: mask |= 1 << int(i)
    return mask


def local_search(T, m, tau, rng, steps):
    """Greedy single-edge flips (within balance) maximizing ||sqrt P||_op."""
    mask = random_balanced_mask(m, rng)
    cur = op_sqrtP(T, mask, tau)
    for _ in range(steps):
        i = int(rng.integers(m))
        new = mask ^ (1 << i)
        k = bin(new).count("1")
        if not (m / 3 < k <= 2 * m / 3): continue
        val = op_sqrtP(T, new, tau)
        if val >= cur: mask, cur = new, val
    return mask, cur


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default="prism_K4")
    ap.add_argument("--cuts", type=int, default=5)
    ap.add_argument("--search", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    n, E = GRAPHS[args.graph]()
    m = len(E)
    deg = [0] * n
    for u, v in E: deg[u] += 1; deg[v] += 1
    T = enum_trees(n, E)
    tau = len(T)
    cycles = short_cycles(n, E, 4)
    print(f"{args.graph}: n={n} m={m} tau={tau} maxdeg={max(deg)} short cycles={len(cycles)}")
    failures = 0
    # local lemma: every short cycle, every outside configuration
    worst_ratio = min(w / b for w, b in (local_lemma_check(n, E, T, c, deg) for c in cycles))
    print(f"local lemma: min_C,e,F Pr[T∩C=C-e | F] / proved bound = {worst_ratio:.3f} (must be >= 1)")
    if not worst_ratio >= 1 - 1e-12:
        print("VIOLATION: local-lemma ratio is below 1")
        failures += 1
    masks = [random_balanced_mask(m, rng) for _ in range(args.cuts)]
    tags = ["random"] * len(masks)
    for _ in range(args.search):
        mk, val = local_search(T, m, tau, rng, 300)
        masks.append(mk); tags.append("search")
    for mk, tag in zip(masks, tags):
        opP = op_sqrtP(T, mk, tau)
        sel = greedy_disjoint_mixed(cycles, E, mk, rng)
        if not sel:
            print(f"[{tag}] |A|={bin(mk).count('1')} opP={opP:.4f} no mixed short cycle"); continue
        op_mu, gam, qth, pq_emp, frame_ok = chain_check(n, E, T, mk, sel, deg)
        prodg = float(np.prod(gam))
        sharp = float(np.prod([np.sqrt((1 + np.sqrt(1 - 4 * q)) / 2) for q in qth]))
        theo = float(np.prod([np.sqrt(1 - q) for q in qth]))
        ok = (opP <= op_mu + 1e-9 and op_mu <= prodg + 1e-9 and prodg <= sharp + 1e-9 and sharp <= theo + 1e-12
              and frame_ok and min(np.array(pq_emp) - np.array(qth)) >= -1e-12)
        if not ok:
            failures += 1
        print(f"[{tag}] |A|={bin(mk).count('1')} r={len(sel)} opP={opP:.4f} <= opMu={op_mu:.4f} <= prod gamma={prodg:.4f} "
              f"<= theory={sharp:.6f} <= {theo:.6f}  min K(1,0)K(0,1)={min(pq_emp):.4f} (q_theory={min(qth):.5f})  "
              f"{'OK' if ok else 'VIOLATION'}")

    if failures:
        print(f"FAILED: {failures} check(s) violated the proved bounds.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

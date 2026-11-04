# /// script
# dependencies = ["numpy", "networkx"]
# ///
"""Checks of Section "Does subtraction help?" on small 2-connected graphs, by enumerating spanning trees.

For each graph and a few random balanced splits A|B (|A| = floor(M/2)), with a greedy family of vertex-disjoint
mixed cycles of length <= 6:
1. Theorem "exact rank": the explicit 2^m x 2^m submatrix of the proof (rows X_s, columns Y_t) contains exactly
   one spanning tree per row and per column, at t = 1 - s; and rank M_{A,B}(P) >= 2^m (numpy matrix_rank on the
   rows and columns that contain a tree; the other rows and columns are zero).
2. Lemma "blocks": with the cycles paired, every block has 2^g trees, its cell set R holds exactly one tree per
   row and column, and the sets R of different blocks are disjoint.
3. Theorem "approximate rank": for Q the best rank-r approximation (truncated SVD) of M_{A,B}(P), r = 1, 2, 4, 8,
   2 eps >= alpha_g (1 - r/2^g) for every g, where eps = 1/2 ||P - Q||_1.  For the normalized Q/sum(Q) (when
   sum(Q) > 0): Lemma lem:fact-n (overlap >= 1 - eps, ||Q||_1 <= 1 + 2 eps), the block overlap bound
   overlap <= delta(g) + sqrt(r ||Q||_1 / 2^g), and Theorem thm:signed-o2, [1 - eps - delta(g)]_+^2 <= r / 2^g.
4. Lemma "good pairs are frequent", exactly (Fractions): Pr[pair good | T cap E(C_l) for all cycles C_l of the
   other pairs] >= 2 p_i p_j with p_k = D^-(|C_k|-1), for every conditioning value of positive probability.
5. Theorem "normalized squares", the core inequality sum_x sqrt(P) ||d(x)|| <= beta sqrt(s) sqrt(Z) for random
   stacked amplitudes of rank <= s and for the truncated SVD of the square-root matrix, beta = ||sqrt M(P)||_op.
6. The counterexample after that theorem: det of the reduced Laplacian of x is 1 on spanning trees and 0 on
   random (N-1)-edge sets that are not trees.
Also Lemma "l1 distance to low rank" on random permutation matrices and rank-r matrices, with the equality case.
Exits with status 1 on any failure.
"""
import sys, random
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np
import networkx as nx
from common import named_graph, edge_list, enum_trees

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL:", msg)


def popcount(x):
    return bin(int(x)).count("1")


def graph(name):
    if name.startswith("prism"):
        G = nx.circular_ladder_graph(int(name[5:]))
        return nx.convert_node_labels_to_integers(G)
    if name.startswith("mobladder"):
        G = nx.circulant_graph(2 * int(name[9:]), [1, int(name[9:])])
        return nx.convert_node_labels_to_integers(G)
    return named_graph(name)


def is_forest_add(parent, u, v):
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    a, b = find(u), find(v)
    if a == b:
        return False
    parent[a] = b
    return True


def disjoint_mixed_cycles(G, E, idx, inA, rng, lmax=6):
    cyc = []
    for c in nx.simple_cycles(G, length_bound=lmax):
        es = [idx[tuple(sorted((c[i], c[(i + 1) % len(c)])))] for i in range(len(c))]
        if any(inA[e] for e in es) and any(not inA[e] for e in es):
            cyc.append((len(c), rng.random(), c, es))
    cyc.sort()
    used, out = set(), []
    for _, _, c, es in cyc:
        if used.isdisjoint(c):
            used.update(c)
            out.append(es)
    return out


def check_l1_lemma(rng):
    worst = np.inf
    for _ in range(3000):
        n = rng.randint(1, 8)
        r = rng.randint(0, n)
        Pi = np.eye(n)[np.random.permutation(n)]
        kind = rng.randint(0, 2)
        if kind == 0:
            Q = np.random.randn(n, r) @ np.random.randn(r, n) * rng.choice([0.01, 0.3, 1.0])
        elif kind == 1:
            U, S, Vt = np.linalg.svd(Pi + 0.2 * np.random.randn(n, n))
            Q = (U[:, :r] * S[:r]) @ Vt[:r]
        else:
            keep = np.zeros(n)
            keep[np.random.permutation(n)[:r]] = 1
            Q = Pi * keep[:, None] + (np.random.randn(n, r) @ np.random.randn(r, n)) * 1e-3 * (r > 0)
            if np.linalg.matrix_rank(Q) > r:
                Q = Pi * keep[:, None]
        val = np.abs(Pi - Q).sum() - (n - r)
        worst = min(worst, val)
        check(val >= -1e-9, f"l1 lemma n={n} r={r} kind={kind} slack={val}")
    # equality case: keep r of the ones
    Pi = np.eye(6)
    Q = np.diag([1, 1, 0, 0, 0, 0.0])
    check(abs(np.abs(Pi - Q).sum() - 4) < 1e-12, "l1 lemma equality case")
    print(f"l1-to-low-rank lemma: 3000 random cases, smallest slack {worst:.3g} (equality case attains 0)")


def run_graph(name, nsplits, seed):
    G = graph(name)
    n, E = G.number_of_nodes(), edge_list(G)
    M = len(E)
    idx = {e: i for i, e in enumerate(E)}
    D = max(d for _, d in G.degree())
    T = enum_trees(n, E)
    tau = len(T)
    Tset = set(int(t) for t in T)
    rng = random.Random(seed)
    np.random.seed(seed)
    # counterexample after the squares theorem
    Lfull = None
    def red_det(mask):
        L = np.zeros((n, n))
        for e, (u, v) in enumerate(E):
            if mask >> e & 1:
                L[u, u] += 1; L[v, v] += 1; L[u, v] -= 1; L[v, u] -= 1
        return np.linalg.det(L[1:, 1:])
    for t in rng.sample(sorted(Tset), min(200, tau)):
        check(abs(red_det(t) - 1) < 1e-8, f"{name}: det on a tree")
    bad = 0
    for _ in range(300):
        s = rng.sample(range(M), n - 1)
        mask = sum(1 << e for e in s)
        if mask in Tset:
            continue
        bad += 1
        check(abs(red_det(mask)) < 1e-8, f"{name}: det on a non-tree")
    stats = []
    for sp in range(nsplits):
        A = set(rng.sample(range(M), M // 2))
        inA = [e in A for e in range(M)]
        mA = sum(1 << e for e in A)
        mB = ((1 << M) - 1) ^ mA
        cyc = disjoint_mixed_cycles(G, E, idx, inA, rng)
        m = len(cyc)
        # ---- 1. explicit permutation submatrix
        ef = []
        for es in cyc:
            a = [e for e in es if inA[e]]; b = [e for e in es if not inA[e]]
            ef.append((rng.choice(a), rng.choice(b)))
        parent = list(range(n))
        T0 = 0
        for es, (e, f) in zip(cyc, ef):
            for x in es:
                if x != e:
                    ok = is_forest_add(parent, *E[x]); check(ok, f"{name}: cycle paths not a forest"); T0 |= 1 << x
        order = list(range(M)); rng.shuffle(order)
        for x in order:
            if is_forest_add(parent, *E[x]):
                T0 |= 1 << x
        check(T0 in Tset, f"{name}: T0 not a tree")
        for e, f in ef:
            check(not T0 >> e & 1 and T0 >> f & 1, f"{name}: e_i in T0 or f_i not in T0")
        Fm = T0
        for e, f in ef:
            Fm &= ~(1 << f)
        for s in range(2 ** m):
            X = (Fm & mA) | sum(1 << ef[i][0] for i in range(m) if s >> i & 1)
            for t in range(2 ** m):
                Y = (Fm & mB) | sum(1 << ef[i][1] for i in range(m) if t >> i & 1)
                check(((X | Y) in Tset) == (t == (2 ** m - 1) ^ s), f"{name}: permutation submatrix cell")
        # ---- rank of M_{A,B}(P)
        ra = {}; cb = {}
        for t in Tset:
            ra.setdefault(t & mA, len(ra)); cb.setdefault(t & mB, len(cb))
        Mat = np.zeros((len(ra), len(cb)))
        for t in Tset:
            Mat[ra[t & mA], cb[t & mB]] = 1.0 / tau
        sv = np.linalg.svd(Mat, compute_uv=False)
        rank = int((sv > sv[0] * 1e-10).sum())
        check(rank >= 2 ** m, f"{name}: rank {rank} < 2^{m}")
        # ---- 2-4. blocks, approximation, good pairs
        npairs = m // 2
        pairs = [(2 * k, 2 * k + 1) for k in range(npairs)]
        nblocks = 0; worst_good = None; approx_min = None; norm_min = None
        if npairs >= 1:
            cmask = [sum(1 << x for x in es) for es in cyc]
            def pattern(t, k):
                i, j = pairs[k]
                ci, cj = cmask[i], cmask[j]
                arcs = (ci & ~(1 << ef[i][0]) & ~(1 << ef[i][1])) | (cj & ~(1 << ef[j][0]) & ~(1 << ef[j][1]))
                if t & arcs != arcs:
                    return 0
                if t & ci == ci & ~(1 << ef[i][1]) and t & cj == cj & ~(1 << ef[j][0]):
                    return 1
                if t & ci == ci & ~(1 << ef[i][0]) and t & cj == cj & ~(1 << ef[j][1]):
                    return 2
                return 0
            blocks = defaultdict(list)
            for t in Tset:
                good = tuple(k for k in range(npairs) if pattern(t, k))
                Ft = t
                for k in good:
                    i, j = pairs[k]
                    for x in (ef[i][0], ef[i][1], ef[j][0], ef[j][1]):
                        Ft &= ~(1 << x)
                blocks[(good, Ft)].append(t)
            allcells = set(); ncell = 0
            gcount = defaultdict(int)
            for (good, Ft), ts in blocks.items():
                g = len(good)
                gcount[g] += len(ts)
                check(len(ts) == 2 ** g, f"{name}: block size {len(ts)} != 2^{g}")
                rows = []; cols = []
                for s in range(2 ** g):
                    X = Ft & mA; Y = Ft & mB
                    for b, k in enumerate(good):
                        i, j = pairs[k]
                        X |= 1 << (ef[i][0] if s >> b & 1 else ef[j][0])
                        Y |= 1 << (ef[i][1] if s >> b & 1 else ef[j][1])
                    rows.append(X); cols.append(Y)
                for X in rows:
                    check(sum((X | Y) in Tset for Y in cols) == 1, f"{name}: block row not permutation")
                for Y in cols:
                    check(sum((X | Y) in Tset for X in rows) == 1, f"{name}: block col not permutation")
                for X in rows:
                    for Y in cols:
                        allcells.add((X, Y)); ncell += 1
            nblocks = len(blocks)
            check(len(allcells) == ncell, f"{name}: block cell sets overlap")
            alpha = [sum(c for gg, c in gcount.items() if gg >= g) / tau for g in range(npairs + 1)]
            U, S, Vt = np.linalg.svd(Mat, full_matrices=False)
            for r in (1, 2, 4, 8):
                Q = (U[:, :r] * S[:r]) @ Vt[:r]
                twoeps = np.abs(Mat - Q).sum()
                for g in range(npairs + 1):
                    slack = twoeps - alpha[g] * (1 - r / 2 ** g)
                    check(slack >= -1e-12, f"{name}: approx thm r={r} g={g} slack={slack}")
                    approx_min = slack if approx_min is None else min(approx_min, slack)
                # normalized output: Fact N, block overlap lemma and the normalized bound
                lam = Q.sum()
                if lam > 0:
                    Qn = Q / lam
                    eps_n = 0.5 * np.abs(Mat - Qn).sum()
                    L1 = np.abs(Qn).sum()
                    om = np.minimum(Mat, np.maximum(Qn, 0)).sum()
                    check(om >= 1 - eps_n - 1e-12 and L1 <= 1 + 2 * eps_n + 1e-12, f"{name}: Fact N r={r}")
                    for g in range(npairs + 1):
                        delta = 1 - alpha[g]
                        check(om <= delta + np.sqrt(r * L1 / 2 ** g) + 1e-12, f"{name}: block overlap r={r} g={g}")
                        nslack = r / 2 ** g - max(0.0, 1 - eps_n - delta) ** 2
                        check(nslack >= -1e-12, f"{name}: normalized approx thm r={r} g={g}")
                        norm_min = nslack if norm_min is None else min(norm_min, nslack)
            # good pairs, conditioned on the intersections with all cycles of the other pairs
            for k in range(npairs):
                i, j = pairs[k]
                p2 = 2 * Fr(1, D ** (len(cyc[i]) - 1)) * Fr(1, D ** (len(cyc[j]) - 1))
                others = [cmask[l] for kk in range(npairs) if kk != k for l in pairs[kk]]
                cnt = defaultdict(lambda: [0, 0])
                for t in Tset:
                    key = tuple(t & c for c in others)
                    cnt[key][0] += 1
                    cnt[key][1] += pattern(t, k) > 0
                for key, (tot, gd) in cnt.items():
                    ratio = Fr(gd, tot) / p2
                    check(ratio >= 1, f"{name}: good-pair bound fails ({gd}/{tot} < {p2})")
                    worst_good = ratio if worst_good is None else min(worst_good, ratio)
        # ---- 5. squares
        W = np.sqrt(Mat)
        beta = np.linalg.svd(W, compute_uv=False)[0]
        sq_max = 0.0
        for s_ in (1, 2, 4, 8):
            for K in (1, 3):
                U_ = np.random.randn(W.shape[0], s_)
                Ds = [U_ @ np.random.randn(s_, W.shape[1]) for _ in range(K)]
                nd = np.sqrt(sum(Dm ** 2 for Dm in Ds))
                lhs = (W * nd).sum(); rhs = beta * np.sqrt(s_) * np.sqrt((nd ** 2).sum())
                check(lhs <= rhs * (1 + 1e-12), f"{name}: squares inequality s={s_} K={K}")
                sq_max = max(sq_max, lhs / rhs)
            Uw, Sw, Vw = np.linalg.svd(W, full_matrices=False)
            Dt = (Uw[:, :s_] * Sw[:s_]) @ Vw[:s_]
            lhs = (W * np.abs(Dt)).sum(); rhs = beta * np.sqrt(s_) * np.linalg.norm(Dt)
            check(lhs <= rhs * (1 + 1e-12), f"{name}: squares inequality (SVD) s={s_}")
            sq_max = max(sq_max, lhs / rhs)
        stats.append((m, rank, nblocks, approx_min, norm_min, worst_good, beta, sq_max))
    print(f"{name}: N={n} M={M} D={D} trees={tau}; non-tree (N-1)-sets tested {bad}")
    for m, rank, nb, am, nm, wg, beta, sq in stats:
        wgs = f"{float(wg):.3g}" if wg is not None else "-"
        ams = f"{am:.3g}" if am is not None else "-"
        nms = f"{nm:.3g}" if nm is not None else "-"
        print(f"   m={m} rank={rank} (>= 2^m={2**m}) blocks={nb} approx-slack-min={ams} normalized-slack-min={nms} "
              f"good-pair ratio min={wgs} beta={beta:.3f} squares lhs/rhs max={sq:.3f}")


if __name__ == "__main__":
    rng = random.Random(0); np.random.seed(0)
    check_l1_lemma(rng)
    for i, name in enumerate(["K33", "cube", "prism5", "petersen", "prism6", "mobladder6"]):
        run_graph(name, 3, 10 + i)
    for i, name in enumerate(["prism8"]):
        run_graph(name, 1, 20 + i)
    if FAIL:
        print(f"{len(FAIL)} FAILURES")
        sys.exit(1)
    print("ALL PASS")

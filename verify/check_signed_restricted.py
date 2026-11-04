# /// script
# dependencies = ["numpy", "networkx"]
# ///
"""Checks for Section "Restricted unstructured classes" and the sign split (sec-signed.tex).

  S1  Lemma lem:digon on K_{t,t} with the placements of Theorems thm:signed-depth and thm:signed-formula:
      the constants of Raz and Yehudayoff (pairs in 2x2 blocks, the other block cells and a matching set to 1)
      give the zero function, while the constants of Lemma lem:digon give prod_i [y_i xor z_i], rank 2^m;
      for Raz-type pairs in a common row or column, completed by spare rows and columns as in the proof of
      Theorem thm:signed-formula, likewise
  S2  Lemma lem:digon on random cubic graphs with vertex-disjoint cycles; signed support-exact values keep rank
      2^m; for ST_G the restriction is prod_i (x_e + x_f + (|C_i| - 2) x_e x_f)
  S3  Lemma lem:minor-restrict: Petersen -> K_{3,3}, and random minors of random cubic graphs
  S4  Lemma lem:frontier on random aligned circuits (rank <= |I| + 1 <= number of nodes, for every vtree node)
  S5  Lemma lem:vtree-refine and the packing in Theorem thm:signed-sums: |P| = |Q| >= t^2 3^-r, P and Q separated
      by a node of every vtree, m >= t 3^-r / 4 - 1 four-cycles, and full rank end to end for small m
  S6  The sign split after Theorem thm:signed-cancel on random signed smooth decomposable circuits with real
      leaf functions: f = g - h, g + h = f_{C^+}, g, h >= 0, scopes unchanged, size <= 6|C|; the overlap step
      BC(P, g / sum g) >= (1 - eps) sqrt(lambda / sum g) for g >= f^+; and ST(triangle)(1,1,1) = 3

Exits 1 on any failure.
"""
import itertools
import random
import sys

import networkx as nx
import numpy as np

FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL:", msg, flush=True)


class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return False
        self.p[a] = b
        return True


def is_spanning_tree(n, E, chosen):
    chosen = list(chosen)
    if len(chosen) != n - 1:
        return False
    d = DSU(n)
    for i in chosen:
        u, v = E[i]
        if u == v or not d.union(u, v):
            return False
    return True


def count_trees(n, E, chosen):
    L = np.zeros((n, n))
    for i in chosen:
        u, v = E[i]
        L[u, u] += 1
        L[v, v] += 1
        L[u, v] -= 1
        L[v, u] -= 1
    return int(round(np.linalg.det(L[1:, 1:])))


def restricted_table(n, E, ones, free, fn=is_spanning_tree):
    out = {}
    for bits in itertools.product((0, 1), repeat=len(free)):
        out[bits] = int(fn(n, E, list(ones) + [free[i] for i, b in enumerate(bits) if b]))
    return out


def digon_restriction(n, E, cycles):
    """Proof of Lemma lem:digon: T contains every C_i - f_i; K = T minus the e_i. Returns (K, pairs)."""
    d = DSU(n)
    T = []
    fs = {f for _, _, f in cycles}
    for cyc, e, f in cycles:
        for i in cyc:
            if i != f:
                assert d.union(*E[i]), "C_i - f_i is not a path"
                T.append(i)
    incyc = {i for cyc, _, _ in cycles for i in cyc}
    for i, (u, v) in enumerate(E):
        if i not in incyc and i not in fs and d.union(u, v):
            T.append(i)
    assert len(T) == n - 1
    es = {e for _, e, _ in cycles}
    return [i for i in T if i not in es], [(e, f) for _, e, f in cycles]


def value_matrix(table, m):
    M = np.zeros((2 ** m, 2 ** m))
    for bits, v in table.items():
        M[int("".join(map(str, bits[:m])), 2), int("".join(map(str, bits[m:])), 2)] = v
    return M


def xor_ok(table, m):
    return all(v == int(all(b[i] != b[m + i] for i in range(m))) for b, v in table.items())


def digon_table(n, E, cycles, fn=is_spanning_tree):
    K, pairs = digon_restriction(n, E, cycles)
    free = [e for e, _ in pairs] + [f for _, f in pairs]
    return restricted_table(n, E, K, free, fn), len(pairs), K


def kbb(t):
    return [(i, t + j) for i in range(t) for j in range(t)]


# ------------------------------------------------------------------ S1
def s1_placements(rng):
    for t, m in ((5, 2), (7, 3), (9, 4), (12, 5)):
        E = kbb(t)
        ix = lambda rc: rc[0] * t + rc[1]  # noqa: E731
        for _ in range(4):
            rows, cols = rng.sample(range(t), 2 * m), rng.sample(range(t), 2 * m)
            blocks = []
            for i in range(m):
                r1, r2, c1, c2 = rows[2 * i], rows[2 * i + 1], cols[2 * i], cols[2 * i + 1]
                if rng.random() < 0.5:
                    y, z, ones = (r1, c1), (r1, c2), [(r2, c1), (r2, c2)]
                else:
                    y, z, ones = (r1, c1), (r2, c1), [(r1, c2), (r2, c2)]
                blocks.append((y, z, ones, (r1, r2, c1, c2)))
            ident = list(zip([r for r in range(t) if r not in rows], [c for c in range(t) if c not in cols]))
            free = [ix(b[0]) for b in blocks] + [ix(b[1]) for b in blocks]
            onesRY = [ix(o) for b in blocks for o in b[2]] + [ix(p) for p in ident]
            tab = restricted_table(2 * t, E, onesRY, free)
            check(all(v == 0 for v in tab.values()), f"S1 t={t}: Raz-Yehudayoff constants should kill the indicator")
            cycles = []
            for y, z, ones, (r1, r2, c1, c2) in blocks:
                cycles.append(([ix((r1, c1)), ix((r1, c2)), ix((r2, c1)), ix((r2, c2))], ix(y), ix(z)))
            tab, mm, _ = digon_table(2 * t, E, cycles)
            rk = np.linalg.matrix_rank(value_matrix(tab, mm))
            check(xor_ok(tab, mm) and rk == 2 ** mm, f"S1 t={t}: digon constants")
        print(f"S1 RY placement t={t} m={m}: their constants give 0; Lemma lem:digon gives prod XOR, rank 2^{m}")
    for t in (9, 12, 15):
        m = -(-t // 4)
        while 3 * m > t:
            m -= 1
        E = kbb(t)
        ix = lambda rc: rc[0] * t + rc[1]  # noqa: E731
        for _ in range(4):
            # m pairs, each in a common row or a common column, using disjoint rows and columns
            free_r, free_c = list(range(t)), list(range(t))
            rng.shuffle(free_r)
            rng.shuffle(free_c)
            pairs = []
            for _ in range(m):
                if rng.random() < 0.5:
                    r, c, c2 = free_r.pop(), free_c.pop(), free_c.pop()
                    pairs.append(("row", (r, c), (r, c2)))
                else:
                    c, r, r2 = free_c.pop(), free_r.pop(), free_r.pop()
                    pairs.append(("col", (r, c), (r2, c)))
            cycles = []
            for kind, (r1, c1), (r2, c2) in pairs:
                if kind == "row":
                    rr = free_r.pop()
                    cyc = [ix((r1, c1)), ix((r1, c2)), ix((rr, c1)), ix((rr, c2))]
                else:
                    cc = free_c.pop()
                    cyc = [ix((r1, c1)), ix((r2, c1)), ix((r1, cc)), ix((r2, cc))]
                cycles.append((cyc, ix((r1, c1)), ix((r2, c2))))
            tab, mm, _ = digon_table(2 * t, E, cycles)
            rk = np.linalg.matrix_rank(value_matrix(tab, mm))
            check(xor_ok(tab, mm) and rk == 2 ** mm, f"S1 Raz t={t}")
        print(f"S1 Raz placement t={t} m={m} (2m <= t - m): spare rows and columns suffice, prod XOR, rank 2^{m}")


# ------------------------------------------------------------------ S2
def disjoint_cycles(G, E, rng, mmax):
    eid = {frozenset(e): i for i, e in enumerate(E)}
    used, out = set(), []
    order = list(range(len(E)))
    rng.shuffle(order)
    for e in order:
        a, b = E[e]
        if a in used or b in used:
            continue
        H = G.copy()
        H.remove_edge(a, b)
        H.remove_nodes_from(used)
        try:
            path = nx.shortest_path(H, a, b)
        except nx.NetworkXNoPath:
            continue
        cyc = [eid[frozenset(p)] for p in zip(path, path[1:])] + [e]
        out.append((cyc, e, rng.choice(cyc[:-1])))
        used |= set(path)
        if len(out) == mmax:
            break
    return out


def s2_digon(rng):
    for N, seed in ((20, 1), (30, 2), (40, 3), (60, 4)):
        G = nx.random_regular_graph(3, N, seed=seed)
        while not nx.is_connected(G):
            seed += 100
            G = nx.random_regular_graph(3, N, seed=seed)
        E = [tuple(sorted(e)) for e in G.edges()]
        cycles = disjoint_cycles(G, E, rng, 5)
        tab, m, K = digon_table(N, E, cycles)
        rk = np.linalg.matrix_rank(value_matrix(tab, m))
        w = {b: (rng.choice([-1, 1]) * rng.uniform(0.1, 3) if v else 0.0) for b, v in tab.items()}
        rk2 = np.linalg.matrix_rank(value_matrix(w, m))
        # ST_G restricted: number of spanning trees inside K + chosen edges
        st = restricted_table(N, E, K, [e for _, e, _ in cycles] + [f for _, _, f in cycles], count_trees)
        lens = [len(c) for c, _, _ in cycles]
        ok_st = all(v == int(np.prod([(b[i] + b[m + i] + (lens[i] - 2) * b[i] * b[m + i]) for i in range(m)]))
                    for b, v in st.items())
        check(xor_ok(tab, m) and rk == 2 ** m and rk2 == 2 ** m and ok_st, f"S2 N={N}")
        print(f"S2 random cubic N={N}: m={m} disjoint cycles of lengths {lens}: indicator = prod XOR (rank {rk}), "
              f"signed support-exact values rank {rk2}, ST_G = prod (x_e + x_f + (|C|-2) x_e x_f): {ok_st}")


# ------------------------------------------------------------------ S3
def quotient(n, E, forest):
    d = DSU(n)
    for i in forest:
        assert d.union(*E[i])
    roots = sorted({d.find(v) for v in range(n)})
    lab = {r: k for k, r in enumerate(roots)}
    return [lab[d.find(v)] for v in range(n)], len(roots)


def minor_ok(n, E, forest, keep):
    comp, k = quotient(n, E, forest)
    EH = [(comp[E[i][0]], comp[E[i][1]]) for i in keep]
    tab = restricted_table(n, E, forest, keep)
    return all(v == int(is_spanning_tree(k, EH, [j for j, b in enumerate(bits) if b])) for bits, v in tab.items())


def s3_minors(rng):
    G = nx.petersen_graph()
    E = [tuple(sorted(e)) for e in G.edges()]
    found = False
    for _ in range(4000):
        order = list(range(len(E)))
        rng.shuffle(order)
        d = DSU(10)
        forest = []
        for i in order:
            if len(forest) == 4:
                break
            if d.union(*E[i]):
                forest.append(i)
        comp, k = quotient(10, E, forest)
        cand = {}
        for i, (u, v) in enumerate(E):
            if comp[u] != comp[v] and i not in forest:
                cand.setdefault(frozenset((comp[u], comp[v])), i)
        for A in itertools.combinations(range(6), 3):
            B = [x for x in range(6) if x not in A]
            if all(frozenset((a, b)) in cand for a in A for b in B):
                check(minor_ok(10, E, forest, [cand[frozenset((a, b))] for a in A for b in B]), "S3 Petersen")
                found = True
                break
        if found:
            break
    check(found, "S3 no K33 minor model found in the Petersen graph")
    cnt = 0
    for seed in range(8):
        G = nx.random_regular_graph(3, 14, seed=seed)
        if not nx.is_connected(G):
            continue
        E = [tuple(sorted(e)) for e in G.edges()]
        order = list(range(len(E)))
        rng.shuffle(order)
        d = DSU(14)
        forest = []
        target = rng.randint(4, 7)
        for i in order:
            if len(forest) == 14 - target:
                break
            if d.union(*E[i]):
                forest.append(i)
        comp, k = quotient(14, E, forest)
        cand = {}
        for i, (u, v) in enumerate(E):
            if comp[u] != comp[v] and i not in forest:
                cand.setdefault(frozenset((comp[u], comp[v])), i)
        keep = list(cand.values())[:14]
        check(minor_ok(14, E, forest, keep), f"S3 random minor seed={seed}")
        cnt += 1
    print(f"S3 Lemma lem:minor-restrict: Petersen -> K_3,3 ({'ok' if found else 'NOT FOUND'}) and {cnt} random minors")


# ------------------------------------------------------------------ S4
def random_vtree(leaves, rng):
    leaves = list(leaves)
    rng.shuffle(leaves)

    def build(lo, hi):
        if hi - lo == 1:
            return leaves[lo]
        mid = rng.randint(lo + 1, hi - 1)
        return (build(lo, mid), build(mid, hi))
    return build(0, len(leaves))


def vtree_nodes(T):
    out = []

    def rec(x):
        s = rec(x[0]) | rec(x[1]) if isinstance(x, tuple) else frozenset([x])
        out.append((x, s))
        return s
    rec(T)
    return out


def random_aligned_circuit(nv, rng):
    V = random_vtree(range(nv), rng)
    gates = []      # ("leaf", var, (a0, a1)) | ("const", c) | ("sum", [(w, child)]) | ("prod", a, b)
    scope = []

    def add(g, s):
        gates.append(g)
        scope.append(s)
        return len(gates) - 1

    def visit(x):
        if not isinstance(x, tuple):
            for _ in range(rng.randint(1, 3)):
                add(("leaf", x, (rng.uniform(-2, 2), rng.uniform(-2, 2))), frozenset([x]))
            return frozenset([x])
        s1, s2 = visit(x[0]), visit(x[1])
        s = s1 | s2
        for _ in range(rng.randint(1, 3)):
            c1 = [i for i in range(len(gates)) if scope[i] <= s1]
            c2 = [i for i in range(len(gates)) if scope[i] <= s2]
            add(("prod", rng.choice(c1), rng.choice(c2)), None)
            scope[-1] = scope[gates[-1][1]] | scope[gates[-1][2]]
        for _ in range(rng.randint(0, 2)):
            pool = [i for i in range(len(gates)) if scope[i] <= s]
            kids = rng.sample(pool, min(len(pool), rng.randint(1, 4)))
            add(("sum", [(rng.uniform(-2, 2), c) for c in kids]), frozenset().union(*[scope[c] for c in kids]))
        if rng.random() < 0.3:
            add(("const", rng.uniform(-1, 1)), frozenset())
        return s
    visit(V)
    pool = list(range(len(gates)))
    kids = rng.sample(pool, min(len(pool), 4))
    add(("sum", [(rng.uniform(-2, 2), c) for c in kids]), frozenset().union(*[scope[c] for c in kids]))
    return V, gates, scope


def eval_gates(gates, x):
    val = []
    for g in gates:
        if g[0] == "leaf":
            val.append(g[2][x[g[1]]])
        elif g[0] == "const":
            val.append(g[1])
        elif g[0] == "prod":
            val.append(val[g[1]] * val[g[2]])
        else:
            val.append(sum(w * val[c] for w, c in g[1]))
    return val


def s4_frontier(rng):
    nv = 8
    worst, cnt = 0, 0
    for _ in range(60):
        V, gates, scope = random_aligned_circuit(nv, rng)
        root = len(gates) - 1
        xs = list(itertools.product((0, 1), repeat=nv))
        f = {x: eval_gates(gates, x)[root] for x in xs}
        parents = [[] for _ in gates]
        for i, g in enumerate(gates):
            kids = [c for _, c in g[1]] if g[0] == "sum" else ([g[1], g[2]] if g[0] == "prod" else [])
            for c in kids:
                parents[c].append(i)
        for _, S in vtree_nodes(V):
            Sl = sorted(S)
            Tl = [v for v in range(nv) if v not in S]
            if not Tl:
                continue
            M = np.zeros((2 ** len(Sl), 2 ** len(Tl)))
            for x in xs:
                M[int("".join(str(x[v]) for v in Sl), 2), int("".join(str(x[v]) for v in Tl), 2)] = f[x]
            rk = np.linalg.matrix_rank(M, tol=1e-8 * max(1.0, np.abs(M).max()))
            strad = [bool(scope[i] & S) and not scope[i] <= S for i in range(len(gates))]
            inside = [scope[i] <= S and bool(scope[i]) for i in range(len(gates))]
            I = [i for i in range(len(gates)) if inside[i] and any(strad[p] for p in parents[i])]
            ok = rk <= len(I) + 1 <= len(gates)
            check(ok, "S4 frontier rank")
            worst = max(worst, rk)
            cnt += 1
    print(f"S4 Lemma lem:frontier: {cnt} (circuit, vtree node) pairs, rank <= |I| + 1 <= number of nodes; "
          f"largest rank {worst}")


# ------------------------------------------------------------------ S5
def balanced_node(scopes, U):
    U = set(U)
    for S in scopes:
        if len(U) / 3 < len(S & U) <= 2 * len(U) / 3:
            return S
    raise AssertionError("no balanced node")


def refine(vts, E0):
    S = balanced_node(vts[0], E0)
    P, Q = sorted(set(E0) & S), sorted(set(E0) - S)
    n = min(len(P), len(Q))
    P, Q = P[:n], Q[:n]
    for scopes in vts[1:]:
        S = balanced_node(scopes, set(P) | set(Q))
        a, b = [x for x in P if x in S], [x for x in P if x not in S]
        c, d = [x for x in Q if x in S], [x for x in Q if x not in S]
        P, Q = (a, d) if min(len(a), len(d)) >= min(len(b), len(c)) else (b, c)
        n = min(len(P), len(Q))
        P, Q = P[:n], Q[:n]
    return P, Q


def separated(scopes, P, Q):
    P, Q = set(P), set(Q)
    return any((P <= S and not Q & S) or (Q <= S and not P & S) for S in scopes)


def pack(t, P, Q):
    E = kbb(t)
    used, out = set(), []
    for e in P:
        r1, c1 = E[e]
        if r1 in used or c1 in used:
            continue
        for f in Q:
            r2, c2 = E[f]
            if {r1, c1, r2, c2} & used or r1 == r2 or c1 == c2:
                continue
            out.append(([e, r2 * t + (c1 - t), f, r1 * t + (c2 - t)], e, f))
            used |= {r1, c1, r2, c2}
            break
    return out


def s5_refine(rng):
    for t in (10, 16, 24):
        E0 = list(range(t * t))
        for r in (1, 2, 3, 4):
            wn, wm = None, None
            for _ in range(5):
                vts = [[s for _, s in vtree_nodes(random_vtree(E0, rng))] for _ in range(r)]
                P, Q = refine(vts, E0)
                check(len(P) == len(Q) >= t * t / 3 ** r and not set(P) & set(Q), f"S5 sizes t={t} r={r}")
                check(all(separated(sc, P, Q) for sc in vts), f"S5 separation t={t} r={r}")
                gad = pack(t, P, Q)
                check(len(gad) >= t / (4 * 3 ** r) - 1, f"S5 packing t={t} r={r}")
                wn = len(P) if wn is None else min(wn, len(P))
                wm = len(gad) if wm is None else min(wm, len(gad))
                if 1 <= len(gad) <= 6:
                    tab, mm, _ = digon_table(2 * t, kbb(t), gad)
                    check(xor_ok(tab, mm) and np.linalg.matrix_rank(value_matrix(tab, mm)) == 2 ** mm,
                          f"S5 end to end t={t} r={r}")
            print(f"S5 t={t} r={r}: min |P| = {wn} (bound {t * t / 3 ** r:.1f}), min m = {wm} "
                  f"(bound {t / (4 * 3 ** r) - 1:.2f})")


# ------------------------------------------------------------------ S6
def random_smooth_dc(nv, rng):
    gates = []

    def build(vs):
        if len(vs) == 1:
            k = []
            for _ in range(rng.randint(1, 2)):
                gates.append(("leaf", vs[0], (rng.uniform(-2, 2), rng.uniform(-2, 2))))
                k.append(len(gates) - 1)
            gates.append(("sum", [(rng.uniform(-2, 2), c) for c in k]))
            return len(gates) - 1
        kids = []
        for _ in range(rng.randint(1, 3)):
            vv = list(vs)
            rng.shuffle(vv)
            cut = rng.randint(1, len(vv) - 1)
            a, b = build(vv[:cut]), build(vv[cut:])
            gates.append(("prod", a, b))
            kids.append(len(gates) - 1)
        gates.append(("sum", [(rng.uniform(-2, 2), c) for c in kids]))
        return len(gates) - 1
    return gates, build(list(range(nv)))


def scopes_of(gates):
    sc = []
    for g in gates:
        if g[0] == "leaf":
            sc.append(frozenset([g[1]]))
        elif g[0] == "prod":
            check(not sc[g[1]] & sc[g[2]], "S6 decomposability")
            sc.append(sc[g[1]] | sc[g[2]])
        else:
            ks = [sc[c] for _, c in g[1]]
            check(len(ks) > 0 and all(k == ks[0] for k in ks), "S6 smoothness")
            sc.append(ks[0])
    return sc


def sign_split(gates):
    new, idx = [], {}
    for i, g in enumerate(gates):
        if g[0] == "leaf":
            a0, a1 = g[2]
            new.append(("leaf", g[1], (max(a0, 0), max(a1, 0))))
            new.append(("leaf", g[1], (max(-a0, 0), max(-a1, 0))))
            idx[i] = (len(new) - 2, len(new) - 1)
        elif g[0] == "prod":
            (ap, an), (bp, bn) = idx[g[1]], idx[g[2]]
            q = []
            for a, b in ((ap, bp), (an, bn), (ap, bn), (an, bp)):
                new.append(("prod", a, b))
                q.append(len(new) - 1)
            new.append(("sum", [(1.0, q[0]), (1.0, q[1])]))
            new.append(("sum", [(1.0, q[2]), (1.0, q[3])]))
            idx[i] = (len(new) - 2, len(new) - 1)
        else:
            pos, neg = [], []
            for w, c in g[1]:
                cp, cn = idx[c]
                pos.append((abs(w), cp if w >= 0 else cn))
                neg.append((abs(w), cn if w >= 0 else cp))
            new.append(("sum", pos))
            new.append(("sum", neg))
            idx[i] = (len(new) - 2, len(new) - 1)
    return new, idx


def abs_circuit(gates):
    out = []
    for g in gates:
        if g[0] == "leaf":
            out.append(("leaf", g[1], (abs(g[2][0]), abs(g[2][1]))))
        elif g[0] == "sum":
            out.append(("sum", [(abs(w), c) for w, c in g[1]]))
        else:
            out.append(g)
    return out


def s6_sign_split(rng):
    for _ in range(40):
        nv = rng.randint(2, 7)
        gates, root = random_smooth_dc(nv, rng)
        sc = scopes_of(gates)
        new, idx = sign_split(gates)
        sc2 = scopes_of(new)
        check(all(sc2[p] == sc2[q] == sc[i] for i, (p, q) in idx.items()), "S6 scopes changed")
        check(len(new) <= 6 * len(gates), "S6 size")
        ab = abs_circuit(gates)
        p, q = idx[root]
        for x in itertools.product((0, 1), repeat=nv):
            f = eval_gates(gates, x)[root]
            fa = eval_gates(ab, x)[root]
            v = eval_gates(new, x)
            tol = 1e-9 * (1 + abs(fa))
            check(abs(v[p] - v[q] - f) < tol and abs(v[p] + v[q] - fa) < tol and min(v[p], v[q]) >= -1e-12,
                  "S6 sign split values")
    print("S6 40 random signed smooth decomposable circuits: f = g - h, g + h = f_{C^+}, g, h >= 0, "
          "scopes unchanged, size <= 6|C|")
    worst = np.inf
    for trial in range(20000):
        k = rng.randint(2, 40)
        g0 = np.random.default_rng(trial)
        P = g0.dirichlet(np.ones(k) * rng.choice([0.2, 1, 5]))
        f = P * np.exp(g0.normal(0, rng.choice([0.1, 1, 2]), k)) + g0.normal(0, rng.choice([0, 0.05, 0.3]), k) * P.max()
        lam = f.sum()
        if lam <= 0:
            continue
        eps = 0.5 * np.abs(P - f / lam).sum()
        if eps >= 1:
            continue
        g = np.maximum(f, 0) + np.abs(g0.normal(0, rng.choice([0, 0.1, 1]), k)) * P.max()
        bc = np.sqrt(P * g / g.sum()).sum()
        rhs = (1 - eps) * np.sqrt(lam / g.sum())
        worst = min(worst, bc - rhs)
        check(bc >= rhs - 1e-12, "S6 overlap inequality")
    print(f"S6 overlap step of Theorem thm:signed-cancel: min(BC - (1 - eps) sqrt(lambda/Z)) = {worst:.2e} (>= 0)")
    E = [(0, 1), (1, 2), (0, 2)]
    st = count_trees(3, E, [0, 1, 2])
    check(st == 3 and not is_spanning_tree(3, E, [0, 1, 2]), "S6 triangle")
    print(f"S6 triangle: ST(1,1,1) = {st}, indicator of spanning trees at (1,1,1) = 0")


def main():
    rng = random.Random(20261104)
    s1_placements(rng)
    s2_digon(rng)
    s3_minors(rng)
    s4_frontier(rng)
    s5_refine(rng)
    s6_sign_split(rng)
    if FAILS:
        print("FAILED:", FAILS)
        sys.exit(1)
    print("check_signed_restricted: all checks passed")


if __name__ == "__main__":
    main()

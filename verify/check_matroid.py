# /// script
# dependencies = ["numpy"]
# ///
"""Exact checks for Section "Matroids: circuits and cocircuits" (sec-matroid.tex) and the matroid
part of Section "Limits". Everything is computed from the explicit list of bases; probabilities are
exact counts, evaluated in floating point (Fractions where the paper quotes an exact value).

  M1  Proposition prop:mat-dual   beta(U(M);X|Y) = beta(U(M*);X|Y) on random splits
  M2  Lemma lem:mat-exchange      lambda(e)|D| >= 1 for every cocircuit and random cuts D containing e (M and M*)
  M3  Lemma lem:mat-frame         random framed circuits/cocircuits, S strictly inside K, random F avoiding the frames
  M4  Theorem thm:mat-local       beta <= ||sqrt law(bits)|| <= prod gamma_k <= prod sqrt((1+sqrt(1-4w))/2)
  M5  Lemma lem:mat-global        exhaustive |F|,|S| <= 2: no violation on balanced matroids, violations on S8
  M6  S8 is not balanced          the positively correlated pair (exact counts)
  M7  Lemma lem:mat-resist        1 - lambda(e; M\\Z) = 1/(1+R_Z(e)) for TU representations (R10 and M(K4))
  M8  Theorem thm:mat-global      chain <= global bound on regular matroids
  M9  U_{k,2k}                    beta^2 = max_j Pr[|B cap X| = j] >= 1/(k+1)
  M10 Fano-decorated K3           (F1)-(F3), gadget lemma, local and lifted frames, Theorem thm:mat-local

Exits 1 on any failure.
"""
import itertools
import math
import random
import sys
from fractions import Fraction

import numpy as np

TOL = 1e-12
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL:", msg)


def pc(x):
    return bin(x).count("1")


def elems(S):
    out, i = [], 0
    while S:
        if S & 1:
            out.append(i)
        S >>= 1
        i += 1
    return out


def gf2_rank(vecs):
    piv, r = {}, 0
    for v in vecs:
        x = v
        while x:
            h = x.bit_length() - 1
            if h in piv:
                x ^= piv[h]
            else:
                piv[h] = x
                r += 1
                break
    return r


class Matroid:
    def __init__(self, n, bases, name, rankf=None):
        self.n, self.name = n, name
        self.bases = sorted(set(bases))
        self.full = (1 << n) - 1
        self.r = pc(self.bases[0])
        assert all(pc(b) == self.r for b in self.bases)
        self._rankf = rankf
        self._circ = None

    def rank(self, S):
        if self._rankf is not None:
            return self._rankf(S)
        return max(pc(b & S) for b in self.bases)

    def corank(self, S):  # rank in the dual
        return pc(S) + self.rank(self.full & ~S) - self.r

    def cl(self, S):
        r, out = self.rank(S), S
        for e in range(self.n):
            if not S >> e & 1 and self.rank(S | 1 << e) == r:
                out |= 1 << e
        return out

    def dual(self):
        return Matroid(self.n, [self.full ^ b for b in self.bases], self.name + "*")

    def is_circuit(self, K):
        k = pc(K)
        return k >= 1 and self.rank(K) == k - 1 and all(self.rank(K & ~(1 << e)) == k - 1 for e in elems(K))

    def is_cocircuit(self, K):
        k = pc(K)
        return k >= 1 and self.corank(K) == k - 1 and all(self.corank(K & ~(1 << e)) == k - 1 for e in elems(K))

    def is_cut(self, D):
        Z = self.full & ~D
        return self.cl(Z) == Z

    def circuits(self):
        if self._circ is None:
            ind = set()
            for b in self.bases:
                s = b
                while True:
                    ind.add(s)
                    if s == 0:
                        break
                    s = (s - 1) & b
            self._circ = [S for S in range(1, 1 << self.n)
                          if S not in ind and all((S & ~(1 << e)) in ind for e in elems(S))]
        return self._circ

    def contract(self, F):
        rF = self.rank(F)
        return Matroid(self.n, [b & ~F for b in self.bases if pc(b & F) == rF], self.name + "/F")

    def delete(self, F):
        rE = self.rank(self.full & ~F)
        return Matroid(self.n, [b & ~F for b in self.bases if pc(b & ~F) == rE], self.name + "\\F")

    def lev(self, e):
        return sum(b >> e & 1 for b in self.bases) / len(self.bases)


def binary_matroid(cols, name):
    n, R = len(cols), gf2_rank(cols)
    bases = [sum(1 << i for i in T) for T in itertools.combinations(range(n), R)
             if gf2_rank([cols[i] for i in T]) == R]
    return Matroid(n, bases, name, rankf=lambda S: gf2_rank([cols[i] for i in elems(S)]))


def standard_form(rows, name):
    r = len(rows)
    cols = [1 << i for i in range(r)]
    for j in range(len(rows[0])):
        cols.append(sum(1 << i for i in range(r) if rows[i][j] == "1"))
    return binary_matroid(cols, name)


def graphic(edges, name):
    return binary_matroid([(1 << u) | (1 << v) for u, v in edges], name)


def uniform(k, n):
    return Matroid(n, [sum(1 << i for i in T) for T in itertools.combinations(range(n), k)], f"U{k},{n}")


K4_EDGES = list(itertools.combinations(range(4), 2))
K33_EDGES = [(u, v) for u in range(3) for v in range(3, 6)]
K5_EDGES = list(itertools.combinations(range(5), 2))
R10_ROWS = ["11001", "11100", "01110", "00111", "10011"]
S8_ROWS = ["0111", "1011", "1101", "1111"]   # [I_4 | .] as printed by Huh, Schroeter and Wang


def library():
    L = {}
    L["M(K4)"] = graphic(K4_EDGES, "M(K4)")
    L["M(K33)"] = graphic(K33_EDGES, "M(K33)")
    L["M(K5)"] = graphic(K5_EDGES, "M(K5)")
    L["M*(K33)"] = L["M(K33)"].dual()
    L["R10"] = standard_form(R10_ROWS, "R10")
    L["F7"] = standard_form(["1101", "1011", "0111"], "F7")
    L["F7*"] = L["F7"].dual()
    L["AG(3,2)"] = standard_form(["0111", "1011", "1101", "1110"], "AG(3,2)")
    L["S8"] = standard_form(S8_ROWS, "S8")
    L["U2,4"] = uniform(2, 4)
    for k, M in L.items():
        M.name = k
    return L


# ------------------------------------------------------------------ affinity, bits, chain
def opnorm_sqrt(counts, total):
    rows, cols = {}, {}
    for (a, c) in counts:
        rows.setdefault(a, len(rows))
        cols.setdefault(c, len(cols))
    A = np.zeros((len(rows), len(cols)))
    for (a, c), m in counts.items():
        A[rows[a], cols[c]] = m / total
    return float(np.linalg.norm(np.sqrt(A), 2))


def beta_exact(M, X):
    Y = M.full & ~X
    cnt = {}
    for b in M.bases:
        key = (b & X, b & Y)
        cnt[key] = cnt.get(key, 0) + 1
    return opnorm_sqrt(cnt, len(M.bases))


def piece_bits(b, K, kind, X, Y):
    KX, KY = K & X, K & Y
    if kind == "c":
        return int((b & KX) == KX), int((b & KY) == KY)
    return int((b & KX) == 0), int((b & KY) == 0)


def bits_and_chain(M, X, pieces):
    Y = M.full & ~X
    recs = [[piece_bits(b, K, kind, X, Y) for (K, kind) in pieces] for b in M.bases]
    cnt = {}
    for rec in recs:
        key = (tuple(u for u, _ in rec), tuple(v for _, v in rec))
        cnt[key] = cnt.get(key, 0) + 1
    bits = opnorm_sqrt(cnt, len(M.bases))
    gammas, eleven = [], 0
    for k in range(len(pieces)):
        groups = {}
        for rec in recs:
            groups.setdefault(tuple(rec[:k]), []).append(rec[k])
        g = 0.0
        for lst in groups.values():
            Kk = np.zeros((2, 2))
            for (x, y) in lst:
                Kk[x, y] += 1
            eleven += int(Kk[1, 1] > 0)
            g = max(g, float(np.linalg.norm(np.sqrt(Kk / Kk.sum()), 2)))
        gammas.append(g)
    return bits, gammas, eleven


def balanced_split(n, rng):
    while True:
        X = sum(1 << e for e in range(n) if rng.random() < 0.5)
        if 3 * pc(X) >= n and 3 * (n - pc(X)) >= n:
            return X


# ------------------------------------------------------------------ frames
def canonical_frame(M, K, order, kind):
    """Frame of Definition def:frame: D_e = E - cl(K - {e, pi(e)}), in M (circuit) or M* (cocircuit)."""
    k = len(order)
    D = {}
    for t, e in enumerate(order):
        f = order[(t + 1) % k]
        rest = K & ~(1 << e) & ~(1 << f)
        if kind == "c":
            D[e] = M.full & ~M.cl(rest)
        else:   # closure in the dual
            r = M.corank(rest)
            Z = rest
            for x in range(M.n):
                if not rest >> x & 1 and M.corank(rest | 1 << x) == r:
                    Z |= 1 << x
            D[e] = M.full & ~Z
    return D


def is_cycle(M, D):  # complement of a flat of M*
    Z = M.full & ~D
    r = M.corank(Z)
    return all(M.corank(Z | 1 << x) > r for x in elems(D))


def frame_ok(M, K, order, D, kind):
    k = len(order)
    for t, e in enumerate(order):
        f = order[(t + 1) % k]
        if (D[e] & K) != ((1 << e) | (1 << f)):
            return False
        if kind == "c" and not M.is_cut(D[e]):
            return False
        if kind == "d" and not is_cycle(M, D[e]):
            return False
    return True


# ------------------------------------------------------------------ checks
def m1_duality(L, rng):
    worst = 0.0
    for M in L.values():
        Md = M.dual()
        for _ in range(6):
            X = balanced_split(M.n, rng)
            worst = max(worst, abs(beta_exact(M, X) - beta_exact(Md, X)))
    print(f"M1 duality: max |beta(M)-beta(M*)| over {6 * len(L)} splits = {worst:.1e}")
    check(worst < 1e-12, "M1")


def m2_exchange(L, rng):
    worst, cases = float("inf"), 0
    for M0 in L.values():
        for M in (M0, M0.dual()):
            Ds = [D for D in M.dual().circuits()]               # cocircuits
            for _ in range(40):                                 # random cuts: complements of closures
                S = sum(1 << e for e in range(M.n) if rng.random() < 0.3)
                Z = M.cl(S)
                if Z != M.full:
                    Ds.append(M.full & ~Z)
            for D in Ds:
                for e in elems(D):
                    worst = min(worst, M.lev(e) * pc(D))
                    cases += 1
    print(f"M2 exchange bound: min lambda(e)*|D| over {cases} (cut, element) pairs = {worst:.4f}  (need >= 1)")
    check(worst >= 1 - TOL, "M2")


def m3_frame(L, rng):
    viol, cases, minratio = 0, 0, float("inf")
    for M in L.values():
        circs = [K for K in M.circuits() if pc(K) >= 2]
        cocs = [K for K in M.dual().circuits() if pc(K) >= 2]
        for _ in range(120):
            kind = rng.choice("cd")
            K = rng.choice(circs if kind == "c" else cocs)
            order = elems(K)
            rng.shuffle(order)
            D = canonical_frame(M, K, order, kind)
            check(frame_ok(M, K, order, D, kind), f"M3 frame construction {M.name}")
            sub = [e for e in order if rng.random() < 0.6]
            if len(sub) == len(order):
                sub = sub[:-1]
            if not sub:
                sub = order[:1]
            S = sum(1 << e for e in sub)
            forbid = K
            for e in sub:
                forbid |= D[e]
            F = sum(1 << x for x in range(M.n) if not forbid >> x & 1 and rng.random() < 0.6)
            bound = math.prod(1 / pc(D[e]) for e in sub)
            groups = {}
            for b in M.bases:
                g = groups.setdefault(b & F, [0, 0])
                g[0] += 1
                hit = (b & S) == S if kind == "c" else (b & S) == 0
                g[1] += int(hit)
            for tot, hit in groups.values():
                cases += 1
                minratio = min(minratio, hit / tot / bound)
                if hit / tot < bound - TOL:
                    viol += 1
    print(f"M3 frame bound: {viol} violations in {cases} (piece, S, F, sigma) cases; min Pr/bound = {minratio:.4f}")
    check(viol == 0, "M3")


def compatible(p, q):
    (K1, _, T1), (K2, _, T2) = p, q
    return not (K1 & T2) and not (K2 & T1)


def m4_local(L, rng):
    n_ok, n_split = 0, 0
    for M in L.values():
        pieces = []
        for kind, fam in (("c", M.circuits()), ("d", M.dual().circuits())):
            for K in fam:
                if pc(K) < 2:
                    continue
                order = elems(K)
                rng.shuffle(order)
                D = canonical_frame(M, K, order, kind)
                T = 0
                for e in order:
                    T |= D[e]
                w = math.prod(1 / pc(D[e]) for e in order)
                pieces.append((K, kind, T, w))
        pieces.sort(key=lambda p: -p[3])
        for _ in range(8):
            X = balanced_split(M.n, rng)
            Y = M.full & ~X
            chosen = []
            for p in pieces:
                if (p[0] & X) and (p[0] & Y) and all(compatible((p[0], p[1], p[2]), (q[0], q[1], q[2])) for q in chosen):
                    chosen.append(p)
            if not chosen:
                continue
            n_split += 1
            beta = beta_exact(M, X)
            bits, gam, eleven = bits_and_chain(M, X, [(p[0], p[1]) for p in chosen])
            chain = math.prod(gam)
            bound = math.prod(math.sqrt((1 + math.sqrt(1 - 4 * p[3])) / 2) for p in chosen)
            ok = eleven == 0 and beta <= bits + TOL and bits <= chain + TOL and chain <= bound + TOL
            n_ok += ok
            check(ok, f"M4 {M.name}")
    print(f"M4 Theorem thm:mat-local: {n_ok}/{n_split} (matroid, split) pairs satisfy beta <= bits <= prod gamma <= bound,"
          f" with no (1,1) cell")


def lemma_global_exhaustive(M, maxF=2, maxS=2):
    viol, cases, worst = 0, 0, Fraction(1)
    for kF in range(1, maxF + 1):
        for Ft in itertools.combinations(range(M.n), kF):
            F = sum(1 << i for i in Ft)
            groups = {}
            for b in M.bases:
                groups.setdefault(b & F, []).append(b)
            MF, MdF = M.contract(F), M.delete(F)
            rest = [i for i in range(M.n) if not F >> i & 1]
            for kS in range(1, maxS + 1):
                for St in itertools.combinations(rest, kS):
                    S = sum(1 << i for i in St)
                    for mode in ("in", "out"):
                        if mode == "in":
                            m = min(Fraction(sum((b & S) == S for b in g), len(g)) for g in groups.values())
                            rhs = Fraction(sum((b & S) == S for b in MF.bases), len(MF.bases))
                        else:
                            m = min(Fraction(sum((b & S) == 0 for b in g), len(g)) for g in groups.values())
                            rhs = Fraction(sum((b & S) == 0 for b in MdF.bases), len(MdF.bases))
                        cases += 1
                        if m < rhs:
                            viol += 1
                            worst = min(worst, m / rhs)
    return viol, cases, worst


def m5_global_lemma(L):
    for name in ("M(K4)", "M(K33)", "M(K5)", "R10", "F7", "F7*", "AG(3,2)", "U2,4", "S8"):
        v, c, w = lemma_global_exhaustive(L[name])
        print(f"M5 Lemma lem:mat-global, |F|,|S|<=2, {name:8s}: {v} violations / {c} cases, worst ratio {w}")
        if name == "S8":
            check(v > 0 and w < 1, "M5 S8 should violate")
        else:
            check(v == 0, f"M5 {name}")


def m6_s8(L):
    M = L["S8"]
    N = len(M.bases)
    pairs = []
    for e, f in itertools.combinations(range(M.n), 2):
        ne = sum(b >> e & 1 for b in M.bases)
        nf = sum(b >> f & 1 for b in M.bases)
        nef = sum((b >> e & 1) & (b >> f & 1) for b in M.bases)
        if nef * N > ne * nf:
            pairs.append((e, f, ne, nf, nef))
    print(f"M6 S8 ({'/'.join(S8_ROWS)} after I_4): {N} bases; positively correlated pairs (e,f,N_e,N_f,N_ef): {pairs}")
    for e, f, ne, nf, nef in pairs:
        print(f"   N_ef * N = {nef * N} > N_e * N_f = {ne * nf}")
    check(N == 48 and len(pairs) == 1 and pairs[0][2:] == (28, 20, 12), "M6")


def is_tu(A):
    r, c = A.shape
    for k in range(1, min(r, c) + 1):
        for R in itertools.combinations(range(r), k):
            sub = A[list(R)]
            for C in itertools.combinations(range(c), k):
                if abs(round(np.linalg.det(sub[:, list(C)]))) > 1:
                    return False
    return True


def tu_signing(rows):
    """Signs on a spanning forest of the support can be fixed by scaling rows and columns,
    so only the remaining entries are searched."""
    A = np.array([[int(ch) for ch in row] for row in rows], dtype=float)
    r = A.shape[0]
    parent = list(range(r + A.shape[1]))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    free = []
    for i in range(r):
        for j in range(A.shape[1]):
            if A[i, j]:
                a, b = find(i), find(r + j)
                if a == b:
                    free.append((i, j))
                else:
                    parent[a] = b
    for signs in itertools.product([1, -1], repeat=len(free)):
        S = A.copy()
        for (i, j), s in zip(free, signs):
            S[i, j] = s
        if is_tu(S):
            return S
    return None


def m7_resistance(L, rng):
    S = tu_signing(R10_ROWS)
    check(S is not None, "M7 no TU signing of R10")
    A_R10 = np.hstack([np.eye(5), S])
    inc = np.zeros((4, len(K4_EDGES)))
    for j, (u, v) in enumerate(K4_EDGES):
        inc[u, j], inc[v, j] = 1, -1
    for name, A, M in (("R10", A_R10, L["R10"]), ("M(K4)", inc[1:], L["M(K4)"])):
        worst, cases, coloops = 0.0, 0, 0
        for _ in range(80):
            Z = [i for i in range(M.n) if rng.random() < 0.35]
            MZ = M.delete(sum(1 << i for i in Z))
            for e in range(M.n):
                if e in Z:
                    continue
                keep = [f for f in range(M.n) if f not in Z and f != e]
                x = np.linalg.lstsq(A[:, keep], A[:, e], rcond=None)[0] if keep else np.zeros(0)
                feas = np.linalg.norm((A[:, keep] @ x if keep else 0) - A[:, e]) < 1e-9
                pred = 1.0 / (1.0 + float(x @ x)) if feas else 0.0
                coloops += int(not feas)
                worst = max(worst, abs(pred - (1 - MZ.lev(e))))
                cases += 1
        print(f"M7 Lemma lem:mat-resist on {name}: max |1/(1+R_Z(e)) - (1-lambda(e;M\\Z))| = {worst:.1e}"
              f" over {cases} (Z,e), {coloops} coloop cases")
        check(worst < 1e-9, f"M7 {name}")


def m8_global(L, rng):
    n_ok, n_split = 0, 0
    for name in ("M(K4)", "M(K33)", "M(K5)", "M*(K33)", "R10"):
        M = L[name]
        fam = [(K, "c") for K in M.circuits()] + [(K, "d") for K in M.dual().circuits()]
        for _ in range(6):
            X = balanced_split(M.n, rng)
            Y = M.full & ~X
            rng.shuffle(fam)
            chosen, used = [], 0
            for K, kind in fam:
                if (K & X) and (K & Y) and not (K & used):
                    chosen.append((K, kind))
                    used |= K
            if not chosen:
                continue
            n_split += 1
            bound, F = 1.0, 0
            for K, kind in chosen:
                phi = 1.0
                for H in (K & X, K & Y):
                    prev = 0
                    for e in elems(H):
                        if kind == "c":
                            phi *= M.contract(F | prev).lev(e)
                        else:
                            phi *= 1 - M.delete(F | prev).lev(e)
                        prev |= 1 << e
                bound *= math.sqrt((1 + math.sqrt(max(0.0, 1 - 4 * phi))) / 2)
                F |= K
            beta = beta_exact(M, X)
            bits, gam, eleven = bits_and_chain(M, X, chosen)
            chain = math.prod(gam)
            ok = eleven == 0 and beta <= bits + TOL and bits <= chain + TOL and chain <= bound + TOL
            n_ok += ok
            check(ok, f"M8 {name}")
    print(f"M8 Theorem thm:mat-global: {n_ok}/{n_split} (regular matroid, split) pairs satisfy beta <= bits <= prod gamma <= bound")


def m9_uniform():
    for k in range(3, 8):
        M = uniform(k, 2 * k)
        for a in (k, k - 1):
            X = (1 << a) - 1
            b2 = beta_exact(M, X) ** 2
            pj = max(math.comb(a, j) * math.comb(2 * k - a, k - j) for j in range(k + 1)) / math.comb(2 * k, k)
            print(f"M9 U_{k},{2 * k} |X|={a}: beta^2 = {b2:.5f}, max_j Pr[|B cap X|=j] = {pj:.5f}, 1/(k+1) = {1 / (k + 1):.5f}")
            check(abs(b2 - pj) < 1e-9 and b2 >= 1 / (k + 1) - TOL, f"M9 U{k},{2 * k}")


# ------------------------------------------------------------------ Fano-decorated K3
P = 7
LINES = [L for L in itertools.combinations(range(1, 8), 3) if L[0] ^ L[1] ^ L[2] == 0]
LINES_P = [L for L in LINES if P in L]
LINES_NP = [L for L in LINES if P not in L]


def fano_decorated(edges, nV):
    """2-sum of M(G) with a Fano plane at every edge (base point p = 7 identified with the edge).
    Element 6*i + (q-1) is the point q of the gadget of edge i, q = 1..6."""
    cols, glue = [], []
    for i, (u, v) in enumerate(edges):
        base = nV + 3 * i
        for q in range(1, 7):
            cols.append(q << base)
        glue.append((1 << u) | (1 << v) | (P << base))
    rZ = gf2_rank(glue)

    def rank(S):
        return gf2_rank([cols[j] for j in elems(S)] + glue) - rZ

    n = len(cols)
    R = rank((1 << n) - 1)
    bases = [sum(1 << j for j in T) for T in itertools.combinations(range(n), R) if rank(sum(1 << j for j in T)) == R]
    return Matroid(n, bases, "Fano-decorated K3", rankf=rank)


def m10_fano(rng):
    edges = [(0, 1), (1, 2), (0, 2)]
    M = fano_decorated(edges, 3)
    el = lambda i, q: 6 * i + q - 1
    pts = lambda i, Q: sum(1 << el(i, q) for q in Q if q != P)
    print(f"M10 {M.name}: n = {M.n}, rank = {M.r}, bases = {len(M.bases)}")
    ok = True
    # (F1) circuits and cocircuits of a gadget avoiding p
    for i in range(3):
        ok &= all(M.is_circuit(pts(i, L)) for L in LINES_NP)
        ok &= all(M.is_cocircuit(pts(i, set(range(1, 8)) - set(L))) for L in LINES_P)
    check(ok, "M10 (F1)")
    # (F2) lifted triangles, all 27 choices of lines through p
    lifts = [sum(pts(i, Ls[i]) for i in range(3)) for Ls in itertools.product(LINES_P, repeat=3)]
    f2 = all(M.is_circuit(K) for K in lifts)
    check(f2, "M10 (F2)")
    # (F3) lifted stars at each vertex, all choices of cocircuits containing p
    cocp = [set(range(1, 8)) - set(L) for L in LINES_NP]
    f3 = True
    for v in range(3):
        star = [i for i, e in enumerate(edges) if v in e]
        for Ds in itertools.product(cocp, repeat=len(star)):
            f3 &= M.is_cocircuit(sum(pts(i, D) for i, D in zip(star, Ds)))
    check(f3, "M10 (F3)")
    # gadget lemma: every nontrivial 2-colouring of the six points has a mixed p-avoiding line
    gad = all(any(len({c >> (q - 1) & 1 for q in L}) == 2 for L in LINES_NP) for c in range(1, 63))
    check(gad, "M10 gadget lemma")
    # local pieces with their frames
    local = []
    for i in range(3):
        for L in LINES_NP:
            K = pts(i, L)
            D = {}
            for t in range(3):
                x, y, z = L[t], L[(t + 1) % 3], L[(t + 2) % 3]
                Lz = next(l for l in LINES_P if z in l)
                D[el(i, x)] = pts(i, set(range(1, 8)) - set(Lz))
            order = [el(i, q) for q in L]
            ok &= frame_ok(M, K, order, D, "c")
            T = 0
            for d in D.values():
                T |= d
            local.append((K, T, math.prod(1 / pc(d) for d in D.values()), i))
    check(ok, "M10 local frames")
    # lifted triangle with frames (one choice of lines)
    ab = {i: [q for q in LINES_P[i] if q != P] for i in range(3)}
    order = []
    for i in range(3):
        order += [el(i, ab[i][0]), el(i, ab[i][1])]
    K = sum(1 << e for e in order)
    D = {}
    for i in range(3):
        a, b = ab[i]
        Lp = next(l for l in LINES_P if a not in l and b not in l)
        D[el(i, a)] = pts(i, set(range(1, 8)) - set(Lp))
        j = (i + 1) % 3
        v = (set(edges[i]) & set(edges[j])).pop()
        cut = 0
        for f, e in enumerate(edges):
            if v not in e:
                continue
            want_in, want_out = (ab[i][1], ab[i][0]) if f == i else (ab[j][0], ab[j][1])
            # complement of a line through want_out avoiding p: contains p and want_in
            Lc = next(l for l in LINES_NP if want_out in l)
            comp = set(range(1, 8)) - set(Lc)
            assert want_in in comp and P in comp
            cut |= pts(f, comp)
        D[el(i, b)] = cut
    lifted_ok = M.is_circuit(K) and frame_ok(M, K, order, D, "c")
    T = 0
    for d in D.values():
        T |= d
    wl = math.prod(1 / pc(d) for d in D.values())
    check(lifted_ok, "M10 lifted frames")
    print(f"M10 (F1) {ok}, (F2) {f2}, (F3) {f3}, gadget lemma {gad}; local pieces framed, weight 4^-3;"
          f" lifted triangle framed: {lifted_ok}, weight {wl:.3e} >= (12*2)^-3 = {24 ** -3:.3e}")
    check(wl >= 24 ** -3 - TOL, "M10 lifted weight")
    # Theorem thm:mat-local on random splits and on gadget-whole splits
    splits = [balanced_split(M.n, rng) for _ in range(4)]
    splits += [pts(0, range(1, 8)) | pts(1, range(1, 8)), pts(2, range(1, 8))]
    for X in splits:
        Y = M.full & ~X
        fams = []
        loc = []
        for (Kp, Tp, w, i) in local:
            if (Kp & X) and (Kp & Y) and all(i != q[3] for q in loc):
                loc.append((Kp, Tp, w, i))
        if loc:
            fams.append(("local", [(p[0], "c") for p in loc], [p[2] for p in loc]))
        if (K & X) and (K & Y):
            fams.append(("lifted", [(K, "c")], [wl]))
        beta = beta_exact(M, X)
        for nm, pcs, ws in fams:
            bits, gam, eleven = bits_and_chain(M, X, pcs)
            chain = math.prod(gam)
            bound = math.prod(math.sqrt((1 + math.sqrt(1 - 4 * w)) / 2) for w in ws)
            good = eleven == 0 and beta <= bits + TOL and bits <= chain + TOL and chain <= bound + TOL
            print(f"   |X|={pc(X):2d} {nm:6s} pieces={len(pcs)}: beta^2={beta ** 2:.4f} chain^2={chain ** 2:.4f} bound^2={bound ** 2:.4f}  {'OK' if good else 'FAIL'}")
            check(good, f"M10 split {nm}")


def main():
    rng = random.Random(20261104)
    L = library()
    print("library:", ", ".join(f"{k} (n={M.n}, {len(M.bases)} bases)" for k, M in L.items()))
    m1_duality(L, rng)
    m2_exchange(L, rng)
    m3_frame(L, rng)
    m4_local(L, rng)
    m5_global_lemma(L)
    m6_s8(L)
    m7_resistance(L, rng)
    m8_global(L, rng)
    m9_uniform()
    m10_fano(rng)
    if FAILS:
        print("FAILED:", FAILS)
        sys.exit(1)
    print("check_matroid: all checks passed")


if __name__ == "__main__":
    main()

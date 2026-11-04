# /// script
# dependencies = ["numpy", "scipy", "networkx", "mpmath"]
# ///
"""Checks for Section "The complete graph and dense graphs" (sec-complete.tex).

  K1  constants of Theorem thm:complete in interval arithmetic: T(rho) = 0.384, Lambda_inf, K, R > 0.0069032 > 1/146,
      C0 < 2.01, |Lambda'| <= 1.094; the ceilings of Remark rem:kn-ceiling (0.00703, 0.1081) and the quasi-star
      exponent (about 0.0135)
  K2  the bound Phi_n of Theorem thm:complete: exact Phi_n(r) (high precision) for 3 <= n <= 300 and the
      Cauchy bound used in the proof for 3 <= n <= 10^6, against 0.0069032 n - ln(n)/2 - 2.01 and n/146 - ln(n)/2 - 3
  K3  Lemma lem:kn-iso (I1)-(I3) exactly at n = 7, 8 (all trees), and Monte Carlo of (I2) at n = 300
  K4  Lemma lem:kn-star is sharp on K_n: smallest ratio Pr[a_w = 1 | z, sigma] / p_A(w) over all pasts, n = 6, 7, 8
  K5  the chain of Lemma lem:kn-master for 6 <= n <= 9 (all trees), with the range of -ln(beta)/n
  K6  Proposition prop:two-cliques for h = 3, 4, 5 (all trees)
  K7  Lemma lem:forest-mgf on small graphs (all trees), and the arithmetic of Corollary cor:dense
  K8  the product-of-marginals estimate of Remark rem:kn-ceiling(c): exact for n = 6, 7, Monte Carlo at n = 40, 160
  K9  Proposition prop:complete-upper: the subset recursion generates every labelled tree once (n <= 7),
      counts n^{n-2} (n <= 12), and its size

Exits 1 on any failure.
"""
import heapq
import itertools
import math
import sys
from functools import lru_cache

import mpmath as mp
import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import gammaln, lambertw

from common import enum_trees

FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL:", msg, flush=True)


def Lam(p):
    return (1 + math.sqrt(max(0.0, 1 - 4 * p))) / 2


def cm(m):
    return (1 + 1 / m) ** (-m)


# ------------------------------------------------------------------ K1 constants
def k1_constants():
    iv = mp.iv
    iv.dps = 50
    alpha = iv.mpf("0.38")
    t = iv.mpf("0.384")
    rho = t * iv.exp(-t)                       # T(rho) = 0.384 exactly, since T(t e^-t) = t for t <= 1
    p_inf = 1 / (9 * iv.exp(1))
    lam_inf = (1 + iv.sqrt(1 - 4 * p_inf)) / 2
    K = t - rho + lam_inf * rho
    R = alpha * (1 + iv.log(rho / alpha)) - (1 - alpha) * K
    C0 = (1 + iv.log(rho / alpha)) + K - iv.log(iv.mpf("0.62")) + iv.log(2 * iv.pi * alpha) / 2 \
        + iv.mpf(1) / 12 + iv.mpf("1.094") * rho / (9 * iv.exp(1))
    dlam = 1 / iv.sqrt(1 - 4 * p_inf)          # |Lambda'| is increasing in p, so this is the max on p <= 1/(9e)
    print(f"K1 (interval arithmetic) rho = {float(rho.a):.10f}, Lambda_inf = {float(lam_inf.a):.8f}, "
          f"K = {float(K.a):.8f}, R in [{float(R.a):.10f}, {float(R.b):.10f}] (1/R = {float((1 / R).a):.4f}), "
          f"C0 <= {float(C0.b):.6f}, max |Lambda'| <= {float(dlam.b):.6f}")
    check(rho.b < 1 / iv.exp(1).b, "K1 rho < 1/e")
    check(R.a > iv.mpf("0.0069032").b and iv.mpf("0.0069032").a > (iv.mpf(1) / 146).b, "K1 R > 0.0069032 > 1/146")
    check(C0.b < 2.01, "K1 C0 < 2.01")
    check(dlam.b <= 1.094, "K1 |Lambda'| <= 1.094")
    check(abs(float(lam_inf.a) - 0.957301) < 1e-6 and abs(float(K.a) - 0.37283) < 1e-5, "K1 printed values")

    # Remark rem:kn-ceiling (a), (b)
    g = lambda p: -math.log(Lam(p))
    a_star = (3 - math.sqrt(5)) / 2
    f_star = a_star * (1 - a_star) * math.exp(-a_star)
    grid = max(a * (1 - a) * math.exp(-a) for a in np.linspace(0, 1, 200001))
    ceil_a = f_star * g(1 / (9 * math.e))
    print(f"K1 max alpha(1-alpha)e^-alpha = {f_star:.5f} at (3-sqrt5)/2 (grid max {grid:.5f}), "
          f"g(1/(9e)) = {g(1 / (9 * math.e)):.5f}, ceiling = {ceil_a:.6f}, R/ceiling = {float(R.a) / ceil_a:.4f}")
    check(abs(grid - f_star) < 1e-9 and abs(f_star - 0.16112) < 5e-6 and abs(g(1 / (9 * math.e)) - 0.04364) < 5e-6,
          "K1 ceiling (a) values")
    check(0.00702 < ceil_a < 0.00704 and 0.98 <= float(R.a) / ceil_a < 0.99, "K1 ceiling (a) 0.00703 and 98%")
    print(f"K1 g(1/(4e)) = {g(1 / (4 * math.e)):.5f}")
    check(abs(g(1 / (4 * math.e)) - 0.1081) < 5e-5, "K1 ceiling (b)")

    # quasi-star witness: same chain with the kernel p_q, optimized over rho
    s = 1 - math.sqrt(2 / 3)
    x = s / (1 - 0.38)
    pq = x * (1 - x) / math.e
    lq = Lam(pq)

    def rate(r_):
        T = float(-lambertw(-r_, 0).real)
        return 0.38 * (1 + math.log(r_ / 0.38)) - 0.62 * (T - r_ + lq * r_)
    rq = max(rate(r_) for r_ in np.linspace(0.05, 1 / math.e - 1e-9, 20001))
    print(f"K1 quasi-star: s = {s:.5f}, p = {pq:.6f}, exponent per vertex = {rq:.5f}")
    check(abs(pq - 0.0766562) < 1e-6 and 0.0134 < rq < 0.0137, "K1 quasi-star 0.0135")
    return float(R.a), float(C0.b)


# ------------------------------------------------------------------ K2 Phi_n
def phi_exact(n, r, lmb):
    """Phi_n(r) at Lambda_0 = lmb, by the explicit alternating sum of (I2), in enough precision."""
    m = n - r
    with mp.workdps(40 + int(r * math.log10(2 / max(lmb, 1e-3)))):
        s = mp.mpf(0)
        ch = 1
        L = mp.mpf(lmb)
        for k in range(r + 1):
            s += ch * (L - 1) ** k * (mp.mpf(m) / n) ** k * (1 - mp.mpf(k) / n) ** (r - k - 1)
            ch = ch * (r - k) // (k + 1)
        return s


def gf_coeff(n, r, lmb):
    """r!/(m n^(r-1)) [x^r] exp(m(T(x)-x) + m lmb x): the right-hand side of (I2)."""
    m = n - r
    with mp.workdps(40):
        a = [mp.mpf(0)] * (r + 1)
        for k in range(2, r + 1):
            a[k] = m * mp.mpf(k) ** (k - 1) / mp.factorial(k)
        a[1] = m * mp.mpf(lmb)
        e = [mp.mpf(0)] * (r + 1)
        e[0] = mp.mpf(1)
        for k in range(1, r + 1):
            e[k] = mp.fsum(j * a[j] * e[k - j] for j in range(1, k + 1)) / k
        return mp.factorial(r) / (m * mp.mpf(n) ** (r - 1)) * e[r]


def k2_phi(R, C0):
    worst = (math.inf, 0)
    worst_clean = (math.inf, 0)
    for n in range(3, 301):
        r = int(math.floor(0.38 * n))
        m = n - r
        l0 = Lam(cm(m) * (1 - 1 / m) / 9)
        v = -float(mp.log(phi_exact(n, r, l0)))
        s1 = v - (0.0069032 * n - 0.5 * math.log(n) - 2.01)
        s2 = v - (n / 146 - 0.5 * math.log(n) - 3)
        worst = min(worst, (s1, n))
        worst_clean = min(worst_clean, (s2, n))
        if n < 40:
            d = abs(gf_coeff(n, r, l0) - phi_exact(n, r, l0)) / phi_exact(n, r, l0)
            check(d < 1e-25, f"K2 (I2) identity n={n}")
    print(f"K2 exact Phi_n(floor(0.38n)), 3 <= n <= 300: smallest slack against 0.0069032n - ln(n)/2 - 2.01 is "
          f"{worst[0]:.3f} (n={worst[1]}), against n/146 - ln(n)/2 - 3 is {worst_clean[0]:.3f} (n={worst_clean[1]})")
    check(worst[0] >= 0 and worst_clean[0] >= 0, "K2 exact Phi_n")
    # the Cauchy bound of the proof at rho = 0.384 e^-0.384, for every n <= 10^6
    ns = np.arange(3, 10 ** 6 + 1, dtype=float)
    r = np.floor(0.38 * ns)
    m = ns - r
    l0 = (1 + np.sqrt(1 - 4 * (1 + 1 / m) ** (-m) * (1 - 1 / m) / 9)) / 2
    rho = 0.384 * math.exp(-0.384)
    F = -(gammaln(r + 1) - np.log(m) - (r - 1) * np.log(ns) + m * (0.384 - rho + l0 * rho) - r * np.log(rho))
    sl = F - (0.0069032 * ns - 0.5 * np.log(ns) - 2.01)
    sl2 = F - (R * ns - 0.5 * np.log(ns) - C0)
    i, j = int(np.argmin(sl)), int(np.argmin(sl2))
    print(f"K2 Cauchy bound, 3 <= n <= 10^6: smallest slack against 0.0069032n - ln(n)/2 - 2.01 is {sl[i]:.4f} "
          f"(n={int(ns[i])}), against Rn - ln(n)/2 - C0 is {sl2[j]:.4f} (n={int(ns[j])}); "
          f"-ln Phi/n at n=10^6 >= {F[-1] / 1e6:.6f}")
    check(sl[i] >= 0 and sl2[j] >= 0, "K2 Cauchy bound")
    for n in (10, 50, 200, 300):
        r = int(0.38 * n)
        m = n - r
        l0 = Lam(cm(m) * (1 - 1 / m) / 9)
        check(float(phi_exact(n, r, l0)) <= math.exp(-F[n - 3]) * (1 + 1e-12), f"K2 Phi <= Cauchy n={n}")


# ------------------------------------------------------------------ helpers on K_n
def kn_edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def bit(k):
    return np.uint64(1) << np.uint64(int(k))


def mask_of(ks):
    m = np.uint64(0)
    for k in ks:
        m |= bit(k)
    return m


def opnorm_sqrt(rows, cols, probs):
    ur, ri = np.unique(rows, return_inverse=True)
    uc, ci = np.unique(cols, return_inverse=True)
    Mp = sp.coo_matrix((probs, (ri, ci)), shape=(len(ur), len(uc))).tocsr()
    Mp.sum_duplicates()
    Ms = Mp.copy()
    Ms.data = np.sqrt(Ms.data)
    if min(Ms.shape) == 1:
        return float(np.sqrt(Mp.sum()))
    if min(Ms.shape) <= 400 and max(Ms.shape) <= 20000:
        return float(np.linalg.svd(Ms.toarray(), compute_uv=False)[0])
    return float(spla.svds(Ms, k=1, tol=1e-12, return_singular_vectors=False)[0])


def beta_exact(T, Am):
    return opnorm_sqrt(T & Am, T & ~Am, np.full(len(T), 1.0 / len(T)))


def make_split(n, E, kind, rng):
    M = len(E)
    lo = -(-M // 3)
    if kind == "random":
        A = set(rng.choice(M, size=M // 2, replace=False).tolist())
    elif kind == "clique":
        s = next(s for s in range(2, n + 1) if s * (s - 1) // 2 >= lo and M - s * (s - 1) // 2 >= lo)
        A = {k for k, (i, j) in enumerate(E) if i < s and j < s}
    elif kind == "star":      # quasi-star: B = edges touching the first s vertices, A = clique on the rest
        s = next(s for s in range(1, n + 1)
                 if lo <= sum(1 for (i, j) in E if i < s or j < s) <= M - lo)
        A = {k for k, (i, j) in enumerate(E) if not (i < s or j < s)}
    elif kind == "bip":
        s = next(s for s in range(1, n) if s * (n - s) >= lo and M - s * (n - s) >= lo)
        A = {k for k, (i, j) in enumerate(E) if (i < s) != (j < s)}
    else:
        raise ValueError(kind)
    inA = np.array([k in A for k in range(M)])
    assert 3 * inA.sum() >= M and 3 * (~inA).sum() >= M
    return inA


def star_kernels(n, E, inA, W):
    """p_A(w), p_B(w) of Lemma lem:kn-star for K_n (lambda_U = d_U = m)."""
    Wset = set(W)
    U = [u for u in range(n) if u not in Wset]
    m = len(U)
    eid = {e: k for k, e in enumerate(E)}
    out = {}
    for w in W:
        dA = sum(1 for u in U if inA[eid[(min(u, w), max(u, w))]])
        dB = m - dA
        out[w] = (dB / m * (1 + 1 / m) ** (-dA), dA / m * (1 + 1 / m) ** (-dB))
    return out


def star_masks(n, E, inA, W):
    Wset = set(W)
    eid = {e: k for k, e in enumerate(E)}
    SA, SB = {}, {}
    for w in W:
        a = b = np.uint64(0)
        for u in range(n):
            if u in Wset:
                continue
            k = eid[(min(u, w), max(u, w))]
            if inA[k]:
                a |= bit(k)
            else:
                b |= bit(k)
        SA[w], SB[w] = a, b
    return SA, SB


def touched_by(z, E, EW):
    t = set()
    for k in EW:
        if int(z) >> k & 1:
            t |= set(E[k])
    return t


# ------------------------------------------------------------------ K3 Lemma lem:kn-iso
def iso_sizes(n, E, T, W):
    Wset = set(W)
    EW = [k for k, (i, j) in enumerate(E) if i in Wset and j in Wset]
    z = T & mask_of(EW)
    uz, zi = np.unique(z, return_inverse=True)
    isos = [frozenset(w for w in W if w not in touched_by(zz, E, EW)) for zz in uz]
    return [isos[g] for g in zi]


def k3_iso(rng):
    for n, r in ((7, 3), (7, 4), (8, 4)):
        E = kn_edges(n)
        T = enum_trees(n, E)
        W = sorted(rng.choice(n, size=r, replace=False).tolist())
        I = iso_sizes(n, E, T, W)
        N = len(T)
        m = n - r
        worst = 0.0
        for s in range(r + 1):
            for S in itertools.combinations(W, s):
                S = frozenset(S)
                emp = sum(1 for x in I if S <= x) / N
                worst = max(worst, abs(emp - (m / n) ** s * (1 - s / n) ** (r - s - 1)))
        sizes = np.array([len(x) for x in I])
        w2 = 0.0
        for lmb in (0.0, 0.3, 0.7, 0.95, 1.3):
            w2 = max(w2, abs(float(np.mean(lmb ** sizes)) - float(phi_exact(n, r, lmb))))
        idx = {w: i for i, w in enumerate(W)}
        w3 = -math.inf
        for t in range(20):
            lw = rng.uniform(0, 1, r) if t % 2 else rng.uniform(0.9, 1.1, r)
            lhs = float(np.mean([np.prod([lw[idx[w]] for w in x]) for x in I]))
            rhs = float(np.mean(lw.mean() ** sizes))
            w3 = max(w3, lhs / rhs - 1)
        print(f"K3 n={n} r={r} ({N} trees): max error (I1) {worst:.1e}, (I2) {w2:.1e}; "
              f"largest (I3) excess {w3:.1e}")
        check(worst < 1e-12 and w2 < 1e-12 and w3 < 1e-12, f"K3 n={n} r={r}")
    n, r, S = 300, 114, 2000
    W = list(range(r))
    sizes = []
    for _ in range(S):
        Tt = prufer_tree(n, rng)
        touched = {x for (a, b) in Tt if a < r and b < r for x in (a, b)}
        sizes.append(r - len(touched))
    sizes = np.array(sizes)
    for lmb in (0.9, 0.97):
        v = lmb ** sizes
        mc, se = float(v.mean()), float(v.std() / math.sqrt(S))
        th = float(phi_exact(n, r, lmb))
        print(f"K3 Monte Carlo n={n} r={r} Lambda={lmb}: E Lambda^|I| = {mc:.4e} +- {se:.1e}, (I2) gives {th:.4e}")
        check(abs(mc - th) < 5 * se + 1e-15, f"K3 MC n={n} lam={lmb}")


# ------------------------------------------------------------------ K4 sharp star deletion
def k4_sharp(rng):
    overall = math.inf
    for n in (6, 7, 8):
        E = kn_edges(n)
        T = enum_trees(n, E)
        for kind in ("random", "star", "bip"):
            inA = make_split(n, E, kind, rng)
            for r in (2, 3):
                perm = rng.permutation(n)
                W = sorted(int(x) for x in perm[:r])
                Wset = set(W)
                U = [u for u in range(n) if u not in Wset]
                EW = [k for k, (i, j) in enumerate(E) if i in Wset and j in Wset]
                z = T & mask_of(EW)
                SA, SB = star_masks(n, E, inA, W)
                pk = star_kernels(n, E, inA, W)
                worst = math.inf
                for zz in np.unique(z):
                    Wp = [w for w in W if w not in touched_by(zz, E, EW)]
                    Tz = T[z == zz]
                    for w in Wp:
                        F = np.uint64(0)
                        for x in Wp:
                            if x != w:
                                F |= SA[x] | SB[x]
                        keys = Tz & F
                        uk, ki = np.unique(keys, return_inverse=True)
                        tot = np.bincount(ki)
                        for side, (S, p) in enumerate(((SA[w], pk[w][0]), (SB[w], pk[w][1]))):
                            if p == 0:
                                continue
                            hit = np.bincount(ki, weights=((Tz & S) == 0).astype(float), minlength=len(uk))
                            worst = min(worst, float((hit / tot).min() / p))
                print(f"K4 n={n} {kind:6s} W={W}: smallest Pr[a_w=1 | z, sigma] / p_A(w) = {worst:.6f}")
                check(worst >= 1 - 1e-12, f"K4 n={n} {kind} r={r}")
                overall = min(overall, worst)
    print(f"K4 smallest ratio over all instances: {overall:.6f}")
    check(overall < 1.0005, "K4 sharpness (smallest ratio is 1.000)")


# ------------------------------------------------------------------ K5 exact chain
def chain_bounds(n, E, T, inA, W):
    tau = len(T)
    Wset = set(W)
    EW = [k for k, (i, j) in enumerate(E) if i in Wset and j in Wset]
    Am = mask_of(np.where(inA)[0])
    SA, SB = star_masks(n, E, inA, W)
    pk = star_kernels(n, E, inA, W)
    lamw = {w: Lam(pk[w][0] * pk[w][1]) for w in W}
    lbar = float(np.mean(list(lamw.values())))
    z = T & mask_of(EW)
    uz, zi = np.unique(z, return_inverse=True)
    B1 = B2 = B3 = B4 = 0.0
    for g, zz in enumerate(uz):
        idx = np.where(zi == g)[0]
        Pz = len(idx) / tau
        Tz = T[idx]
        bz = beta_exact(Tz, Am)
        Wp = [w for w in W if w not in touched_by(zz, E, EW)]
        nb, b3 = 1.0, 1.0
        if Wp:
            ra = np.zeros(len(idx), np.int64)
            rb = np.zeros(len(idx), np.int64)
            for t, w in enumerate(Wp):
                ia = (Tz & SA[w]) == 0
                ib = (Tz & SB[w]) == 0
                check(not np.any(ia & ib), "K5 forbidden cell")
                ra |= ia.astype(np.int64) << t
                rb |= ib.astype(np.int64) << t
            nb = opnorm_sqrt(ra, rb, np.full(len(idx), 1.0 / len(idx)))
            for w in Wp:
                b3 *= lamw[w]
        check(bz <= nb * (1 + 1e-9) and nb ** 2 <= b3 * (1 + 1e-9), f"K5 per-z chain n={n}")
        B1 += Pz * bz ** 2
        B2 += Pz * nb ** 2
        B3 += Pz * b3
        B4 += Pz * lbar ** len(Wp)
    return math.sqrt(B1), math.sqrt(B2), math.sqrt(B3), math.sqrt(B4)


def local_search(n, E, T, rng, iters):
    M = len(E)
    inA = make_split(n, E, "random", rng)
    cur = beta_exact(T, mask_of(np.where(inA)[0]))
    for _ in range(iters):
        cand = inA.copy()
        i = int(rng.integers(M))
        cand[i] = ~cand[i]
        if 3 * cand.sum() < M or 3 * (~cand).sum() < M:
            j = int(rng.choice(np.where(cand == cand[i])[0]))
            cand[j] = ~cand[j]
            if 3 * cand.sum() < M or 3 * (~cand).sum() < M:
                continue
        b = beta_exact(T, mask_of(np.where(cand)[0]))
        if b >= cur:
            inA, cur = cand, b
    return inA


def best_W(n, E, inA, r):
    best, bv = None, -1.0
    for W in itertools.combinations(range(n), r):
        pk = star_kernels(n, E, inA, W)
        v = sum(a * b for a, b in pk.values())
        if v > bv:
            best, bv = W, v
    return list(best), bv / r


def k5_chain():
    plan = {6: (["random", "clique", "star", "bip", "worst"], (1, 2, 3), 200),
            7: (["random", "clique", "star", "bip", "worst"], (1, 2, 3), 100),
            8: (["random", "clique", "star", "bip"], (2, 3), 0),
            9: (["random", "clique", "star", "bip"], (3,), 0)}
    rates = []
    for n, (kinds, rs, iters) in plan.items():
        rng = np.random.default_rng(n)
        E = kn_edges(n)
        T = enum_trees(n, E)
        check(len(T) == n ** (n - 2), f"K5 tree count n={n}")
        for kind in kinds:
            inA = local_search(n, E, T, rng, iters) if kind == "worst" else make_split(n, E, kind, rng)
            beta = beta_exact(T, mask_of(np.where(inA)[0]))
            rates.append(-math.log(beta) / n)
            for r in rs:
                W, pmean = best_W(n, E, inA, r)
                B1, B2, B3, B4 = chain_bounds(n, E, T, inA, W)
                m = n - r
                p0 = cm(m) * (1 - 1 / m) / 9
                phi = float(phi_exact(n, r, Lam(p0)))
                ok = (beta <= B1 * (1 + 1e-9) and B1 <= B2 * (1 + 1e-9) and B2 <= B3 * (1 + 1e-9)
                      and B3 <= B4 * (1 + 1e-9) and B4 <= 1 + 1e-12)
                if pmean >= p0:
                    ok = ok and B4 ** 2 <= phi * (1 + 1e-9)
                print(f"K5 n={n} {kind:6s} |A|={int(inA.sum()):2d} r={r}: beta={beta:.4f} <= {B1:.4f} <= {B2:.4f}"
                      f" <= {B3:.5f} <= {B4:.5f}; Phi_n(r)^(1/2)={math.sqrt(phi):.5f} (mean p_w >= p_0: {pmean >= p0})"
                      f"  -ln(beta)/n={-math.log(beta) / n:.4f}  {'OK' if ok else 'FAIL'}", flush=True)
                check(ok, f"K5 chain n={n} {kind} r={r}")
    lo, hi = min(rates), max(rates)
    print(f"K5 -ln(beta)/n over all splits: min {lo:.4f}, max {hi:.4f}")
    check(0.06 <= lo and hi <= 0.15, "K5 rate range quoted in the paper")


# ------------------------------------------------------------------ K6 two cliques
def k6_two_cliques():
    for h in (3, 4, 5):
        E = ([(i, j) for i in range(h) for j in range(i + 1, h)] + [(i, h + i) for i in range(h)]
             + [(h + i, h + j) for i in range(h) for j in range(i + 1, h)])
        M = len(E)
        nA = h * (h - 1) // 2 + h
        inA = np.array([k < nA for k in range(M)])
        T = enum_trees(2 * h, E)
        tau = len(T)
        match = mask_of(range(h * (h - 1) // 2, nA))
        cnt = np.array([bin(int(x)).count("1") for x in (T & match)])
        P1 = float(np.mean(cnt == 1))
        beta = beta_exact(T, mask_of(np.where(inA)[0]))
        G = nx.Graph(E)
        lam2 = sorted(nx.laplacian_spectrum(G))[1]
        print(f"K6 h={h}: tau={tau} (formula {h ** (h - 2) * (h + 2) ** (h - 1)}), P1={P1:.6f} "
              f"(formula {(h / (h + 2)) ** (h - 1):.6f}), beta^2={beta ** 2:.6f}, e^-2={math.exp(-2):.6f}, "
              f"lambda_2={lam2:.6f}, |A|={nA}, |B|={M - nA}, M/3={M / 3:.2f}")
        check(tau == h ** (h - 2) * (h + 2) ** (h - 1) and abs(P1 - (h / (h + 2)) ** (h - 1)) < 1e-12
              and beta ** 2 >= P1 - 1e-12 and P1 > math.exp(-2) and abs(lam2 - 2) < 1e-9
              and 3 * (M - nA) >= M and 3 * nA >= M, f"K6 h={h}")


# ------------------------------------------------------------------ K7 Lemma lem:forest-mgf, Corollary cor:dense
def k7_dense(rng):
    graphs = {"K5": nx.complete_graph(5), "K33": nx.complete_bipartite_graph(3, 3),
              "Petersen": nx.petersen_graph(), "K44": nx.complete_bipartite_graph(4, 4),
              "W6": nx.wheel_graph(7)}
    worst = 0.0
    cnt = 0
    for name, G in graphs.items():
        G = nx.convert_node_labels_to_integers(G)
        E = [tuple(sorted(e)) for e in G.edges()]
        T = enum_trees(G.number_of_nodes(), E)
        inT = np.array([[(int(t) >> k) & 1 for k in range(len(E))] for t in T], bool)
        pe = inT.mean(0)
        for _ in range(30):
            F = rng.random(len(E)) < rng.uniform(0.2, 1)
            y = np.where(F, rng.uniform(1, 4, len(E)), 1.0)
            lhs = float(np.mean(np.prod(np.where(inT, y, 1.0), axis=1)))
            rhs = math.exp(float(np.sum((y - 1) * pe)))
            worst = max(worst, lhs / rhs)
            cnt += 1
    print(f"K7 Lemma lem:forest-mgf on {cnt} instances: largest ratio of the two sides {worst:.4f}")
    check(worst <= 1 + 1e-12, "K7 forest-mgf")
    l2 = math.log(2)
    ok = (l2 / 336 >= 1 / 485 and 1 / (1 - l2 / 6) <= 1.14 and 2 + 20 / 19 <= 3.5
          and all((n - math.floor(l2 / 6 * n)) * (n - math.floor(l2 / 6 * n) - 1) / ((n - 1) * (n - 2)) >= 0.5
                  for n in range(20, 5000)) and 2 * (l2 / 6) * 3 / (2 * l2) <= 0.5 + 1e-15)
    print(f"K7 Corollary cor:dense arithmetic: ln2/336 = {l2 / 336:.6f} >= 1/485 = {1 / 485:.6f}, "
          f"1/(1 - ln2/6) = {1 / (1 - l2 / 6):.4f} <= 1.14: {ok}")
    check(ok, "K7 cor:dense arithmetic")


# ------------------------------------------------------------------ K8 product of marginals
def prufer_tree(n, rng):
    seq = rng.integers(0, n, size=n - 2)
    deg = np.ones(n, int)
    for x in seq:
        deg[x] += 1
    leaves = [i for i in range(n) if deg[i] == 1]
    heapq.heapify(leaves)
    edges = []
    for x in seq:
        leaf = heapq.heappop(leaves)
        edges.append((min(leaf, int(x)), max(leaf, int(x))))
        deg[x] -= 1
        if deg[x] == 1:
            heapq.heappush(leaves, int(x))
    u, v = heapq.heappop(leaves), heapq.heappop(leaves)
    edges.append((min(u, v), max(u, v)))
    return edges


def minor_logtau(n, undecided, contracted):
    """log of the number of spanning trees of (V, undecided) / contracted."""
    g = nx.utils.UnionFind(range(n))
    for u, v in contracted:
        if g[u] == g[v]:
            return -math.inf
        g.union(u, v)
    roots = sorted({g[x] for x in range(n)})
    idx = {r: i for i, r in enumerate(roots)}
    k = len(roots)
    if k == 1:
        return 0.0
    L = np.zeros((k, k))
    for u, v in undecided:
        a, b = idx[g[u]], idx[g[v]]
        if a != b:
            L[a, a] += 1
            L[b, b] += 1
            L[a, b] -= 1
            L[b, a] -= 1
    s, ld = np.linalg.slogdet(L[1:, 1:])
    return ld if s > 0 else -math.inf


def k8_marginals(rng):
    for n in (6, 7):
        E = kn_edges(n)
        T = enum_trees(n, E)
        inA = make_split(n, E, "star", rng)
        Am = mask_of(np.where(inA)[0])
        N = len(T)
        _, ia, ca = np.unique(T & Am, return_inverse=True, return_counts=True)
        _, ib, cb = np.unique(T & ~Am, return_inverse=True, return_counts=True)
        bc = float(np.sum(np.sqrt((ca[ia] / N) * (cb[ib] / N) / N)))
        beta = beta_exact(T, Am)
        print(f"K8 n={n} quasi-star: BC(mu, mu_A x mu_B) = {bc:.5f} <= beta = {beta:.5f}; "
              f"-2 ln BC / n = {-2 * math.log(bc) / n:.4f}")
        check(bc <= beta * (1 + 1e-12), f"K8 exact n={n}")
    for n, samples in ((40, 300), (160, 300)):
        E = kn_edges(n)
        s = next(s for s in range(1, n) if (n * (n - 1) // 2 - (n - s) * (n - s - 1) // 2) * 3 >= len(E))
        inA = {e: not (e[0] < s or e[1] < s) for e in E}
        Aed = [e for e in E if inA[e]]
        Bed = [e for e in E if not inA[e]]
        lt = []
        for _ in range(samples):
            Tt = prufer_tree(n, rng)
            x = [e for e in Tt if inA[e]]
            y = [e for e in Tt if not inA[e]]
            lt.append(0.5 * (minor_logtau(n, Bed, x) + minor_logtau(n, Aed, y) - (n - 2) * math.log(n)))
        lt = np.array(lt)
        logbc = float(lt.max() + np.log(np.mean(np.exp(lt - lt.max()))))
        rate = -2 * logbc / n
        print(f"K8 n={n} quasi-star, {samples} Pruefer samples: -2 ln BC / n ~ {rate:.3f} (Monte Carlo, not certified)")
        if n == 160:
            check(0.16 <= rate <= 0.19, "K8 Monte Carlo value about 0.175 at n=160")


# ------------------------------------------------------------------ K9 exact upper bound
def gen_trees(S, v):
    S = frozenset(S)
    if len(S) == 1:
        yield frozenset()
        return
    rest = sorted(S - {v})
    s0, others = rest[0], rest[1:]
    for k in range(len(others) + 1):
        for extra in itertools.combinations(others, k):
            Sp = frozenset((s0,) + extra)
            for u in Sp:
                for t1 in gen_trees(Sp, u):
                    for t2 in gen_trees(S - Sp, v):
                        yield t1 | t2 | {tuple(sorted((u, v)))}


@lru_cache(None)
def count(S, v):
    if len(S) == 1:
        return 1
    rest = sorted(S - {v})
    s0, others = rest[0], rest[1:]
    tot = 0
    for k in range(len(others) + 1):
        for extra in itertools.combinations(others, k):
            Sp = frozenset((s0,) + extra)
            tot += sum(count(Sp, u) for u in Sp) * count(S - Sp, v)
    return tot


def circuit_size(n):
    """Nodes + wires: shared literal leaves and one constant; per state (S,v) one sum node; per child (S',u)
    a product node, its wire to the sum, wires from the two subproblems and the attachment literal, and
    one wire from the zero literal of every other edge between S' and S - S'."""
    size = n * (n - 1) + 1
    for s in range(2, n + 1):
        states = math.comb(n, s) * s
        children = sum(math.comb(s - 2, t - 1) * t for t in range(1, s))
        cut_lits = sum(math.comb(s - 2, t - 1) * t * (t * (s - t) - 1) for t in range(1, s))
        size += states * (1 + 5 * children + cut_lits)
    return size


def k9_upper():
    for n in range(3, 8):
        E = kn_edges(n)
        gen = list(gen_trees(range(n), 0))
        allT = {frozenset(E[k] for k in range(len(E)) if int(t) >> k & 1) for t in enum_trees(n, E)}
        ok = len(gen) == len(set(gen)) and set(gen) == allT
        print(f"K9 n={n}: recursion generates {len(gen)} trees, no duplicates and exactly the labelled trees: {ok}")
        check(ok, f"K9 generation n={n}")
    for n in range(2, 13):
        check(count(frozenset(range(n)), 0) == n ** (n - 2), f"K9 count n={n}")
    print("K9 counts f(V,0) = n^(n-2) for 2 <= n <= 12: " + ("OK" if not any("K9 count" in f for f in FAILS) else "FAIL"))
    for n in (8, 10, 12):
        s = circuit_size(n)
        print(f"K9 circuit size (nodes + wires) n={n}: {s:,}  (n^4 3^n = {n ** 4 * 3 ** n:,})")
        check(s <= n ** 4 * 3 ** n, f"K9 size n={n}")


def main():
    rng = np.random.default_rng(2026)
    R, C0 = k1_constants()
    k2_phi(R, C0)
    k3_iso(rng)
    k4_sharp(rng)
    k5_chain()
    k6_two_cliques()
    k7_dense(rng)
    k8_marginals(rng)
    k9_upper()
    if FAILS:
        print("FAILED:", FAILS)
        sys.exit(1)
    print("check_complete: all checks passed")


if __name__ == "__main__":
    main()

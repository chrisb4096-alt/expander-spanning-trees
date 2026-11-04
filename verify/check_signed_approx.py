# /// script
# dependencies = ["numpy", "networkx"]
# ///
"""Checks of the approximation results for circuits with real weights (Section "Does subtraction help?").
Ported from toy/r5/approx-o2-positive, toy/r5/approx-restrictions and toy/r5/approx-robust-measures.

[1] Lemma lem:kfh (Ky Fan-Hoffman) and Lemma lem:block-one (one block): on random scaled matrices of rank r
    (some with zero rows), sum_i |Q_ii| / sum_j |Q_ij| <= rank, sum_i min(s, Q_ii^+) <= sqrt(s r ||Q||_1) and
    sum_i sqrt(Q_ii^+) <= sqrt(r ||Q||_1); all sign patterns of 3x3 rank-one matrices; equality for r diagonal
    blocks J_m / (r m^2) inside an n x n zero matrix, n = r m^2.
[2] Lemma lem:fact-n on random signed Q with sum 1: overlap = 1 - eps + t and ||Q||_1 = 1 + 2t with t <= eps,
    where t is the negative mass and eps the TV distance.
[3] Pendant blocks with one-vertex gadgets (Lemmas lem:pendant-blocks, lem:pendant-prob, lem:pendant-pairs),
    exactly on small 2-connected graphs: for a balanced split and an independent set W of mixed vertices with
    G - W connected, every block has 2^g trees, its rectangle is a permutation, rectangles are disjoint, the
    conditional pendant bound (1/d)(1+kappa)^-(d-1) holds for every configuration of the other stars and every
    b in S(w), and the pair bound 2 p p' holds.  For truncated SVDs Q of M_{A,B}(P) (and Q plus a random rank-one
    term): Theorem thm:block-approx (b) unnormalized, and for Q / sum(Q) part (a), the block overlap lemma and
    Lemma lem:block-one on every block.
[4] Positive part (remark after Theorem thm:block-approx): four triangles joined by three bridges, uniform Q of
    rank 1; the bound without the positive part would claim rank >= 3.59 (>= 1.2 with the factor 1 + 2 eps),
    the stated bound holds.
[5] Proposition prop:dv: Psi(x) = prod_k (1 - prod_{a in C_k} x_a) equals tau P on every cell of every
    cycle-pair block rectangle (including the zero cells), on prisms with square and hexagon gadgets.
[6] Remark rem:padding: on K4, every assignment slice with positive UST mass has padded cancellation ratio
    1 + 2 c u(rho) / mu(rho) >= 1 + 2^(N+1) for c = 2^(M+N); zero-mass slices are skipped.

Exits nonzero if any check fails.
"""
import itertools
import math
import random
import sys
from collections import defaultdict
from fractions import Fraction

import networkx as nx
import numpy as np

from common import algebraic_connectivity, edge_list, enum_trees, laplacian, mixed_vertices, named_graph, star

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL:", msg, flush=True)
    return cond


# ----------------------------------------------------------------------------- [1] matrices
def kfh_ratio(Q):
    rows = np.abs(Q).sum(axis=1)
    return float(np.divide(np.abs(np.diag(Q)), rows, out=np.zeros(len(Q)), where=rows != 0).sum())


def check_matrices():
    print("[1] Ky Fan-Hoffman and the one-block lemma", flush=True)
    rng = np.random.default_rng(431867)
    tested, worst_kfh, worst_o, worst_op = 0, 0.0, 0.0, 0.0
    for n in range(1, 13):
        for r in range(1, n + 1):
            for _ in range(25):
                Q = rng.normal(size=(n, r)) @ rng.normal(size=(r, n))
                if rng.random() < 0.3:  # nearly diagonal: the extremal regime
                    Q = Q * 0.05 + np.diag(rng.uniform(0, 1, n)) @ (np.eye(n) if r == n else np.zeros((n, n)))
                Q *= 10.0 ** rng.uniform(-6, 6)
                if rng.random() < 0.25:
                    Q[rng.integers(n), :] = 0
                rk = np.linalg.matrix_rank(Q)
                mass = np.abs(Q).sum()
                d = np.maximum(np.diag(Q), 0)
                k = kfh_ratio(Q)
                check(k <= rk + 1e-9, f"Ky Fan-Hoffman n={n} r={rk}: {k}")
                s = 10.0 ** rng.uniform(-8, 8) * (mass / n if mass else 1)
                lhs = np.minimum(s, d).sum()
                rhs = math.sqrt(s * rk * mass)
                check(lhs <= rhs * (1 + 1e-9) + 1e-300, f"one-block lemma n={n} r={rk}")
                lhs2 = np.sqrt(d).sum()
                check(lhs2 <= math.sqrt(rk * mass) * (1 + 1e-9) + 1e-300, f"parameter-free form n={n} r={rk}")
                if rk:
                    worst_kfh = max(worst_kfh, k / rk)
                if mass:
                    worst_o = max(worst_o, lhs / rhs)
                    worst_op = max(worst_op, lhs2 / math.sqrt(rk * mass))
                tested += 1
    # permuted diagonal: apply the lemma to Q Pi^T
    for _ in range(200):
        n = int(rng.integers(2, 9)); r = int(rng.integers(1, n + 1))
        Q = rng.normal(size=(n, r)) @ rng.normal(size=(r, n))
        pi = rng.permutation(n)
        s = float(rng.uniform(0.01, 2))
        lhs = sum(min(s, max(Q[i, pi[i]], 0)) for i in range(n))
        check(lhs <= math.sqrt(s * np.linalg.matrix_rank(Q) * np.abs(Q).sum()) + 1e-9, "one-block lemma, permuted")
        tested += 1
    for u in itertools.product((-1, 0, 1), repeat=3):
        for v in itertools.product((-1, 0, 1), repeat=3):
            Q = np.outer(u, v).astype(float)
            check(kfh_ratio(Q) <= 1 + 1e-12, "Ky Fan-Hoffman, rank-one sign patterns")
            check(np.sqrt(np.maximum(np.diag(Q), 0)).sum() ** 2 <= np.abs(Q).sum() + 1e-12, "rank-one sign patterns")
            tested += 1
    eq = 0
    for r, m in ((1, 3), (2, 3), (3, 2), (4, 2), (2, 4)):
        n = r * m * m
        Q = np.zeros((n, n))
        Q[: r * m, : r * m] = np.kron(np.eye(r), np.ones((m, m))) / (r * m * m)
        s = 1.0 / n
        lhs = np.minimum(s, np.maximum(np.diag(Q), 0)).sum()
        check(np.linalg.matrix_rank(Q) == r and abs(np.abs(Q).sum() - 1) < 1e-12, "equality example shape")
        check(abs(lhs - math.sqrt(r / n)) < 1e-12, f"equality in the one-block lemma r={r} m={m}")
        eq += 1
    print(f"  {tested} matrices; max KFH/rank {worst_kfh:.6f}, max lemma ratio {worst_o:.6f}, "
          f"max parameter-free ratio {worst_op:.6f}; {eq} equality cases exact", flush=True)


# ----------------------------------------------------------------------------- [2] Fact N
def check_fact_n():
    print("[2] normalized approximations have small mass and large overlap", flush=True)
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(3000):
        k = int(rng.integers(2, 40))
        P = rng.exponential(size=k) * (rng.random(k) < 0.6)
        if P.sum() == 0:
            P[0] = 1
        P /= P.sum()
        Q = P + rng.normal(scale=10.0 ** rng.uniform(-3, 1), size=k)
        Q += (1 - Q.sum()) / k
        eps = 0.5 * np.abs(P - Q).sum()
        om = np.minimum(P, np.maximum(Q, 0)).sum()
        t = np.maximum(-Q, 0).sum()
        check(abs(om - (1 - eps + t)) < 1e-12, "overlap identity: omega = 1 - eps + t")
        check(abs(np.abs(Q).sum() - (1 + 2 * t)) < 1e-12 and t <= eps + 1e-12, "||Q||_1 = 1 + 2t with t <= eps")
        check(om >= 1 - eps - 1e-12, "overlap >= 1 - eps")
        check(np.abs(Q).sum() <= 1 + 2 * eps + 1e-12, "||Q||_1 <= 1 + 2 eps")
        worst = max(worst, np.abs(Q).sum() - 1 - 2 * eps)
    print(f"  3000 random signed Q with sum 1: identity and both inequalities hold (max ||Q||_1 - 1 - 2eps = {worst:.2e})", flush=True)


# ----------------------------------------------------------------------------- [3] pendant blocks
def pendant_blocks(name, G, seed, svd_max=11):
    rng = random.Random(seed)
    n = G.number_of_nodes()
    E = edge_list(G)
    M = len(E)
    S = star(n, E)
    perm = list(range(M)); rng.shuffle(perm)
    inA = [False] * M
    for i in perm[: M // 2]:
        inA[i] = True
    mixed = mixed_vertices(n, E, inA)
    rng.shuffle(mixed)
    W = []
    for v in mixed:
        if all(not G.has_edge(v, w) for w in W):
            H = G.copy(); H.remove_nodes_from(W + [v])
            if H.number_of_nodes() >= 2 and nx.is_connected(H):
                W.append(v)
    if len(W) < 2:
        print(f"  {name}: fewer than 2 usable mixed vertices, skipped")
        return
    pairs = [(W[2 * i], W[2 * i + 1]) for i in range(len(W) // 2)]  # an odd vertex stays unpaired
    a = {w: rng.choice([i for i in S[w] if inA[i]]) for w in W}
    b = {w: rng.choice([i for i in S[w] if not inA[i]]) for w in W}
    smask = {w: sum(1 << i for i in S[w]) for w in W}
    keep = [x for x in range(n) if x not in set(W)]
    lam = algebraic_connectivity(laplacian(n, E, keep_vertices=set(keep)), keep)
    kappa = 1.0 / lam
    pb = {w: (1.0 / len(S[w])) * (1 + kappa) ** (-(len(S[w]) - 1)) for w in W}
    trees = [int(t) for t in enum_trees(n, E)]
    tau = len(trees)
    treeset = set(trees)

    def pattern(t, q):
        w, w2 = q
        sw, sw2 = t & smask[w], t & smask[w2]
        if sw == 1 << a[w] and sw2 == 1 << b[w2]:
            return 1
        if sw == 1 << b[w] and sw2 == 1 << a[w2]:
            return 2
        return 0

    blocks = defaultdict(list)
    gof = {}
    for t in trees:
        good = tuple(qi for qi, q in enumerate(pairs) if pattern(t, q))
        F = t
        for qi in good:
            F &= ~(smask[pairs[qi][0]] | smask[pairs[qi][1]])
        blocks[(good, F)].append(t)
        gof[t] = len(good)
    owner = {}
    rects = []
    ok_cube = ok_perm = ok_disj = True
    for key, ts in blocks.items():
        good, F = key
        g = len(good)
        ok_cube &= len(ts) == 2 ** g
        rows = [sum(1 << a[pairs[qi][c]] for qi, c in zip(good, bits)) for bits in itertools.product((0, 1), repeat=g)]
        cols = [sum(1 << b[pairs[qi][c]] for qi, c in zip(good, bits)) for bits in itertools.product((0, 1), repeat=g)]
        Mat = np.zeros((len(rows), len(cols)), dtype=int)
        for i, rs in enumerate(rows):
            for j, cs in enumerate(cols):
                cell = F | rs | cs
                Mat[i, j] = cell in treeset
                if cell in owner:
                    ok_disj = False
                owner[cell] = key
        ok_perm &= bool((Mat.sum(0) == 1).all() and (Mat.sum(1) == 1).all())
        rects.append((g, ts, [[F | rs | cs for cs in cols] for rs in rows]))
    check(ok_cube, f"{name}: pendant blocks are cubes")
    check(ok_perm, f"{name}: pendant rectangles are permutations")
    check(ok_disj, f"{name}: pendant rectangles are disjoint")

    worst_leaf = np.inf
    for w in W:
        others = [x for x in W if x != w]
        groups = defaultdict(lambda: defaultdict(int))
        for t in trees:
            groups[tuple(t & smask[x] for x in others)][t & smask[w]] += 1
        for cnt in groups.values():
            tot = sum(cnt.values())
            for e in S[w]:
                worst_leaf = min(worst_leaf, cnt.get(1 << e, 0) / tot / pb[w])
    check(worst_leaf >= 1 - 1e-9, f"{name}: conditional pendant bound (worst ratio {worst_leaf:.3f})")
    worst_pair = np.inf
    for q in pairs:
        others = [x for x in W if x not in q]
        groups = defaultdict(lambda: [0, 0])
        for t in trees:
            c = groups[tuple(t & smask[x] for x in others)]
            c[0] += 1; c[1] += pattern(t, q) > 0
        for tot, gd in groups.values():
            worst_pair = min(worst_pair, gd / tot / (2 * pb[q[0]] * pb[q[1]]))
    if pairs:
        check(worst_pair >= 1 - 1e-9, f"{name}: conditional pair bound (worst ratio {worst_pair:.3f})")

    gl = np.array([gof[t] for t in trees])
    gmax = int(gl.max())
    delta = {g: float((gl < g).mean()) for g in range(gmax + 2)}
    line = (f"  {name}: N={n} M={M} tau={tau} |W|={len(W)} pairs={len(pairs)} lambda2(G-W)={lam:.3f} "
            f"blocks={len(blocks)} max g={gmax} pendant ratio={worst_leaf:.2f} pair ratio={worst_pair:.2f}")
    Aidx = [i for i in range(M) if inA[i]]
    Bidx = [i for i in range(M) if not inA[i]]
    if len(Aidx) > svd_max or len(Bidx) > svd_max:
        print(line + "  (no SVD: split too large)", flush=True)
        return
    posA = {i: k for k, i in enumerate(Aidx)}
    posB = {i: k for k, i in enumerate(Bidx)}

    def rc(x):
        return (sum(1 << posA[i] for i in Aidx if x >> i & 1), sum(1 << posB[i] for i in Bidx if x >> i & 1))

    Pm = np.zeros((2 ** len(Aidx), 2 ** len(Bidx)))
    for t in trees:
        Pm[rc(t)] = 1.0 / tau
    U, sv, Vt = np.linalg.svd(Pm, full_matrices=False)
    rects_rc = [(g, [rc(t) for t in ts], [[rc(x) for x in row] for row in cells]) for g, ts, cells in rects]
    nrng = np.random.default_rng(seed)
    min_a = min_b = np.inf
    ntests = 0
    for r in range(1, min(len(sv), 2 ** (gmax + 1)) + 1):
        Qr = (U[:, :r] * sv[:r]) @ Vt[:r]
        cands = [(r, Qr)]
        u = nrng.normal(size=Pm.shape[0]); v = nrng.normal(size=Pm.shape[1])
        cands.append((r + 1, Qr + 0.3 * np.abs(Qr).max() * np.outer(u, v) / np.abs(np.outer(u, v)).max()))
        for rk, Q in cands:
            eps = 0.5 * np.abs(Pm - Q).sum()
            for g in range(gmax + 2):
                if delta[g] < 1:
                    sl = rk - 2 ** g * (1 - 2 * eps / (1 - delta[g]))
                    check(sl >= -1e-9, f"{name}: block theorem (b) r={rk} g={g}")
                    min_b = min(min_b, sl)
            lamQ = Q.sum()
            if lamQ <= 0:
                continue
            Qn = Q / lamQ
            eps = 0.5 * np.abs(Pm - Qn).sum()
            L1 = np.abs(Qn).sum()
            om = np.minimum(Pm, np.maximum(Qn, 0)).sum()
            check(om >= 1 - eps - 1e-12 and L1 <= 1 + 2 * eps + 1e-12, f"{name}: Fact N")
            for g, tcells, cells in rects_rc:
                ell = sum(abs(Qn[x]) for row in cells for x in row)
                omb = sum(min(1.0 / tau, max(Qn[x], 0.0)) for x in tcells)
                check(omb <= math.sqrt(rk * ell / tau) + 1e-12, f"{name}: one-block lemma on a block")
            for g in range(gmax + 2):
                check(om <= delta[g] + math.sqrt(rk * L1 / 2 ** g) + 1e-12, f"{name}: block overlap g={g}")
                sl = rk / 2 ** g - max(0.0, 1 - eps - delta[g]) ** 2
                check(sl >= -1e-12, f"{name}: block theorem (a) r={rk} g={g}")
                min_a = min(min_a, sl)
            ntests += 1
    print(line + f"  SVD tests={ntests} min slack (a)={min_a:.3g} (b)={min_b:.3g}", flush=True)


# ----------------------------------------------------------------------------- [4] positive part
def check_positive_part():
    print("[4] the positive part is needed: four triangles joined by three bridges", flush=True)
    # triangle c has edges 3c (e, in A), 3c+1 (f, in B), 3c+2 (arc, in A); bridges 12, 13, 14 in B
    trees = []
    for omitted in itertools.product(range(3), repeat=4):
        t = sum(1 << (12 + i) for i in range(3))
        for c, om in enumerate(omitted):
            t |= sum(1 << (3 * c + j) for j in range(3) if j != om)
        trees.append(t)
    tau = len(trees)
    blocks = {}
    for t in trees:
        good = [p for p in range(2) if (((t >> 6 * p) & 7), ((t >> (6 * p + 3)) & 7)) in ((5, 6), (6, 5))]
        special = sum((3 << (6 * p)) | (3 << (6 * p + 3)) for p in good)
        blocks.setdefault((tuple(good), t & ~special), []).append(t)
    delta2 = sum(len(ts) for (good, _), ts in blocks.items() if len(good) < 2) / tau
    eps = 1 - tau / 2 ** 15  # TV distance from UST to the uniform distribution, which has rank 1
    naive = 2 ** 2 * (1 - eps - delta2) ** 2
    naive_weak = naive / (1 + 2 * eps)
    stated = 2 ** 2 * max(0.0, 1 - eps - delta2) ** 2
    check(tau == 81 and abs(delta2 - 77 / 81) < 1e-12, "four triangles: tau = 81, delta(2) = 77/81")
    check(naive > 1 and naive_weak > 1, "without the positive part both forms would exceed the true rank 1")
    check(stated <= 1, "the stated bound holds")
    print(f"  tau={tau}, delta(2)={delta2:.4f}, eps={eps:.5f}: unsquared-margin bounds {naive:.3f} and "
          f"{naive_weak:.3f} (with 1+2eps) > rank 1; stated bound {stated:.3f} <= 1", flush=True)


# ----------------------------------------------------------------------------- [5] DV
def prism(k):
    E = [(i, (i + 1) % k) for i in range(k)] + [(k + i, k + (i + 1) % k) for i in range(k)] + [(i, k + i) for i in range(k)]
    return 2 * k, [tuple(sorted(e)) for e in E]


def squares(k, E, orient):
    idx = {e: i for i, e in enumerate(E)}
    out = []
    for j in range(0, k - 1, 2):
        top = idx[tuple(sorted((j, j + 1)))]; bot = idx[tuple(sorted((k + j, k + j + 1)))]
        r1 = idx[tuple(sorted((j, k + j)))]; r2 = idx[tuple(sorted((j + 1, k + j + 1)))]
        if orient == "tb":
            out.append(dict(cyc=[top, bot, r1, r2], e=top, f=bot))
        else:
            out.append(dict(cyc=[top, bot, r1, r2], e=r1, f=bot))
    return out


def dv_system(n, E, cycles, pairs, tag):
    trees = set(int(t) for t in enum_trees(n, E))
    tau = len(trees)
    cmask = [sum(1 << j for j in c["cyc"]) for c in cycles]

    def state(t, c):
        arcs = [j for j in c["cyc"] if j not in (c["e"], c["f"])]
        if not all(t >> j & 1 for j in arcs):
            return 0
        e, f = t >> c["e"] & 1, t >> c["f"] & 1
        return 1 if e and not f else (2 if f and not e else 0)

    blocks = set()
    for t in trees:
        good = tuple(q for q, (i, j) in enumerate(pairs) if (state(t, cycles[i]), state(t, cycles[j])) in ((1, 2), (2, 1)))
        free = sum(1 << cycles[k][x] for q in good for k in pairs[q] for x in ("e", "f"))
        blocks.add((good, t & ~free))
    seen = set()
    cells = zeros = 0
    for good, base in blocks:
        bits = list(itertools.product((0, 1), repeat=len(good)))
        for row in bits:
            for col in bits:
                z = base
                for q, x, y in zip(good, row, col):
                    ci, cj = (cycles[k] for k in pairs[q])
                    z |= 1 << (ci["e"] if x == 0 else cj["e"])
                    z |= 1 << (cj["f"] if y == 0 else ci["f"])
                target = z in trees
                psi = all((z & cm) != cm for cm in cmask)
                check(psi == target, f"DV {tag}: Psi = tau P on a block cell")
                check(target == (row == col), f"DV {tag}: permutation")
                check(z not in seen, f"DV {tag}: disjoint rectangles")
                seen.add(z)
                cells += 1; zeros += not target
    print(f"  {tag}: tau={tau} blocks={len(blocks)} cells={cells} (zero cells {zeros}): Psi = tau P on all", flush=True)


def check_dv():
    print("[5] one product-depth-2 circuit fills every cycle-pair block", flush=True)
    for k in (6, 8):
        n, E = prism(k)
        for orient, pairs in (("tb", [(0, 1)] if k == 6 else [(0, 2), (1, 3)]),
                              ("rb", [(0, 2)] if k == 6 else [(0, 3), (1, 2)])):
            dv_system(n, E, squares(k, E, orient), pairs, f"prism{k}-squares-{orient}")
    n, E = prism(6)
    for e, f in ((0, 1), (0, 3)):
        cs = [dict(cyc=list(range(sh, sh + 6)), e=e + sh, f=f + sh) for sh in (0, 6)]
        dv_system(n, E, cs, [(0, 1)], f"prism6-hexagons-{e}{f}")


# ----------------------------------------------------------------------------- [6] padding
def check_padding():
    print("[6] padding makes every slice highly cancelling (K4)", flush=True)
    G = nx.complete_graph(4)
    N = G.number_of_nodes()
    E = edge_list(G)
    M = len(E)
    ts = [int(t) for t in enum_trees(N, E)]
    tau = len(ts)
    c = 2 ** (M + N)
    ratios, zero = [], 0
    for pat in itertools.product((-1, 0, 1), repeat=M):
        fixed = [(i, v) for i, v in enumerate(pat) if v >= 0]
        cnt = sum(all((t >> i & 1) == v for i, v in fixed) for t in ts)
        if cnt == 0:
            zero += 1
            continue
        ratio = 1 + 2 * c * Fraction(1, 2 ** len(fixed)) / Fraction(cnt, tau)
        ratios.append(ratio)
        check(ratio >= 1 + 2 ** (N + 1), f"padded slice ratio {pat}")
    check(bool(ratios) and zero > 0, "positive- and zero-mass slices exercised")
    print(f"  {len(ratios)} positive-mass slices, {zero} zero-mass slices skipped; min ratio {float(min(ratios)):.1f} "
          f">= 1 + 2^(N+1) = {1 + 2 ** (N + 1)}", flush=True)


def main():
    check_matrices()
    check_fact_n()
    print("[3] pendant blocks with one-vertex gadgets, exact enumeration", flush=True)
    graphs = ["K33", "cube", "petersen", "heawood", "mobius_kantor", "rr3_14", "rr3_16"]
    for gi, name in enumerate(graphs):
        G = named_graph(name, seed=3 + gi)
        for seed in range(3):
            pendant_blocks(f"{name}/s{seed}", G, seed)
    check_positive_part()
    check_dv()
    check_padding()
    if FAIL:
        print(f"\n{len(FAIL)} CHECK(S) FAILED")
        sys.exit(1)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()

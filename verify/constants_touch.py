# /// script
# dependencies = []
# ///
"""The constants of Theorem "2-connected expanders" and Proposition "parametric bounds", and Table "constants".

All quantities are explicit formulas; no graph is involved. Columns of the table (entries are ln(1/c),
lower is better; c is the constant in beta^2 <= e^{-cN} for balanced splits, so kappa = 3c):
  old     : the constant of the previous version, 2h/(3 D^{s+1}(D+h)) (2+16D^3/h^2)^{-D}, s = 2 ceil(2D/h)+1;
  closed  : Theorem part (a), c = h1^{5D+4} (D+1)^{-4(D+1)} e^{-10(D+1)}, h1 = min(h,1);
  P'      : Proposition (i) at h' = h1/(D+1), q0 = (8 D Lam^3)^{-(D+1)}, Lam = e(1+1/h1), with
            lambda_2(H) >= h'^2/(2D) (the exact formula behind the closed form);
  P opt   : Proposition (i) with mu(D,h') and the largest q with S(q) <= 1/2 found by bisection in
            floating point, maximized over h' = (h/(D+1)) i/12, i = 1..12. NOT interval-certified;
  spectral: Theorem part (b), c = lam (D-1)(1-1/lam)^D / (3 D^3 (D+lam/2)(D^2+1)).
Checks (exit status 1 on failure):
  * the closed form: ln(3/kappa_P') <= (5D+4) ln(1/h) + 4(D+1) ln(D+1) + 10(D+1) for D = 2..300 and
    400 values of h in (0,1] (the paper proves it; this is a numerical cross-check);
  * S(q0) <= 1/3 for several (D, h), including D = 2 (proved in the paper; cross-check);
  * the 6-regular Ramanujan value c = 8.011e-8 at lam = 6 - 2 sqrt 5.
"""
import math, sys

LN = math.log
FAIL = []


def check(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        FAIL.append(msg)


def lmix(D, h):            # ln of |Z|/N lower bound / 1: 2h/(3D(D+h)), using m0 = M/3 >= N/3
    return LN(2 * h / (3 * D * (D + h)))


def lquality(D, kappa):    # Lemma "telescoped touch bound": pA pB >= ((D-1)/D^2)(1+kappa)^{-D}
    return LN((D - 1) / D**2) - D * LN(1 + kappa)


def mu(D, i):
    return i * i / (D + math.sqrt(D * D - i * i))   # = D - sqrt(D^2 - i^2), stable form


def ln_old(D, h):
    s = 2 * math.ceil(2 * D / h) + 1
    return -(lmix(D, h) - s * LN(D) - D * LN(2 + 16 * D**3 / h**2))


def ln_closed(D, h):
    h1 = min(h, 1.0)
    return (5 * D + 4) * LN(1 / h1) + 4 * (D + 1) * LN(D + 1) + 10 * (D + 1)


def ln_Pprime(D, h):
    h1 = min(h, 1.0); hp = h1 / (D + 1)
    lnLam = 1 + LN(1 + 1 / h1)
    lnq = -(D + 1) * (LN(8 * D) + 3 * lnLam)
    return -(lmix(D, h1) + lnq - LN(2 * (D + 1)) + lquality(D, 2 * D / hp**2))


def ln_spec(D, lam):
    return -(LN(lam * (D - 1) / (3 * D**3 * (D + lam / 2) * (D * D + 1))) + D * LN(1 - 1 / lam))


# ---------------------------------------------------------------- the series S(q) of Lemma "random sparse deletion"
def lbinom(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def lse(a, b):
    if a == -math.inf:
        return b
    m = max(a, b)
    return m + LN(math.exp(a - m) + math.exp(b - m))


def ln_S(D, h, hp, lnq, B=3000):
    """ln S(q), S(q) = D sum_{b>=2} F(b) binom(b-1, m(b)-1) q^{m(b)-1}, F(b) <= binom(K+b+1, b+1), K = floor(b/h).
    Terms b <= B are summed; the tail is bounded by a geometric series using F(b) <= Lam^{b+1},
    binom(b-1, m-1) <= 2^b and m(b)-1 >= theta b - 1."""
    theta = (1 - hp / h) / D
    tot = -math.inf
    for b in range(2, B + 1):
        K = math.floor(b / h)
        if K < 1:
            continue
        m = max(2, math.floor(b * theta) + 1)
        tot = lse(tot, lbinom(K + b + 1, b + 1) + lbinom(b - 1, m - 1) + (m - 1) * lnq)
    lnLam = 1 + LN(1 + 1 / h)
    lr = lnLam + LN(2) + theta * lnq
    if lr >= 0:
        return math.inf
    tot = lse(tot, lnLam - lnq + (B + 1) * lr - LN(1 - math.exp(lr)))
    return LN(D) + tot


def best_lnq(D, h, hp):
    lo, hi = -2000.0, 0.0
    while ln_S(D, h, hp, lo) > LN(0.5):
        lo *= 2
        assert lo > -1e7, "no feasible q"
    for _ in range(50):
        mid = (lo + hi) / 2
        if ln_S(D, h, hp, mid) <= LN(0.5):
            lo = mid
        else:
            hi = mid
    return lo


def ln_Popt(D, h, grid=12):
    best = math.inf
    for i in range(1, grid + 1):
        hp = h / (D + 1) * i / grid
        val = -(lmix(D, h) + best_lnq(D, h, hp) - LN(2 * (D + 1)) + lquality(D, 1 / mu(D, hp)))
        best = min(best, val)
    return best


# ---------------------------------------------------------------- checks
worst = -math.inf
for D in range(2, 301):
    for t in range(400):
        h = math.exp(-12 * t / 399)          # h from 1 down to e^-12
        worst = max(worst, ln_Pprime(D, h) - ln_closed(D, h))
check(worst <= 0, f"closed form: ln(3/kappa_P') - [(5D+4)ln(1/h) + 4(D+1)ln(D+1) + 10(D+1)] <= 0 on the grid "
      f"(largest value {worst:.2f})")
ok, big = True, -math.inf
for D in (2, 3, 4, 6, 10, 32):
    for h in (0.01, 0.05, 0.2, 0.5, 1.0):
        lnLam = 1 + LN(1 + 1 / h)
        lnq0 = -(D + 1) * (LN(8 * D) + 3 * lnLam)
        v = ln_S(D, h, h / (D + 1), lnq0)
        ok &= v <= LN(1 / 3); big = max(big, v)
check(ok, f"S(q0) <= 1/3 at h' = h/(D+1) for D in (2,3,4,6,10,32), h in (0.01,...,1) (largest S = {math.exp(big):.3g})")
lam = 6 - 2 * math.sqrt(5)
cR = math.exp(-ln_spec(6, lam))
check(abs(cR - 8.0113e-8) < 1e-11, f"6-regular Ramanujan, lambda_2 >= 6 - 2 sqrt 5: c = {cR:.6g}")
cs = {D: math.exp(-ln_spec(D, D - 2 * math.sqrt(D - 1))) for D in range(6, 201)}
low = min(cs[D] for D in range(6, 65))
check(low >= 2.7e-9, f"Ramanujan bound lambda = D - 2 sqrt(D-1): min over 6 <= D <= 64 of c_lambda(D) = {low:.4g} >= 2.7e-9")
print("      c_lambda at the Ramanujan bound: " + ", ".join(f"D={D}: {cs[D]:.3g}" for D in (6, 7, 8, 16, 64, 100, 200))
      + f"; D^4 c_lambda(D) for D in [64, 200] lies in [{min(D**4 * cs[D] for D in range(64, 201)):.3f}, {max(D**4 * cs[D] for D in range(64, 201)):.3f}]")
c35 = math.exp(-ln_closed(3, 0.2))
print(f"      D=3, h=1/5: closed form c = {c35:.3g}, P' c = {math.exp(-ln_Pprime(3, 0.2)):.3g}, old c = {math.exp(-ln_old(3, 0.2)):.3g}")

# ---------------------------------------------------------------- table
print("\nTable: ln(1/c); 'P opt' is a floating-point optimization, not interval-certified.")
print(f"{'D':>3} {'h':>6} {'lambda':>7} | {'old':>7} {'closed':>7} {'Pprime':>7} {'P opt':>7} {'spectral':>8} | best c")
rows = [(3, 0.05, None), (3, 0.1, None), (3, 0.2, None), (4, 0.1, None), (4, 0.2, None), (4, 0.4, None),
        (6, 0.2, None), (6, 0.5, None), (6, 0.76, None)]
for D in (6, 8, 10, 16, 64):
    lamD = D - 2 * math.sqrt(D - 1)
    rows.append((D, lamD / 2, lamD))
for D, h, lamD in rows:
    vals = [ln_old(D, h), ln_closed(D, h), ln_Pprime(D, h), ln_Popt(D, h)]
    sp = ln_spec(D, lamD) if lamD else math.inf
    best = min(vals[1:] + [sp])
    print(f"{D:>3} {h:>6.3f} {'' if lamD is None else f'{lamD:.3f}':>7} | " + " ".join(f"{v:7.1f}" for v in vals)
          + f" {'' if lamD is None else f'{sp:8.1f}':>8} | {math.exp(-best):.2e}", flush=True)

print("ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES")
sys.exit(1 if FAIL else 0)

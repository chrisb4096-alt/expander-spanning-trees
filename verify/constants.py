# /// script
# dependencies = []
# ///
"""Exact-arithmetic checks of the constants: the cycle-bit constants of the main-results section
(Corollary "closed-form constants" and Example "cubic prisms"), and the vertex-touch constant c(D,h)
of the previous version of the paper, for comparison (the current constant of Theorem 1.1 is checked in
constants_touch.py).

Everything is a rational computation except the two printed decimal values of the per-cycle factors,
which use the standard library's float functions only for display. The final exponent is certified
with the exact inequality -log(1-x) >= x. Exits nonzero (AssertionError) on any failure.
"""
from fractions import Fraction as F
import math


def t_lower(h, D, ell, N, M):
    """Lemma "many mixed cycles" with m0 = M/3: t >= min{ hN/(4 D ell), 2 m0 / (ell D (1 + D/h)) }."""
    return min(h * N / (4 * D * ell), 2 * (F(M) / 3) / (ell * D * (1 + D / h)))


# ---------- Example "cubic prisms": prism over a cubic H, as a function of h' = min(1, h_H) ----------
for hp in [F(1, 100), F(1, 20), F(1, 10), F(1, 5), F(1, 2), F(1)]:
    N = F(1)                      # work per vertex: N = 1, M = 2N for a 3-regular base
    t = t_lower(hp / 3, 4, 4, N, 2 * N)
    assert t == hp / 192, (hp, t)  # the first term is the minimum whenever h' <= 4
    r = t / (4 * 3)                # ell * kappa = 12
    assert r == hp / 2304
    w = F(1, 256)                  # star-framed weight D^-ell = 4^-4 (Theorem "linked expanders")
    c_lower = r * w                # -log(1-w) >= w
    assert c_lower == hp / 589824
print("cubic prism: t/N >= h'/192, r/N >= h'/2304, c >= h'/589824 for h' in {1/100,...,1}")

hp = F(1, 5)
t = t_lower(hp / 3, 4, 4, F(1), F(2))
assert t == F(1, 960)
second = 2 * (F(2) / 3) / (4 * 4 * (1 + 4 / (hp / 3)))
assert second == F(1, 732)
r = t / 12
assert r == F(1, 11520)
c_lower = r * F(1, 256)
assert c_lower == F(1, 2949120) and c_lower > F(3, 10**7)
print(f"h_H = 1/5: t/N >= {t} (other term {second}), r/N >= {r}, "
      f"c >= 1/(256*11520) = {float(c_lower):.6e} > 3e-7")
print(f"  per-cycle factor sqrt(1-1/256) = {math.sqrt(1 - 1/256):.10f}; "
      f"sharper sqrt((1+sqrt(1-4/256))/2) = {math.sqrt((1 + math.sqrt(1 - 4/256)) / 2):.10f}; "
      f"-log(1-1/256)/11520 = {-math.log1p(-1/256)/11520:.6e}")
closed = F(1, 5) / (24 * 3 * 4 ** 5 * 13)
assert closed == F(1, 4792320) and closed <= c_lower
print(f"  closed form of Corollary 'closed-form constants' at Delta = 3, h' = 1/5: 1/4792320 = {float(closed):.4e}")
# first version of the paper used D^-2(ell-1) = 4^-6: 16 times weaker
assert F(1, 256) / F(1, 4096) == 16

# ---------- Corollary "closed-form constants": general closed forms ----------
# Per-cycle weight D^-ell (Theorem "affinity bound for a split" with the star-framed bound).
# prism: D <= Delta+1, ell = 4, kappa <= Delta, h >= h'/3, M >= N; c >= h'/(24 Delta (Delta+1)^5 (3 Delta + 4)).
# square: D <= Delta^2, ell = 3, kappa <= 3 Delta (Delta-1)/2, h >= h_H, M >= 2N/3;
#         c >= 8 h_H / (243 Delta^9 (Delta-1) (h_H + Delta^2)).
checked = 0
for Delta in range(2, 13):
    for hH in [F(1, 1000), F(1, 50), F(1, 10), F(1, 3), F(1), F(2), F(Delta)]:
        hp = min(F(1), hH)
        N = F(1)
        # prism
        D, ell, kappa = Delta + 1, 4, Delta
        t = t_lower(hp / 3, D, ell, N, N)
        r = t / (ell * kappa)
        exact = r * F(1, D ** ell)
        closed = hp / (24 * Delta * (Delta + 1) ** 5 * (3 * Delta + 4))
        assert closed <= exact, ("prism", Delta, hH, closed, exact)
        # square
        D, ell, kappa = Delta ** 2, 3, F(3 * Delta * (Delta - 1), 2)
        t = t_lower(hH, D, ell, N, F(2, 3) * N)
        r = t / (ell * kappa)
        exact = r * F(1, D ** ell)
        closed = 8 * hH / (243 * Delta ** 9 * (Delta - 1) * (hH + Delta ** 2))
        assert closed <= exact, ("square", Delta, hH, closed, exact)
        checked += 1
print(f"general closed forms for prisms and squares hold at {checked} (Delta, h_H) pairs, Delta = 2..12")


# ---------- The vertex-touch constant c(D,h) of the previous version (column "previous" of the table) ----------
# s = 2*ceil(2D/h)+1, q = (2 + 16 D^3/h^2)^(-D), kappa = 2 h q / (D^(s+1) (D+h)), c = kappa/3.
def c_touch(D, h):
    t = -((-2 * D) // h)          # ceil(2D/h) for a Fraction h
    assert t - 1 < F(2 * D) / h <= t
    s = 2 * t + 1
    q = 1 / (2 + 16 * F(D) ** 3 / h ** 2) ** D
    kappa = 2 * h * q / (D ** (s + 1) * (D + h))
    return s, q, kappa, kappa / 3

def log10(x):
    return math.log10(x.numerator) - math.log10(x.denominator)

for D, h in [(3, F(1, 5)), (3, F(18, 100)), (3, F(1)), (4, F(1)), (10, F(1))]:
    s_, q, kappa, c = c_touch(D, h)
    assert c > 0 and kappa == 3 * c
    print(f"previous c(D={D}, h={h}) = 10^{log10(c):.3f}   (s = {s_}, q = 10^{log10(q):.3f})")
s_, q, kappa, c = c_touch(3, F(1, 5))
assert s_ == 61 and q == F(1, 10802 ** 3)
assert F(85, 10 ** 45) < c < F(95, 10 ** 45), float(c)   # "about 9e-44" in the paper
print(f"previous constant, D = 3, h = 1/5: s = 61, q = 1/10802^3, c = {float(c):.4e} (paper: about 9e-44)")

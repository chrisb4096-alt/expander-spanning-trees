# /// script
# dependencies = []
# ///
"""Exact-arithmetic checks of the constants in Section 6 (Corollary 6.2 and Example 6.3).

Everything is a rational computation except the two printed decimal values of the per-cycle factors,
which use the standard library's float functions only for display. The final exponent is certified
with the exact inequality -log(1-x) >= x. Exits nonzero (AssertionError) on any failure.
"""
from fractions import Fraction as F
import math


def t_lower(h, D, ell, N, M):
    """Lemma 5.3 with m0 = M/3: t >= min{ hN/(4 D ell), 2 m0 / (ell D (1 + D/h)) }."""
    return min(h * N / (4 * D * ell), 2 * (F(M) / 3) / (ell * D * (1 + D / h)))


# ---------- Example 6.3: prism over a cubic H, as a function of h' = min(1, h_H) ----------
for hp in [F(1, 100), F(1, 20), F(1, 10), F(1, 5), F(1, 2), F(1)]:
    N = F(1)                      # work per vertex: N = 1, M = 2N for a 3-regular base
    t = t_lower(hp / 3, 4, 4, N, 2 * N)
    assert t == hp / 192, (hp, t)  # the first term is the minimum whenever h' <= 4
    r = t / (4 * 3)                # ell * kappa = 12
    assert r == hp / 2304
    p = F(1, 64)                   # every cycle vertex has degree 4: 4^3 = 64
    c_lower = r * p * p            # -log(1-p^2) >= p^2
    assert c_lower == hp / 9437184
print("cubic prism: t/N >= h'/192, r/N >= h'/2304, c >= h'/9437184 for h' in {1/100,...,1}")

hp = F(1, 5)
t = t_lower(hp / 3, 4, 4, F(1), F(2))
assert t == F(1, 960)
second = 2 * (F(2) / 3) / (4 * 4 * (1 + 4 / (hp / 3)))
assert second == F(1, 732)
r = t / 12
assert r == F(1, 11520)
c_lower = r * F(1, 4096)
assert c_lower == F(1, 47185920) and c_lower > F(2, 10**8)
print(f"h_H = 1/5: t/N >= {t} (other term {second}), r/N >= {r}, "
      f"c >= 1/(4096*11520) = {float(c_lower):.6e} > 2e-8")
print(f"  per-cycle factor sqrt(1-1/4096) = {math.sqrt(1 - 1/4096):.10f}; "
      f"sharper sqrt((1+sqrt(1-4/4096))/2) = {math.sqrt((1 + math.sqrt(1 - 4/4096)) / 2):.10f}; "
      f"-log(1-1/4096)/11520 = {-math.log1p(-1/4096)/11520:.6e}")

# ---------- Corollary 6.2: general closed forms ----------
# prism: D <= Delta+1, ell = 4, kappa <= Delta, h >= h'/3, M >= N; c >= h'/(24 Delta (Delta+1)^7 (3 Delta + 4)).
# square: D <= Delta^2, ell = 3, kappa <= 3 Delta (Delta-1)/2, h >= h_H, M >= 2N/3;
#         c >= 8 h_H / (243 Delta^11 (Delta-1) (h_H + Delta^2)).
checked = 0
for Delta in range(2, 13):
    for hH in [F(1, 1000), F(1, 50), F(1, 10), F(1, 3), F(1), F(2), F(Delta)]:
        hp = min(F(1), hH)
        N = F(1)
        # prism
        D, ell, kappa = Delta + 1, 4, Delta
        t = t_lower(hp / 3, D, ell, N, N)
        r = t / (ell * kappa)
        exact = r * F(1, D ** (2 * (ell - 1)))
        closed = hp / (24 * Delta * (Delta + 1) ** 7 * (3 * Delta + 4))
        assert closed <= exact, ("prism", Delta, hH, closed, exact)
        # square
        D, ell, kappa = Delta ** 2, 3, F(3 * Delta * (Delta - 1), 2)
        t = t_lower(hH, D, ell, N, F(2, 3) * N)
        r = t / (ell * kappa)
        exact = r * F(1, D ** (2 * (ell - 1)))
        closed = 8 * hH / (243 * Delta ** 11 * (Delta - 1) * (hH + Delta ** 2))
        assert closed <= exact, ("square", Delta, hH, closed, exact)
        checked += 1
print(f"general closed forms for prisms and squares hold at {checked} (Delta, h_H) pairs, Delta = 2..12")

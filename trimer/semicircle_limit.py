"""Large-K limit of the pair thresholds at fixed t* = sqrt(K) t (semicircle DOS of half-width 2 t*).
Uniform: U2^u/t* -> 2 (from 2(K-1)/sqrt K). Square-summable: only d = 0 survives in
1/U2 = sum_d b_d P_d(2 sqrt K) since b_d P_d ~ K^{-d/2}; b0 = int int rho rho /(e+e'+4).
Also: fits a + b/sqrt K of the Table I ratios, calibrated on the pair ratio whose limit is exact."""
import numpy as np
from scipy import integrate
rho = lambda e: np.sqrt(max(4 - e * e, 0)) / (2 * np.pi)
b0, _ = integrate.dblquad(lambda e2, e1: rho(e1) * rho(e2) / (e1 + e2 + 4), -2, 2, -2, 2, epsabs=1e-10)
print(f"b0 = {b0:.6f}  U2/t* -> {1/b0:.4f}  U2u/t* -> 2  ratio U2/U2u -> {1/(2*b0):.4f}")
K = np.array([2, 3, 4, 5, 6, 10.]); x = 1 / np.sqrt(K)
U2 = np.array([1.966, 3.251, 4.269, 5.136, 5.904, 8.402]); U2u = 2 * (K - 1) / np.sqrt(K)
U3 = np.array([1.369, 2.159, 2.741, 3.218, 3.633, 4.949]); U3u = np.array([1.310, 2.073, 2.637, 3.103, 3.508, 4.800])
print("cross-sector U3/U2u:", np.round(U3 / U2u, 3))
for name, r in (("pair L2/u (exact limit %.3f)" % (1 / (2 * b0)), U2 / U2u), ("trimer L2", U3 / U2), ("trimer uniform", U3u / U2u), ("trimer cross", U3 / U2u)):
    a, b = np.polyfit(x, r, 1)[::-1]; print(f"{name}: a + b/sqrt K -> a = {a:.3f}")

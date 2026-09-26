"""Square-summable pair threshold on the clean tree via the radial kernel.
B_ij = G2(E; (i,i),(j,j)) = b_{d(i,j)},  B = sum_d b_d A_d = f(A),
f(mu) = sum_d b_d P_d(mu), ||B|| = f(2 sqrt K) (positive kernel, spherical
functions bounded by their edge value). U_2^{L2} = 1/f(2 sqrt K).
Check: uniform-sector value sum_d n_d b_d = 1/U_2^{unif} = sqrt(K)/(2(K-1)).
b_d = int int rho(l) rho(l') phi_d(l) phi_d(l') / (4 sqrt K - l - l') dl dl',
phi_d = P_d / n_d, rho = Kesten-McKay (t = 1)."""
import numpy as np, sys
def thresholds(K, n=1500, D=120):
    r = 2 * np.sqrt(K)
    # nodes clustered at the corner theta -> 0 (lambda -> +2 sqrt K): lambda = r cos(theta), theta = pi s^2
    s, w = np.polynomial.legendre.leggauss(n); s = 0.5 * (s + 1); w = 0.5 * w
    th = np.pi * s ** 2; dth = 2 * np.pi * s * w
    lam = r * np.cos(th)
    rho = (K + 1) / (2 * np.pi) * np.sqrt(np.maximum(4 * K - lam ** 2, 0)) / ((K + 1) ** 2 - lam ** 2)
    W = rho * dth * r * np.sin(th)          # rho(lambda) dlambda
    P = np.zeros((D + 1, n)); P[0] = 1; P[1] = lam; P[2] = lam ** 2 - (K + 1)
    for d in range(2, D): P[d + 1] = lam * P[d] - K * P[d - 1]
    nd = np.array([1] + [(K + 1) * K ** (d - 1) for d in range(1, D + 1)], float)
    phi = P / nd[:, None]
    h = np.sin(th / 2) ** 2
    M = np.outer(W, W) / (2 * r * (h[:, None] + h[None, :]))
    b = np.einsum("di,ij,dj->d", phi, M, phi)
    Pe = np.array([1] + [K ** (d / 2) * ((K + 1) + d * (K - 1)) / K for d in range(1, D + 1)])  # P_d(2 sqrt K)
    fL2 = np.sum(b * Pe); funi = np.cumsum(b * nd)
    return 1 / fL2, funi, b
if __name__ == "__main__":
    for K in (2, 3, 4, 5, 6, 10):
        for n, D in ((800, 80), (1500, 120), (2500, 160)):
            U2, funi, b = thresholds(K, n, D)
            # uniform sum converges like 1/D: Richardson on the last two partial sums
            uni = 1 / funi[-1]; uni_x = 1 / (2 * funi[-1] - funi[D // 2])
            print(f"K={K} n={n} D={D}: U2_L2={U2:.6f}   uniform check: 1/sum={uni:.4f} (Richardson {uni_x:.4f}) exact {2*(K-1)/np.sqrt(K):.4f}", flush=True)

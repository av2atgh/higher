"""T=0 BCS mean field of two colours with on-site attraction U on the Bethe lattice
(Kesten-McKay DOS, t=1), density n per colour per site.
 gap:    1 = (U/2) int rho(e) / E(e),  E = sqrt((e-mu)^2 + D^2)
 number: n = int rho(e) v^2,  v^2 = (1 - (e-mu)/E)/2
 energy per site: E_SF = 2 int rho e v^2 - U n^2 - D^2/U   (Hartree -U n^2 included)
Dilute limit: 2 mu -> uniform-sector pair energy E_2^u (the condensate is the Perron mode)."""
import numpy as np
from scipy.optimize import brentq
from tree3 import pair_energy
def dos(K, n=4000):
    r = 2 * np.sqrt(K)
    th = (np.arange(n) + 0.5) * np.pi / n; e = r * np.cos(th)
    rho = (K + 1) / (2 * np.pi) * np.sqrt(np.maximum(4 * K - e ** 2, 0)) / ((K + 1) ** 2 - e ** 2)
    w = rho * r * np.sin(th) * np.pi / n
    return e, w
def solve(K, U, n, e=None, w=None):
    if e is None: e, w = dos(K)
    def gap_res(mu, D): return (U / 2) * np.sum(w / np.sqrt((e - mu) ** 2 + D ** 2)) - 1
    def num(mu, D): return np.sum(w * 0.5 * (1 - (e - mu) / np.sqrt((e - mu) ** 2 + D ** 2)))
    # for given mu find D from the gap equation, then adjust mu for the density
    def D_of_mu(mu):
        f = lambda D: gap_res(mu, D)
        if f(1e-12) < 0: return 0.0
        return brentq(f, 1e-12, 50.0)
    g = lambda mu: num(mu, D_of_mu(mu)) - n
    mu = brentq(g, -30.0, 2 * np.sqrt(K) + U + 1)
    D = D_of_mu(mu)
    v2 = 0.5 * (1 - (e - mu) / np.sqrt((e - mu) ** 2 + D ** 2))
    E = 2 * np.sum(w * e * v2) - U * n ** 2 - D ** 2 / U
    return mu, D, E
def free(K, n, e=None, w=None):
    """energy per site of one colour of free fermions at density n."""
    if e is None: e, w = dos(K)
    c = np.cumsum(w[::-1])[::-1]  # descending e -> ascending fill from the bottom: reorder
    order = np.argsort(e); es, ws = e[order], w[order]
    cum = np.cumsum(ws); i = np.searchsorted(cum, n)
    return np.sum(ws[:i] * es[:i]) + (n - (cum[i - 1] if i > 0 else 0)) * es[i], es[i]
if __name__ == "__main__":
    for K in (2, 3):
        e, w = dos(K)
        print(f"K={K}: dilute-limit check, 2 mu vs E_2^u:")
        for U in (2.0, 3.0, 4.0):
            for n in (1e-2, 1e-3, 1e-4):
                mu, D, E = solve(K, U, n, e, w)
                print(f"   U={U} n={n:.0e}: 2mu={2*mu:.5f}  E2u={pair_energy(K, U):.5f}  E/n={E/n:.5f}  D={D:.4f}")
        for U in (1.0, 1.5):
            mu, D, E = solve(K, U, 0.05, e, w); print(f"   U={U} n=0.05: mu={mu:.4f} D={D:.4f} E={E:.5f}  free(2n colours)={2*free(K,0.05,e,w)[0]:.5f}")

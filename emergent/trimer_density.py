"""One-body radial densities of the composites, for the effective disorder of a composite
(uniform sector, tree3) and for the support radius of a site-bound trimer (root-fixed sector, tree4).
 s_m^2 = sum_d rho(d)^2 / n_d  (variance factor of the site-energy shift of the composite; compact trimer: 3, compact pair: 2)
 r90   = smallest radius holding 90% of the one-body density of the site-bound trimer."""
import numpy as np, scipy.sparse.linalg as spl, scipy.sparse as sp, json, sys, time
sys.path.insert(0, "../trimer")
import tree3, tree4
from tree3 import pair_energy

def n_d(K, d): return 1 if d == 0 else (K + 1) * K ** (d - 1)

def trimer_uniform_density(K, U, L=60):
    H, S = tree3.build(K, U, L)
    w, v = spl.eigsh(H, k=1, which="SA", tol=1e-10, maxiter=20000)
    g2 = v[:, 0] ** 2; g2 /= g2.sum()
    rho = {}
    for s, q in zip(S, g2):
        for d in s: rho[d] = rho.get(d, 0.0) + q
    dmax = max(rho); r = np.array([rho.get(d, 0.0) for d in range(dmax + 1)])
    s2 = sum(r[d] ** 2 / n_d(K, d) for d in range(dmax + 1))
    return float(w[0]), r, float(np.sqrt(s2))

def pair_uniform_density(K, U, L=200):
    """half line in d: site energy -U at d=0, hoppings J0 then J; density of the second particle relative to the first."""
    J = 2 * np.sqrt(K); J0 = 2 * np.sqrt(K + 1)
    diag = np.zeros(L + 1); diag[0] = -U
    off = np.full(L, -J); off[0] = -J0
    H = sp.diags([off, diag, off], [-1, 0, 1]).tocsr()
    w, v = spl.eigsh(H, k=1, which="SA", tol=1e-12)
    g2 = v[:, 0] ** 2; g2 /= g2.sum()
    # density: reference particle at the origin (1) plus the partner at distance d with prob g2[d]
    s2 = (1 + g2[0]) ** 2 + sum(g2[d] ** 2 / n_d(K, d) for d in range(1, L + 1))
    return float(w[0]), g2, float(np.sqrt(s2))

class SiteTrimer:
    """root-fixed sector at fixed (K, L): build the kinetic part once, vary (U, eps0)."""
    def __init__(self, K, L):
        t0 = time.time(); self.K, self.L = K, L
        H0, self.S = tree4.build(3, K, 0.0, 0.0, L)
        self.Hk = H0
        self.nroot = np.array([sum(1 for x in s if len(x) == 0) for s in self.S], float)
        self.npair = np.array([sum(1 for a in range(3) for b in range(a + 1, 3) if s[a] == s[b]) for s in self.S], float)
        self.depths = [[len(x) for x in s] for s in self.S]
        print(f"K={K} L={L}: {len(self.S)} states ({time.time()-t0:.0f}s)", flush=True)
    def density(self, U, eps0):
        H = self.Hk + sp.diags(eps0 * self.nroot - U * self.npair)
        w, v = spl.eigsh(H, k=1, which="SA", tol=1e-9, maxiter=50000)
        g2 = v[:, 0] ** 2; g2 /= g2.sum()
        rho = np.zeros(self.L + 1)
        for q, ds in zip(g2, self.depths):
            for d in ds: rho[d] += q
        cum = np.cumsum(rho) / 3
        r90 = int(np.argmax(cum >= 0.9)); r99 = int(np.argmax(cum >= 0.99))
        return float(w[0]), rho, r90, r99

if __name__ == "__main__":
    what = sys.argv[1]
    if what == "uniform":
        out = {}
        for K in (2, 3):
            Us = [1.38, 1.4, 1.45, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 4.0] if K == 2 else [2.17, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 3.0, 3.5, 4.0, 5.0]
            for U in Us:
                e3, r3, s3 = trimer_uniform_density(K, U); e2, g2, s2 = pair_uniform_density(K, U)
                out[f"{K},{U}"] = dict(E3u=e3, rho3=r3.tolist(), s3=s3, E2u=e2, s2=s2, rho3_0=float(r3[0]), g2_0=float(g2[0]))
                print(f"K={K} U={U:.2f}: trimer rho(0)={r3[0]:.3f} rho(1)={r3[1]:.3f} s3={s3:.3f} | pair g0^2={g2[0]:.3f} s2={s2:.3f}", flush=True)
        json.dump(out, open("uniform_density.json", "w"))
    elif what == "site":
        K = int(sys.argv[2]); L = int(sys.argv[3]); st = SiteTrimer(K, L)
        ec = tree4.impurity_threshold(K)
        r = np.array(json.load(open(f"../trimer/scan_bs_K{K}_L60.json")))  # eps0, frac, L, u2, u3
        fr, u2, u3 = r[:, 1], r[:, 3], r[:, 4]
        out = []
        for f in (0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95):
            eps0 = -f * ec; U3 = np.interp(f, fr, u3); U2 = np.interp(f, fr, u2)
            for x in (0.1, 0.5, 0.9):   # position inside the Borromean window (U3, U2)
                U = U3 + x * (U2 - U3)
                e, rho, r90, r99 = st.density(U, eps0)
                out.append(dict(K=K, f=f, eps0=eps0, U=U, x=x, E=e, edge=-6 * np.sqrt(K), rho=rho[:12].tolist(), r90=r90, r99=r99))
                print(f"K={K} eps0={-f:.2f}ec U={U:.3f} (x={x}): E-edge={e+6*np.sqrt(K):+.4f}  rho(0..4)={np.round(rho[:5],3)}  r90={r90} r99={r99}", flush=True)
                json.dump(out, open(f"site_density_K{K}.json", "w"))

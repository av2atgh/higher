"""Dilute-limit competition on the Bethe lattice at density n per colour:
 trion liquid: free fermions of mass set by t_3(U), band e_3 + [-2 sqrt K t_3, 2 sqrt K t_3],
   t_3 = (E3_L2 - E3_u)/(K+1-2 sqrt K), e_3 = E3_u + (K+1) t_3 (E3_L2 extrapolated 1/L^2)
 colour superfluid: BCS of two colours (bcs.solve) + free third colour (bcs.free)
 mixed: E(n) = min_x [E_tri(x n) + E_SF((1-x) n)]"""
import numpy as np, json, sys
from bcs import solve, free, dos
def load(K):
    r = np.array(json.load(open(f"l2_energies_K{K}.json")))
    U = r[:, 0]; e3 = r[:, 2] + (r[:, 2] - r[:, 1]) * (1 / 2025) / (1 / 900 - 1 / 2025)
    e2 = r[:, 4] + (r[:, 4] - r[:, 3]) * (1 / 8100) / (1 / 3600 - 1 / 8100)
    e3u, e2u = r[:, 5], r[:, 6]
    c = K + 1 - 2 * np.sqrt(K)
    t3 = (e3 - e3u) / c; t2 = (e2 - e2u) / c
    return U, e3, e2, e3u, e2u, t3, t2
def E_tri(K, m, e3, t3, e, w):
    """energy per site of trions at density m: band e3 + t3*lambda, KM DOS in lambda."""
    if m <= 0: return 0.0
    order = np.argsort(e); es, ws = e[order], w[order]
    cum = np.cumsum(ws); i = np.searchsorted(cum, m)
    filled = np.sum(ws[:i] * es[:i]) + (m - (cum[i - 1] if i > 0 else 0)) * es[i]
    return m * e3 + t3 * filled
if __name__ == "__main__":
    K = int(sys.argv[1]); e, w = dos(K)
    U, e3, e2, e3u, e2u, t3, t2 = load(K)
    edge = -2 * np.sqrt(K)
    print(f"K={K}: U, E3_L2, E3_u, t3, e3 | E2_L2, E2_u, t2 | trion bound (E3_L2 < E2_u + edge)?")
    for i in range(len(U)):
        print(f"  {U[i]:.2f}: {e3[i]:.4f} {e3u[i]:.4f} t3={t3[i]:.4f} e3={e3[i]-edge*t3[i]:.4f} | {e2[i]:.4f} {e2u[i]:.4f} t2={t2[i]:.4f} | {e3[i] < min(3*edge, e2u[i] + edge) - 1e-6}")
    out = []
    ns = np.logspace(-4, np.log10(0.3), 40)
    xs = np.linspace(0, 1, 101)
    print("phase boundary: n1 (superfluid enters, x<1), n2 (trions gone, x=0)")
    for i in range(len(U)):
        if not e3[i] < min(3 * edge, e2u[i] + edge) - 1e-6: continue
        Esf_cache = {}
        def Esf(np_):
            if np_ <= 1e-9: return 0.0
            key = round(np_, 9)
            if key not in Esf_cache:
                mu, D, E = solve(K, U[i], np_, e, w); Esf_cache[key] = E + free(K, np_, e, w)[0]
            return Esf_cache[key]
        xstar = []
        for n in ns:
            Es = [E_tri(K, x * n, e3[i] - edge * t3[i], t3[i], e, w) + Esf((1 - x) * n) for x in xs]
            xstar.append(xs[int(np.argmin(Es))])
        xstar = np.array(xstar)
        n1 = ns[np.argmax(xstar < 0.999)] if np.any(xstar < 0.999) else np.nan
        n2 = ns[np.argmax(xstar < 0.001)] if np.any(xstar < 0.001) else np.nan
        out.append([U[i], n1, n2] + xstar.tolist())
        print(f"  U={U[i]:.2f}: n1={n1:.4g}  n2={n2:.4g}   x*(n) at n=1e-4,1e-3,1e-2,1e-1: {xstar[0]:.2f} {xstar[np.argmin(abs(ns-1e-3))]:.2f} {xstar[np.argmin(abs(ns-1e-2))]:.2f} {xstar[np.argmin(abs(ns-1e-1))]:.2f}", flush=True)
    json.dump({"ns": ns.tolist(), "rows": out}, open(f"phase_K{K}.json", "w"))

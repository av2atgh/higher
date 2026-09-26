"""Independent determination of t_3/rho_3(0), the combination that governs the composite's response to a site
potential, from the three-body root-fixed sector itself: for U above the clean square-summable threshold the trimer
is a normalisable composite with a delocalised centre of mass; a single site of energy eps0 binds that centre of mass
when |eps0| rho_3(0) > (K-1) t_3 / sqrt(K), i.e. at f* = |eps0*|/eps_c = t_3/rho_3(0). The onset is quadratic
(square-root band edge), so sqrt(-Delta E) is linear in f near f*: fit and extrapolate. L fixed, kinetic part built once."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, json, sys, time
sys.path.insert(0, "../trimer")
import tree4
K = int(sys.argv[1]); L = int(sys.argv[2]); Us = [float(u) for u in sys.argv[3:]]
t0 = time.time(); H0, S = tree4.build(3, K, 0.0, 0.0, L)
nroot = np.array([sum(1 for x in s if len(x) == 0) for s in S], float)
npair = np.array([sum(1 for a in range(3) for b in range(a + 1, 3) if s[a] == s[b]) for s in S], float)
print(f"K={K} L={L}: {len(S)} states ({time.time()-t0:.0f}s)", flush=True)
ec = tree4.impurity_threshold(K)
def E(U, eps0):
    H = H0 + sp.diags(eps0 * nroot - U * npair)
    return spl.eigsh(H, k=1, which="SA", tol=1e-9, maxiter=50000)[0][0]
out = {}
for U in Us:
    E0 = E(U, 0.0); rows = []
    for f in (0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.25, 0.30):
        dE = E(U, -f * ec) - E0; rows.append((f, dE))
        print(f"K={K} U={U}: f={f:.2f} dE={dE:+.5f}", flush=True)
    fs = np.array([r[0] for r in rows]); dE = np.array([r[1] for r in rows])
    m = dE < -2e-3
    if m.sum() >= 3:
        y = np.sqrt(-dE[m]); a, b = np.polyfit(fs[m], y, 1); fstar = -b / a
    else: fstar = float("nan")
    out[U] = dict(E0=E0, rows=rows, fstar=fstar)
    print(f"K={K} U={U}: f* = {fstar:.4f}  (t3/rho0 from the sector definition: see localization.json)", flush=True)
    json.dump(out, open(f"composite_impurity_K{K}.json", "w"), indent=1)

"""Independent t_3 from the trimer band in the square-summable (root-fixed) sector: where the whole trimer band lies below
the three-body continuum (binding > 2 sqrt K t_3), the eigenvalues below the edge at truncation L are the box-quantised
centre-of-mass band, spanning [E_int - 2 sqrt K t_3, E_int + 2 sqrt K t_3]; its width, extrapolated in L, is 4 sqrt K t_3.
No uniform-sector energy enters."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, json, sys, time
sys.path.insert(0, "../trimer")
import tree4
K = int(sys.argv[1]); Ls = [int(x) for x in sys.argv[2].split(",")]; Us = [float(u) for u in sys.argv[3:]]
e2 = {round(r[0], 3): r[4] + (r[4] - r[3]) * (1 / 8100) / (1 / 3600 - 1 / 8100) for r in json.load(open(f"../trimer/l2_energies_K{K}.json"))}  # extrapolated E2_L2
out = {}
for L in Ls:
    t0 = time.time(); H0, S = tree4.build(3, K, 0.0, 0.0, L)
    npair = np.array([sum(1 for a in range(3) for b in range(a + 1, 3) if s[a] == s[b]) for s in S], float)
    print(f"K={K} L={L}: {len(S)} states ({time.time()-t0:.0f}s)", flush=True)
    for U in Us:
        H = H0 - U * sp.diags(npair)
        w = np.sort(spl.eigsh(H, k=40, which="SA", tol=1e-8, maxiter=100000)[0])
        edge = min(-6 * np.sqrt(K), e2[round(U, 3)] - 2 * np.sqrt(K))   # lowest continuum: three free, or bound pair + free particle
        below = w[w < edge - 1e-6]
        width = below[-1] - below[0] if len(below) else float("nan")
        out[f"{K},{L},{U}"] = dict(edge=float(edge), eigs=w.tolist(), n_below=int(len(below)), bottom=float(below[0]) if len(below) else None, top=float(below[-1]) if len(below) else None, width=float(width), t3=float(width / (4 * np.sqrt(K))))
        print(f"K={K} L={L} U={U}: {len(below)} states below the edge {edge:.3f}; band [{below[0]:.4f}, {below[-1]:.4f}] width {width:.4f} -> t3 = {width/(4*np.sqrt(K)):.4f}   (next state {w[len(below)] if len(below) < len(w) else float('nan'):.4f})", flush=True)
        json.dump(out, open(f"trimer_bandwidth_K{K}.json", "w"), indent=1)

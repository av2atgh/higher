"""Onset density n_1(U) of the colour superfluid (trimer Fermi level reaches the binding energy) and the
superfluid fraction above it; writes phase_onset_K{K}.json. Same construction as phase.py, analytic onset."""
import json, numpy as np, sys
from bcs import dos, solve, free
from phase import load, E_tri
K = int(sys.argv[1])
e, w = dos(K); order = np.argsort(e); es, ws = e[order], w[order]; cum = np.cumsum(ws)
U, e3, e2, e3u, e2u, t3, t2 = load(K); edge = -2 * np.sqrt(K)
out = []
for i in range(len(U)):
    thr = min(3 * edge, e2u[i] + edge); b = thr - e3[i]; lamF = edge + b / t3[i]
    n1 = float(np.interp(lamF, es, cum)) if lamF < -edge else 1.0
    fr = []
    if n1 < 0.3:
        for n in (0.01, 0.03, 0.1, 0.3):
            if n <= n1: fr.append(0.0); continue
            xs = np.arange(0.80, 1.0001, 0.0005); best = (1e9, 1.0)
            for x in xs:
                npd = (1 - x) * n
                Esf = 0.0 if npd < 1e-9 else solve(K, U[i], npd, e, w)[2] + U[i] * npd ** 2 + free(K, npd, e, w)[0]
                Et = E_tri(K, x * n, e3[i] - edge * t3[i], t3[i], e, w) + Esf
                if Et < best[0]: best = (Et, x)
            fr.append(1 - best[1])
    print(f"K={K} U={U[i]:.2f}: binding {b:.4f} n1 {n1:.5f} fractions {fr}", flush=True)
    out.append([float(U[i]), float(b), float(t3[i]), n1, fr])
json.dump(out, open(f"phase_onset_K{K}.json", "w"))

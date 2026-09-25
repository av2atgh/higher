import numpy as np, json, sys, time
from tree4 import E, impurity_threshold
from tree3 import pair_threshold
K = int(sys.argv[1]); L = int(sys.argv[2])
ec = impurity_threshold(K)
def thr(m, eps0, Lm, hi, tol=2e-3):
    edge = -2 * m * np.sqrt(K)
    lo = 0.0
    if E(m, K, lo, eps0, Lm) < edge - 1e-9: return 0.0
    while E(m, K, hi, eps0, Lm) >= edge - 1e-9: hi *= 1.5
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if E(m, K, mid, eps0, Lm) < edge - 1e-9: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)
out = []
hi2 = hi3 = 0.5
for frac in (0.8, 0.6, 0.9, 0.4, 0.2, 0.95, 0.1, 0.0):
    eps0 = -frac * ec
    t0 = time.time()
    u2 = thr(2, eps0, 2 * L, hi2 + 0.3); u3 = thr(3, eps0, L, hi3 + 0.3)
    hi2, hi3 = max(u2, 0.2), max(u3, 0.2)
    out.append([eps0, L, u2, u3])
    print(f"K={K} L={L} eps0={eps0:8.4f} ({frac:.2f} eps_c): U2={u2:.4f}  U3={u3:.4f}  U3/U2={u3/u2 if u2>0 else float('nan'):.4f}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(f"scan_imp_K{K}_L{L}.json", "w"))

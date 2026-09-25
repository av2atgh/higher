import numpy as np, json, sys, time
from bs import impurity
from tree4 import impurity_threshold
K = int(sys.argv[1]); L = int(sys.argv[2])
ec = impurity_threshold(K)
out = []
for frac in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99):
    eps0 = -frac * ec; t0 = time.time()
    u2 = impurity(2, K, eps0, 2 * L); u3 = impurity(3, K, eps0, L)
    out.append([eps0, frac, L, u2, u3])
    print(f"K={K} L={L} eps0={eps0:8.4f} ({frac:.2f} eps_c): U2={u2:.4f}  U3={u3:.4f}  U3/U2={u3/u2:.4f}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(f"scan_bs_K{K}_L{L}.json", "w"))

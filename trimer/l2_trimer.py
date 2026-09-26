import numpy as np, json, sys, time
from bs import impurity
K = int(sys.argv[1]); out = []
for L in (30, 40, 50, 60, 70):
    t0 = time.time(); u3 = impurity(3, K, 0.0, L); u2 = impurity(2, K, 0.0, 2 * L)
    out.append([L, u2, u3]); json.dump(out, open(f"l2_trimer_K{K}.json", "w"))
    print(f"K={K} L={L}: U2_L2(2L)={u2:.5f}  U3_L2={u3:.5f}  ({time.time()-t0:.0f}s)", flush=True)

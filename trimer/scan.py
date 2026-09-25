import numpy as np, json, sys, time
from tree3 import build, ground, pair_threshold, pair_energy

def E3(K, U, L, sector="su3"):
    H, S = build(K, U, L, sector=sector)
    return ground(H)[0]

def bisect_threshold(K, L, sector="su3", lo=0.0, hi=None, tol=2e-3):
    """Smallest U with E3 < true three-particle edge (truncated continuum lies above it)."""
    edge = -6 * np.sqrt(K)
    if hi is None: hi = 2 * pair_threshold(K) + 2
    assert E3(K, hi, L, sector) < edge - 1e-9
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if E3(K, mid, L, sector) < edge - 1e-9: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)

if __name__ == "__main__":
    out = {}
    for K in (2, 3):
        edge = -6 * np.sqrt(K); U2 = pair_threshold(K)
        print(f"\n=== K={K} (z={K+1}): pair threshold U2 = {U2:.5f}, 3-body edge {edge:.5f}", flush=True)
        res = {"U2": U2, "edge": edge, "thr": {}, "grid": []}
        for L in (60, 120, 180):
            t0 = time.time()
            u3 = bisect_threshold(K, L)
            res["thr"][L] = u3
            print(f"  L={L}: U3 <= {u3:.4f}   ({time.time()-t0:.0f}s)", flush=True)
        L = 120
        print(f"  U      E3(L={L})     edge3      E2+edge1     bound?   E3_anti21")
        for U in np.arange(0.5, 6.01, 0.5):
            e3 = E3(K, U, L)
            e2 = pair_energy(K, U) - 2 * np.sqrt(K)
            ea = E3(K, U, L, sector="anti21")
            res["grid"].append([U, e3, e2, ea])
            print(f"  {U:4.2f}  {e3:10.5f}  {edge:9.5f}  {e2:10.5f}   {'T' if e3 < min(edge, e2) - 1e-9 else '-'}     {ea:10.5f}", flush=True)
        out[K] = res
    json.dump(out, open("scan.json", "w"), indent=1, default=float)

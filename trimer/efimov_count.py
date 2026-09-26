"""Number of trimer states below the three-body edge at the pair threshold."""
import numpy as np, scipy.sparse.linalg as spl
import tree3, tree4
from tree3 import pair_threshold
from bs import impurity
for K in (2, 3):
    U2 = pair_threshold(K); edge = -6 * np.sqrt(K)
    for U in (U2, U2 + 0.05):
        for L in (60, 120, 180):
            H, S = tree3.build(K, U, L)
            w = np.sort(spl.eigsh(H, k=8, which="SA", tol=1e-9)[0])
            print(f"uniform K={K} U={U:.4f} L={L}: below edge {np.sum(w < edge - 1e-7)}  lowest {np.round(w[:4] - edge, 5)}", flush=True)
    for L in (40, 60):
        U2l = impurity(2, K, 0.0, 2 * L)
        for U in (U2l, U2l + 0.05):
            H, S = tree4.build(3, K, U, 0.0, L)
            w = np.sort(spl.eigsh(H, k=8, which="SA", tol=1e-9)[0])
            print(f"L2 K={K} U={U:.4f} (U2_L2 at 2L={U2l:.4f}) L={L}: below edge {np.sum(w < edge - 1e-7)}  lowest {np.round(w[:4] - edge, 5)}", flush=True)

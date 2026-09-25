"""Binding thresholds by Birman-Schwinger at the continuum edge.
H = H0 - U V with V = number of coincident pairs (diagonal, >= 0).
A bound state below E exists iff U * lambda_max[ V^{1/2} (H0 - E)^{-1} V^{1/2} ] >= 1,
so U_c = 1 / lambda_max at E = edge. In the truncated space H0 - edge is positive
definite (truncated continuum above the edge); the truncation error decreases with L."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import tree3, tree4

def threshold_from(H0, V, edge, tol=1e-6):
    A = (H0 - edge * sp.identity(H0.shape[0], format="csr")).tocsr()
    sup = np.flatnonzero(V > 0)
    sq = np.sqrt(V[sup])
    n = len(sup)
    # Jacobi preconditioner for CG
    dinv = 1.0 / A.diagonal()
    Mpre = spl.LinearOperator(A.shape, matvec=lambda x: dinv * x)
    def matvec(v):
        rhs = np.zeros(A.shape[0]); rhs[sup] = sq * v
        x, info = spl.cg(A, rhs, rtol=1e-9, maxiter=20000, M=Mpre)
        assert info == 0, info
        return sq * x[sup]
    B = spl.LinearOperator((n, n), matvec=matvec, dtype=float)
    lam = spl.eigsh(B, k=1, which="LA", tol=tol, maxiter=2000)[0][0]
    return 1.0 / lam

def clean_trimer(K, L):
    H, S = tree3.build(K, 0.0, L)
    V = np.array([sum(1 for d in s if d == 0) * (sum(1 for d in s if d == 0) - 1) // 2 for s in S], float)
    return threshold_from(H, V, -6 * np.sqrt(K))

def impurity(m, K, eps0, L):
    H, S = tree4.build(m, K, 0.0, eps0, L)
    V = np.array([sum(1 for a in range(m) for b in range(a + 1, m) if s[a] == s[b]) for s in S], float)
    return threshold_from(H, V, -2 * m * np.sqrt(K))

if __name__ == "__main__":
    import time
    for K in (2, 3):
        t0 = time.time()
        print(f"clean K={K}: BS U3(L=60)={clean_trimer(K, 60):.4f}  bisection {'1.3100' if K==2 else '2.0728'}  ({time.time()-t0:.0f}s)", flush=True)
    for K, eps0, u2b, u3b in ((2, -0.5657, 1.2612, 0.7586), (3, -0.9238, 1.9083, 1.0986)):
        t0 = time.time()
        u2 = impurity(2, K, eps0, 80); u3 = impurity(3, K, eps0, 40)
        print(f"impurity K={K} eps0={eps0}: BS U2={u2:.4f} U3={u3:.4f}   bisection U2={u2b} U3={u3b}  ({time.time()-t0:.0f}s)", flush=True)

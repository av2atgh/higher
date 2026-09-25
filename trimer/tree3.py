"""Pair and trimer binding on the clean Bethe lattice (branching K, degree z=K+1).

Hopping -t on every edge, on-site attraction -U for each pair of particles
on the same site. Particles are distinguishable (three flavours); the
2+1 sector (two identical fermions + one distinguishable) is also provided.

Automorphism-invariant sector (tree analogue of zero total momentum):
three points of a tree have a unique median m; the orbit of (x1,x2,x3)
under Aut(tree) is labelled by the depths (a,b,c) of the particles below m,
with positive-depth particles in distinct branches. Multiplicity of an
orbit per median: n(a,b,c) = z^(k) (z-1)^(S-k), k = #positive depths,
S = a+b+c, z^(k) the falling factorial. The reduced kinetic operator M is
symmetrised by g = sqrt(n) f, giving H(s,s') = -t sqrt(M(s,s') M(s',s)).

Continuum edges on the infinite tree: E1 = -2 sqrt(K) t per free particle.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from itertools import permutations


def pair_threshold(K, t=1.0):
    """Exact: U_2 = 2 t (K-1)/sqrt(K)."""
    return 2 * t * (K - 1) / np.sqrt(K)


def pair_energy(K, U, t=1.0):
    """Exact pair energy in the invariant sector, or the continuum edge if unbound.
    Relative motion = half line in d with hopping J=2t sqrt(K) (d>=1),
    J0=2t sqrt(K+1) between d=0 and d=1, site energy -U at d=0.
    U = J(l+1/l) - J0^2 l/J, E = -J(l+1/l), l<1 bound."""
    J = 2 * t * np.sqrt(K)
    J0 = 2 * t * np.sqrt(K + 1)
    edge = -2 * J
    if U <= pair_threshold(K, t):
        return edge
    # solve for l in (0,1): J l^2 - (U) l + J - J0^2 l^2 / J ... rearrange:
    # U = J l + J/l - J0^2 l / J  ->  (J - J0^2/J) l^2 - U l + J = 0
    a = J - J0 ** 2 / J
    disc = U * U - 4 * a * J
    roots = [(U + s * np.sqrt(disc)) / (2 * a) for s in (+1, -1)]
    l = [r for r in roots if 0 < r < 1]
    assert len(l) == 1, (roots, U)
    l = l[0]
    return -J * (l + 1 / l)


def states(L, z, nparticles=3):
    """All depth tuples with sum <= L and nonzero multiplicity."""
    out = []
    if nparticles == 3:
        for a in range(L + 1):
            for b in range(L + 1 - a):
                for c in range(L + 1 - a - b):
                    k = (a > 0) + (b > 0) + (c > 0)
                    if k > z:
                        continue
                    ff = 1
                    for j in range(k):
                        ff *= (z - j)
                    if ff == 0:
                        continue
                    out.append((a, b, c))
    else:
        raise ValueError
    return out


def moves(s, z):
    """Yield (s', weight) for one hop of any particle from orbit s (unsymmetrised)."""
    n = len(s)
    for p in range(n):
        d = s[p]
        others = [q for q in range(n) if q != p]
        if d >= 1:
            s2 = list(s); s2[p] = d + 1; yield tuple(s2), z - 1
            s2 = list(s); s2[p] = d - 1; yield tuple(s2), 1
        else:
            occupied = sum(1 for q in others if s[q] > 0)
            if z - occupied > 0:
                s2 = list(s); s2[p] = 1; yield tuple(s2), z - occupied
            for q in others:
                if s[q] > 0:
                    s2 = list(s)
                    s2[q] -= 1
                    for r in others:
                        if r != q:
                            s2[r] += 1
                    yield tuple(s2), 1


def build(K, U, L, t=1.0, sector="su3"):
    """Sparse symmetrised Hamiltonian in the invariant sector, truncated at a+b+c<=L.
    sector: 'su3' (three flavours, all pairs attract) or
            'anti21' (particles 1,2 identical fermions, only pairs (1,3),(2,3) attract)."""
    z = K + 1
    S = states(L, z)
    idx = {s: i for i, s in enumerate(S)}
    rows, cols, vals = [], [], []
    for s in S:
        i = idx[s]
        for s2, w in moves(s, z):
            j = idx.get(s2)
            if j is not None:
                rows.append(i); cols.append(j); vals.append(float(w))
    M = sp.csr_matrix((vals, (rows, cols)), shape=(len(S), len(S)))
    Hk = -t * np.sqrt(M.multiply(M.T).toarray()) if len(S) < 3000 else -t * M.multiply(M.T).sqrt()
    Hk = sp.csr_matrix(Hk)
    if sector == "su3":
        V = np.array([-U * (sum(1 for d in s if d == 0) * (sum(1 for d in s if d == 0) - 1) // 2) for s in S])
    elif sector == "anti21":
        V = np.array([-U * (((s[0] == 0) and (s[2] == 0)) + ((s[1] == 0) and (s[2] == 0))) for s in S], float)
    H = Hk + sp.diags(V)
    if sector == "anti21":
        # antisymmetric combination |a,b,c> - |b,a,c>, a<b
        keep = [i for i, s in enumerate(S) if s[0] < s[1]]
        swap = np.array([idx[(s[1], s[0], s[2])] for s in S])
        Hs = H[:, swap]  # H(s, swap(s'))
        HA = (H - Hs)[keep][:, keep]
        return HA.tocsr(), [S[i] for i in keep]
    return H.tocsr(), S


def ground(H, k=1):
    if H.shape[0] < 400:
        w = np.linalg.eigvalsh(H.toarray())
        return w[:k]
    w = spl.eigsh(H, k=k, which="SA", tol=1e-10, maxiter=20000)[0]
    return np.sort(w)


if __name__ == "__main__":
    import sys
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    for L in (20, 40, 80):
        H, S = build(K, 0.0, L)
        print(K, L, len(S), ground(H)[0], -6 * np.sqrt(K))

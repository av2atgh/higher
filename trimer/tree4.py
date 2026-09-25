"""Few mobile particles on the clean Bethe lattice (branching K, degree z=K+1)
with a single impurity of site energy eps0 at the root.

Orbits under the root-preserving automorphism group (independent relabelling
of the children of every vertex): canonical form = relabel child branches in
order of first use by the particle paths (x1, x2, ...). States are canonical
configurations with total depth <= L. Measure n(s) = product over Steiner
vertices v of the falling factorial deg(v)^(j_v), j_v = number of child
branches of v used, deg(root) = z, deg(v) = K otherwise. Kinetic term
symmetrised as -t sqrt(M(s,s') M(s',s)). Everything here is square-summable
on the infinite tree (radial sector), so bound states are ordinary ones.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def canonical(cfg):
    relabel = {}
    counter = {}
    out = []
    for path in cfg:
        newp = []
        for i, lab in enumerate(path):
            key = (path[:i], lab)
            if key not in relabel:
                c = counter.get(path[:i], 0)
                relabel[key] = c
                counter[path[:i]] = c + 1
            newp.append(relabel[key])
        out.append(tuple(newp))
    return tuple(out)


def measure(cfg, z, K):
    used = {}
    for path in cfg:
        for i, lab in enumerate(path):
            used.setdefault(path[:i], set()).add(lab)
    n = 1
    for v, labs in used.items():
        deg = z if len(v) == 0 else K
        for j in range(len(labs)):
            n *= (deg - j)
    return n


def moves(cfg, z, K):
    """Yield (canonical target, weight) for single hops of any particle."""
    m = len(cfg)
    for p in range(m):
        x = cfg[p]
        d = len(x)
        others = [cfg[q] for q in range(m) if q != p]
        if d >= 1:
            new = list(cfg); new[p] = x[:-1]
            yield canonical(tuple(new)), 1
        deg = z if d == 0 else K
        used = sorted({y[d] for y in others if len(y) > d and y[:d] == x})
        for lab in used:
            new = list(cfg); new[p] = x + (lab,)
            yield canonical(tuple(new)), 1
        fresh = deg - len(used)
        if fresh > 0:
            lab = 0
            while lab in used:
                lab += 1
            new = list(cfg); new[p] = x + (lab,)
            yield canonical(tuple(new)), fresh


def states(m, L, z, K):
    """All canonical configurations of m particles with total depth <= L, by BFS."""
    start = canonical(tuple(() for _ in range(m)))
    seen = {start}
    frontier = [start]
    while frontier:
        nxt = []
        for s in frontier:
            for s2, w in moves(s, z, K):
                if sum(len(x) for x in s2) <= L and s2 not in seen:
                    seen.add(s2); nxt.append(s2)
        frontier = nxt
    return sorted(seen)


def build(m, K, U, eps0, L, t=1.0):
    z = K + 1
    S = states(m, L, z, K)
    idx = {s: i for i, s in enumerate(S)}
    rows, cols, vals = [], [], []
    for s in S:
        i = idx[s]
        acc = {}
        for s2, w in moves(s, z, K):
            j = idx.get(s2)
            if j is not None:
                acc[j] = acc.get(j, 0) + w
        for j, w in acc.items():
            rows.append(i); cols.append(j); vals.append(float(w))
    M = sp.csr_matrix((vals, (rows, cols)), shape=(len(S), len(S)))
    Hk = -t * M.multiply(M.T).sqrt()
    V = np.zeros(len(S))
    for s in S:
        i = idx[s]
        V[i] += eps0 * sum(1 for x in s if len(x) == 0)
        for a in range(m):
            for b in range(a + 1, m):
                if s[a] == s[b]:
                    V[i] -= U
    return (Hk + sp.diags(V)).tocsr(), S


def ground(H):
    if H.shape[0] < 400:
        return np.linalg.eigvalsh(H.toarray())[0]
    return spl.eigsh(H, k=1, which="SA", tol=1e-10, maxiter=50000)[0][0]


def E(m, K, U, eps0, L):
    H, S = build(m, K, U, eps0, L)
    return ground(H)


def impurity_threshold(K, t=1.0):
    """Single particle binds to the root impurity when eps0 < -(K-1) t / sqrt(K)."""
    return (K - 1) * t / np.sqrt(K)


if __name__ == "__main__":
    # checks
    for K in (2, 3):
        z = K + 1
        for m in (2, 3):
            S = states(m, 16, z, K)
            bad_db = bad_w = 0
            for s in S:
                fwd = {}
                tot = 0
                for s2, w in moves(s, z, K):
                    fwd[s2] = fwd.get(s2, 0) + w; tot += w
                if tot != m * z: bad_w += 1
                for s2, w in fwd.items():
                    if sum(len(x) for x in s2) > 16: continue
                    back = sum(w2 for s3, w2 in moves(s2, z, K) if s3 == s)
                    if measure(s, z, K) * w != measure(s2, z, K) * back: bad_db += 1
            print(f"K={K} m={m}: {len(S)} states, out-weight errors {bad_w}, detailed-balance errors {bad_db}")
        # single particle: impurity threshold and edge
        ec = impurity_threshold(K)
        for eps0 in (-0.5 * ec, -ec * 1.001, -1.5 * ec):
            e = E(1, K, 0.0, eps0, 400)
            print(f"K={K} m=1 eps0={eps0:.4f} (eps_c={ec:.4f}): E={e:.6f} edge={-2*np.sqrt(K):.6f}")

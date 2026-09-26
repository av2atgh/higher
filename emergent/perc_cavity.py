"""Percolation of emergent hyperedges on the Bethe lattice (branching K).

Borromean centres are iid with probability p; each hosts a localised trimer whose
support is the ball of radius r around it. A site is *activated* when it lies
within distance r of a centre, and the emergent hyperedges percolate when the
activated sites form an infinite connected subgraph of the tree (two supports
overlap or touch). Exact cavity on the tree, in the style of the book's dependent
layer: the message from a child c to its parent is the law of
  D_c = distance from c to the nearest centre inside the subtree T_c (0..r, or M=r+1 for none),
and the conditional probability Q[u][d] that c is activated and connected to
infinity inside T_c, given D_c = d and U_c = u, the distance from c to the nearest
centre outside T_c (1..M). Activation of a site below c depends on centres above
it only through U, which is why the message is conditioned on u.
r = 0 is ordinary site percolation (p_c = 1/K), the check.
"""
import numpy as np, itertools, sys, json


def fixed_point(K, r, p, iters=5000, tol=1e-12):
    M = r + 1
    combos = list(itertools.product(range(M + 1), repeat=K))
    # law of D_c is independent of Q: iterate it first
    P = np.zeros(M + 1); P[0] = p; P[M] = 1 - p
    for _ in range(M + 2):   # the law of D converges in r+1 steps; renormalise (the map amplifies rounding by K per step)
        Pn = np.zeros(M + 1)
        for ds in combos:
            w = float(np.prod([P[d] for d in ds]))
            Pn[0] += w * p
            Pn[min(M, 1 + min(ds))] += w * (1 - p)
        P = Pn / Pn.sum()
    Q = np.ones((M + 1, M + 1))
    for it in range(iters):
        acc = np.zeros((M + 1, M + 1))
        for ds in combos:
            w = float(np.prod([P[d] for d in ds]))
            if w == 0.0:
                continue
            for chi, wc in ((1, p), (0, 1 - p)):
                d_c = 0 if chi else min(M, 1 + min(ds))
                for u in range(1, M + 1):
                    if min(d_c, u) >= M:
                        continue
                    prob_none = 1.0
                    for k in range(K):
                        others = [ds[l] for l in range(K) if l != k]
                        dmk = 0 if chi else (min(M, 1 + min(others)) if others else M)
                        uk = min(M, 1 + min(u, dmk))
                        prob_none *= (1.0 - Q[uk][ds[k]])
                    acc[u][d_c] += w * wc * (1.0 - prob_none)
        Qn = np.zeros((M + 1, M + 1))
        for d in range(M + 1):
            if P[d] > 0:
                Qn[:, d] = acc[:, d] / P[d]
        Qn[0, :] = 0.0
        Qn = np.clip(Qn, 0.0, 1.0)
        delta = np.abs(Qn - Q).max()
        Q = Qn
        if delta < tol:
            break
    return P, Q


def order_parameter(K, r, p, **kw):
    """Probability that a centre belongs to the infinite activated cluster (root has K+1 children, D^{(-k)} = 0 -> U_k = 1)."""
    M = r + 1
    P, Q = fixed_point(K, r, p, **kw)
    S = 0.0
    for ds in itertools.product(range(M + 1), repeat=K + 1):
        w = np.prod([P[d] for d in ds])
        S += w * (1 - np.prod([1 - Q[1][d] for d in ds]))
    return S


def threshold(K, r, lo=1e-4, hi=1.0, iters=30):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if order_parameter(K, r, mid) > 1e-6: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    out = {}
    for K in (2, 3):
        for r in (0, 1, 2, 3):
            pc = threshold(K, r); out[f"{K},{r}"] = pc
            print(f"K={K} r={r}: p_c = {pc:.5f}" + ("  (site percolation 1/K = %.5f)" % (1 / K) if r == 0 else ""), flush=True)
    json.dump(out, open("perc_thresholds.json", "w"), indent=1)
    # order parameter curves for the MC check
    for K, r in ((2, 1), (3, 1), (2, 2)):
        print(f"K={K} r={r}: S(p):", " ".join(f"{p:.3f}->{order_parameter(K, r, p):.4f}" for p in (0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3)), flush=True)


# ---------------------------------------------------------------- heterogeneous radii
def fixed_point_radii(K, p, radii, iters=5000, tol=1e-12):
    """Centres iid(p); a centre has radius r with probability radii[r]. State d = min over centres
    of (distance - radius): values -R..0 (activating) and 1 (none within reach; cannot activate later)."""
    R = max(radii); off = R  # index = d + R, d in -R..1
    nv = R + 2
    def cap(d): return min(d, 1)
    P = np.zeros(nv); P[1 + off] = 1 - p
    for r_, pr in radii.items(): P[-r_ + off] += p * pr
    P0 = P.copy()
    combos = list(itertools.product(range(nv), repeat=K))
    dvals = np.arange(nv) - off
    for _ in range(R + 3):
        Pn = np.zeros(nv)
        for ds in combos:
            w = float(np.prod([P[i] for i in ds]))
            Pn += w * P0 * 0  # placeholder
            # own centre with radius r_: d_c = -r_ ; no centre: d_c = 1 + min children
            for r_, pr in radii.items(): Pn[-r_ + off] += w * p * pr
            Pn[cap(1 + min(dvals[list(ds)])) + off] += w * (1 - p)
        P = Pn / Pn.sum()
    Q = np.ones((nv, nv))  # Q[u+off][d+off], u in -R+1..1
    for it in range(iters):
        acc = np.zeros((nv, nv))
        for ds in combos:
            w = float(np.prod([P[i] for i in ds]))
            if w == 0.0: continue
            dk = dvals[list(ds)]
            owns = [(-r_, p * pr) for r_, pr in radii.items()] + [(None, 1 - p)]
            for own, wc in owns:
                d_c = own if own is not None else cap(1 + min(dk))
                for ui in range(nv):
                    u = dvals[ui]
                    if u == -R and False: pass
                    if min(d_c, u) > 0: continue
                    prob_none = 1.0
                    for k in range(K):
                        others = [dk[l] for l in range(K) if l != k]
                        dmk = own if own is not None else (cap(1 + min(others)) if others else 1)
                        uk = cap(1 + min(u, dmk))
                        prob_none *= (1.0 - Q[uk + off][ds[k]])
                    acc[ui][d_c + off] += w * wc * (1.0 - prob_none)
        Qn = np.zeros((nv, nv))
        for i in range(nv):
            if P[i] > 0: Qn[:, i] = acc[:, i] / P[i]
        Qn = np.clip(Qn, 0.0, 1.0)
        delta = np.abs(Qn - Q).max(); Q = Qn
        if delta < tol: break
    return P, Q, off, dvals


def order_parameter_radii(K, p, radii, **kw):
    """fraction of centres in the infinite activated cluster (root = centre of radius r: children see u = 1 - r)."""
    P, Q, off, dvals = fixed_point_radii(K, p, radii, **kw)
    nv = len(P); S = 0.0
    for r_, pr in radii.items():
        u = min(1 - r_, 1)
        for ds in itertools.product(range(nv), repeat=K + 1):
            w = float(np.prod([P[i] for i in ds]))
            S += pr * w * (1 - np.prod([1 - Q[u + off][i] for i in ds]))
    return S


def threshold_radii(K, radii, lo=1e-5, hi=1.0, iters=25):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if order_parameter_radii(K, mid, radii) > 1e-6: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)

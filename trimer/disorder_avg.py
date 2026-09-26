"""Single-site (Lifshitz-tail) approximation for box disorder eps_i in [-W/2, W/2]:
a site of energy eps is classified by the clean-environment thresholds U_2(eps), U_3(eps)
(scan_bs, L=60) and eps_c. Fractions of sites: hold one particle (|eps|>eps_c),
hold a trimer but no pair (Borromean), hold a pair (and trimer), hold nothing."""
import numpy as np, json
from tree4 import impurity_threshold
def curves(K):
    r = np.array(json.load(open(f"scan_bs_K{K}_L60.json")))  # eps0, frac, L, u2, u3
    fr, u2, u3 = r[:, 1], r[:, 3], r[:, 4]
    fr = np.append(fr, 1.0); u2 = np.append(u2, 0.0); u3 = np.append(u3, 0.0)
    return fr, u2, u3
def fractions(K, U, W, n=20001):
    ec = impurity_threshold(K); fr, u2, u3 = curves(K)
    eps = np.linspace(-W / 2, W / 2, n)
    f = np.clip(-eps / ec, 0, None)
    one = f > 1
    U2 = np.interp(f, fr, u2); U3 = np.interp(f, fr, u3)
    bor = (~one) & (eps < 0) & (U3 < U) & (U < U2)
    pair = (~one) & (eps < 0) & (U >= U2)
    return one.mean(), bor.mean(), pair.mean()
if __name__ == "__main__":
    for K in (2, 3):
        ec = impurity_threshold(K); fr, u2, u3 = curves(K)
        print(f"K={K}: eps_c={ec:.4f}  U2_L2(0)={u2[0]:.4f} U3_L2(0)={u3[0]:.4f}")
        print("   W/eps_c   U/U2(0)   one   Borromean   pair")
        for Wf in (0.5, 1.0, 2.0, 4.0):
            for Uf in (0.5, 0.7, 0.8, 0.9, 1.0, 1.1):
                one, bor, pair = fractions(K, Uf * u2[0], Wf * ec)
                print(f"   {Wf:5.2f}    {Uf:5.2f}   {one:.3f}   {bor:.3f}     {pair:.3f}")

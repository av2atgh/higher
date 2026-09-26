import numpy as np, json
from scipy.optimize import curve_fit
def fit(L, U):
    f = lambda L, a, c, p: a + c * L ** (-p)
    p0 = [U[-1] - 0.01, 30.0, 2.0]
    (a, c, p), cov = curve_fit(f, L, U, p0=p0, maxfev=20000)
    # bracket with fixed p = 1 and 2 on the last two points
    a1 = U[-1] - (U[-2] - U[-1]) / (L[-2] / L[-1] - 1) * 1.0  # p=1: U(L)=a+c/L
    c1 = (U[-2] - U[-1]) / (1 / L[-2] - 1 / L[-1]); a1 = U[-1] - c1 / L[-1]
    c2 = (U[-2] - U[-1]) / (1 / L[-2] ** 2 - 1 / L[-1] ** 2); a2 = U[-1] - c2 / L[-1] ** 2
    return a, p, a1, a2
exact = {2: 1.965677, 3: 3.251185, 4: 4.268903, 5: 5.136053, 6: 5.903985, 10: 8.402470}
print("pair (validation of the extrapolation against the exact kernel value):")
res = {}
for K in (2, 3, 4, 5, 6, 10):
    r = np.array(json.load(open(f"l2_trimer_K{K}.json")))
    L, u2, u3 = r[:, 0], r[:, 1], r[:, 2]
    a, p, a1, a2 = fit(2 * L, u2)
    print(f"  K={K}: fit {a:.4f} (p={p:.2f})  [p=1: {a1:.4f}, p=2: {a2:.4f}]  exact {exact[K]:.4f}")
print("trimer:")
for K in (2, 3, 4, 5, 6, 10):
    r = np.array(json.load(open(f"l2_trimer_K{K}.json")))
    L, u2, u3 = r[:, 0], r[:, 1], r[:, 2]
    a, p, a1, a2 = fit(L, u3)
    res[K] = [a, p, a1, a2]
    print(f"  K={K}: U3 fit {a:.4f} (p={p:.2f})  [p=1: {a1:.4f}, p=2: {a2:.4f}]  L=60 value {u3[-1]:.4f}  ratio {a/exact[K]:.3f}  window {100*(1-a/exact[K]):.1f}%")
json.dump(res, open("l2_extrap.json", "w"))

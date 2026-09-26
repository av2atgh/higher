"""Percolation of the emergent hyperedges in the (U, W) plane, K = 2, 3.
Sites are Borromean (host a trimer and no pair) with probability f_B(U, W) from the isolated-site
classification (trimer/disorder_avg.py); the trimer at a site of energy eps0 has support radius
r90(eps0, U) from site_density_K{K}.json (interpolated on the (f = -eps0/eps_c, x = (U-U3)/(U2-U3)) grid,
capped at RMAX). S = fraction of Borromean sites in the infinite cluster of overlapping/touching supports,
from the heterogeneous-radius cavity (perc_cavity.fixed_point_radii). Two variants: all Borromean sites,
and compact ones only (r90 <= 3), the range in which the isolated-site picture is safe."""
import numpy as np, json, sys, os, time
from scipy.interpolate import RegularGridInterpolator
from perc_cavity import order_parameter_radii
sys.path.insert(0, "../trimer")
from tree4 import impurity_threshold
RMAX = 6
def radius_interp(K):
    rows = json.load(open(f"site_density_K{K}.json"))
    fs = sorted({r["f"] for r in rows}); xs = sorted({r["x"] for r in rows})
    grid = np.array([[next(r["r90"] for r in rows if r["f"] == f and r["x"] == x) for x in xs] for f in fs], float)
    return RegularGridInterpolator((fs, xs), grid, bounds_error=False, fill_value=None)
def curves(K):
    r = np.array(json.load(open(f"../trimer/scan_bs_K{K}_L60.json")))
    fr, u2, u3 = r[:, 1], r[:, 3], r[:, 4]
    return np.append(fr, 1.0), np.append(u2, 0.0), np.append(u3, 0.0)
def borromean_radii(K, U, W, interp, n=4001):
    ec = impurity_threshold(K); fr, u2, u3 = curves(K)
    eps = np.linspace(-W / 2, W / 2, n); f = np.clip(-eps / ec, 0, None)
    U2 = np.interp(f, fr, u2); U3 = np.interp(f, fr, u3)
    bor = (f <= 1) & (eps < 0) & (U3 < U) & (U < U2)
    fB = bor.mean()
    if fB == 0: return 0.0, {}
    fb, xb = f[bor], (U - U3[bor]) / (U2[bor] - U3[bor])
    r = np.rint(interp(np.stack([np.clip(fb, 0.1, 0.95), np.clip(xb, 0.1, 0.9)], 1))).astype(int)
    r = np.clip(r, 0, RMAX)
    vals, cnt = np.unique(r, return_counts=True)
    return fB, {int(v): c / len(r) for v, c in zip(vals, cnt)}
if __name__ == "__main__":
    K = int(sys.argv[1]); interp = radius_interp(K); ec = impurity_threshold(K); fr, u2, u3 = curves(K)
    out = []
    for Wf in (0.5, 1.0, 2.0, 4.0):
        for Uf in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
            t0 = time.time(); U = Uf * u2[0]; W = Wf * ec
            fB, radii = borromean_radii(K, U, W, interp)
            S = order_parameter_radii(K, fB, radii) if fB > 0 else 0.0
            compact = {r: p for r, p in radii.items() if r <= 3}; fc = fB * sum(compact.values())
            Sc = order_parameter_radii(K, fc, {r: p / sum(compact.values()) for r, p in compact.items()}) if fc > 0 else 0.0
            mean_r = sum(r * p for r, p in radii.items()) if radii else float("nan")
            out.append(dict(K=K, Wf=Wf, Uf=Uf, U=U, W=W, fB=fB, radii=radii, mean_r=mean_r, S=S, f_compact=fc, S_compact=Sc))
            print(f"K={K} W={Wf}ec U={Uf}U2: f_B={fB:.3f} <r90>={mean_r:.2f} radii={radii}  S={S:.3f} | compact f={fc:.3f} S_c={Sc:.3f}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump(out, open(f"perc_map_K{K}.json", "w"))

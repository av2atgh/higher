"""Anderson localisation of the composites, from the book's cavity results (statmech package, Ch. 18).
A composite with hopping t_m and one-body density rho sees site disorder with variance factor
s_m^2 = sum_d rho(d)^2/n_d (first order), so its Anderson problem is the single-particle one with
hop = t_m and W_eff = s_m W: critical disorder W_c^(m) = W_c^(1) t_m / s_m, and the mobility edge
E_c^(m)(W) = E_3 + t_m [E_c^(1)(s_m W / t_m) - band bottom]. Inputs: uniform_density.json (s_m),
l2_energies (t_m via phase.load), anderson_wc.csv and anderson_lines.csv (book, P = 1e6)."""
import numpy as np, json, csv, sys
sys.path.insert(0, "../trimer")
import os; _cwd = os.getcwd(); os.chdir("../trimer"); from phase import load; os.chdir(_cwd)
BOOK = "/Users/avazquez/av2atg/chygraph/statmech/probe/results/"
wc = {(int(r["s"]), int(r["t"])): float(r["Wc"]) for r in csv.DictReader(open(BOOK + "anderson_wc.csv"))}
lines = [r for r in csv.DictReader(open(BOOK + "anderson_lines.csv"))]
def tail(s, t, Wp):
    """E_c^loc - band bottom (hop = 1) at disorder Wp, interpolated in the book's table."""
    rows = sorted([(float(r["W"]), float(r["Ec_loc"]) - float(r["band_bottom"])) for r in lines if int(r["s"]) == s and int(r["t"]) == t])
    Ws, Es = zip(*rows)
    if Wp < Ws[0]: return Es[0] * Wp / Ws[0]
    if Wp > Ws[-1]: return float("nan")
    return float(np.interp(Wp, Ws, Es))
dens = json.load(open("uniform_density.json"))
out = {}
for K in (2, 3):
    st = (K + 1, 0); Wc1 = wc[st]
    os.chdir("../trimer"); U, e3, e2, e3u, e2u, t3, t2 = load(K); os.chdir(_cwd)
    print(f"K={K}: single-particle W_c = {Wc1:.2f} t (book, ({st[0]},0))")
    print("   U      t3     s3    W_c^(3)   t2     s2    W_c^(2)  | tail width at W = 0.3t (trimer units of t): E_c - bottom")
    rows = []
    for i in range(len(U)):
        key = f"{K},{U[i]}"
        if key not in dens: continue
        s3 = dens[key]["s3"]; s2 = dens[key]["s2"]
        Wc3 = Wc1 * t3[i] / s3; Wc2 = Wc1 * t2[i] / s2 if t2[i] > 1e-6 else float("nan")
        W = 0.3; Wp = s3 * W / t3[i]; tw = tail(*st, Wp) * t3[i]
        rows.append(dict(U=float(U[i]), t3=float(t3[i]), s3=s3, Wc3=Wc3, t2=float(t2[i]), s2=s2, Wc2=Wc2, tail_W0p3=tw, Wp_W0p3=Wp))
        print(f"  {U[i]:.2f}  {t3[i]:.3f}  {s3:.3f}  {Wc3:6.3f}   {t2[i]:.3f}  {s2:.3f}  {Wc2:6.3f}  |  W'={Wp:6.2f}  tail={tw:+.4f}")
    out[K] = rows
json.dump(out, open("localization.json", "w"), indent=1)

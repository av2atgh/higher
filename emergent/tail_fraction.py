"""Integrated density of states below the mobility edge on the (K+1,0) Bethe ensemble at reduced disorder W' (hop = 1):
n_tail(W') = fraction of single-particle states localised in the Lifshitz tail. E_c from the book's table (anderson_lines.csv,
interpolated in W'); the integrated DOS from exact diagonalisation of random regular instances built with statmech.resolvent
(N sites, several samples). For a composite with hopping t_m and disorder factor s_m, W' = s_m W / t_m, and a Fermi liquid of
the composites at density n (per site) is fully localised when n < n_tail(W')."""
import sys, json, csv, numpy as np
sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/statmech/src"); sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/percolation/src")
from statmech import resolvent as R
BOOK = "/Users/avazquez/av2atg/chygraph/statmech/probe/results/"
lines = [r for r in csv.DictReader(open(BOOK + "anderson_lines.csv"))]
def Ec(s, Wp):
    rows = sorted([(float(r["W"]), float(r["Ec_loc"])) for r in lines if int(r["s"]) == s and int(r["t"]) == 0])
    Ws, Es = zip(*rows); bb = float([r for r in lines if int(r["s"]) == s][0]["band_bottom"])
    if Wp < Ws[0]: return bb + (Es[0] - bb) * Wp / Ws[0]
    return float(np.interp(Wp, Ws, Es))
out = {}
for K in (2, 3):
    s = K + 1; rng = np.random.default_rng(0); N = 3000
    print(f"K={K} ((s,t)=({s},0)), N={N}, 4 samples:  W'   E_c   n_tail (fraction of states below E_c)")
    for Wp in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0):
        ec = Ec(s, Wp); fr = []
        for smp in range(4):
            complexes = R.regular_instance(N, [2], [s], rng); eps = rng.uniform(-Wp / 2, Wp / 2, N)
            H = R.hamiltonian(N, complexes, eps); w = np.linalg.eigvalsh(np.asarray(H.real)); fr.append(np.mean(w < ec))
        out[f"{K},{Wp}"] = dict(Ec=ec, n_tail=float(np.mean(fr)), err=float(np.std(fr) / 2))
        print(f"   {Wp:5.1f}  {ec:8.4f}   {np.mean(fr):.5f} +- {np.std(fr)/2:.5f}", flush=True)
json.dump(out, open("tail_fraction.json", "w"), indent=1)

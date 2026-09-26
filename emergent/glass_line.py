"""Tail fractions at the extended mobility edges (mobility_ext.json) and the glass line:
W_glass(U) = smallest disorder at which every state of the trimer liquid (n < n_1(U)) is localised, i.e. n_tail(s3 W / t3) >= n_1(U);
also n_tail at W = 0.5 t and 1 t, and W_c^(3). Instances built with statmech.resolvent (N = 3000, 4 samples)."""
import sys, json, csv, numpy as np
sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/statmech/src"); sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/percolation/src")
from statmech import resolvent as R
tail = json.load(open("tail_fraction.json")); ext = json.load(open("mobility_ext.json"))
rng = np.random.default_rng(1); N = 3000
for key, v in ext.items():
    s, Wp = key.split(","); K = int(s) - 1; Wp = float(Wp)
    if f"{K},{Wp}" in tail: continue
    fr = []
    for smp in range(4):
        complexes = R.regular_instance(N, [2], [int(s)], rng); eps = rng.uniform(-Wp / 2, Wp / 2, N)
        H = R.hamiltonian(N, complexes, eps); w = np.linalg.eigvalsh(np.asarray(H.real)); fr.append(np.mean(w < v["Ec"]))
    tail[f"{K},{Wp}"] = dict(Ec=v["Ec"], n_tail=float(np.mean(fr)), err=float(np.std(fr) / 2))
    print(f"K={K} W'={Wp}: E_c={v['Ec']:.4f} n_tail={np.mean(fr):.4f}", flush=True)
json.dump(tail, open("tail_fraction.json", "w"), indent=1)
loc = json.load(open("localization.json")); wc = {2: 17.39, 3: 32.26}
onset = {K: {round(r[0], 3): r[3] for r in json.load(open(f"../trimer/phase_onset_K{K}.json"))} for K in (2, 3)}
out = {}
for K in (2, 3):
    Ws = sorted(float(k.split(",")[1]) for k in tail if k.startswith(f"{K},")); nt = [tail[f"{K},{w}"]["n_tail"] for w in Ws]
    Ws.append(wc[K]); nt.append(1.0)   # every state localised at W_c
    print(f"K={K}: U   n_1   n_tail(W=0.5t)  n_tail(W=1t)  W_glass  W_c^(3)")
    rows = []
    for r in loc[str(K)]:
        ratio = r["s3"] / r["t3"]; n1 = onset[K].get(round(r["U"], 3), float("nan"))
        ntail = lambda W: float(np.interp(ratio * W, Ws, nt))
        Wg = float("nan")
        if n1 == n1 and n1 < 1:
            grid = np.linspace(0.02, 5, 500); idx = next((i for i, W in enumerate(grid) if ntail(W) >= n1), None)
            Wg = float(grid[idx]) if idx is not None else float("nan")
        rows.append(dict(U=r["U"], n1=n1, nt05=ntail(0.5), nt1=ntail(1.0), W_glass=Wg, Wc3=r["Wc3"]))
        print(f"  {r['U']:.2f}  {n1:.3f}   {ntail(0.5):.4f}   {ntail(1.0):.4f}   {Wg:.2f}   {r['Wc3']:.2f}")
    out[K] = rows
json.dump(out, open("glass_line.json", "w"), indent=1)

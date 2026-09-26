"""Resume l2_energies.py: skip U values already in l2_energies_K{K}.json, append the rest."""
import numpy as np, json, sys, time, os
from tree4 import E
import tree3
from tree3 import pair_energy
K = int(sys.argv[1]); fn = f"l2_energies_K{K}.json"
out = json.load(open(fn)) if os.path.exists(fn) else []
done = {round(r[0], 6) for r in out}
Us = [float(u) for u in sys.argv[2:]] if len(sys.argv) > 2 else ([1.4, 1.5, 1.7, 2.0, 2.5, 3.0, 4.0] if K == 2 else [2.2, 2.4, 2.7, 3.0, 3.5, 4.0, 5.0])
for U in Us:
    if round(U, 6) in done: continue
    t0 = time.time()
    e3 = {L: E(3, K, U, 0.0, L) for L in (30, 45)}
    e2 = {L: E(2, K, U, 0.0, 2 * L) for L in (30, 45)}
    H, S = tree3.build(K, U, 60); e3u = tree3.ground(H)[0]
    out.append([U, e3[30], e3[45], e2[30], e2[45], e3u, pair_energy(K, U)])
    out.sort(key=lambda r: r[0])
    print(f"K={K} U={U:.2f}: E3_L2 {e3[30]:.5f} {e3[45]:.5f}  E2_L2 {e2[30]:.5f} {e2[45]:.5f}  E3_u {e3u:.5f}  E2_u {pair_energy(K, U):.5f}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump(out, open(fn, "w"))

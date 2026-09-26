"""Direct mobility-edge run with the statmech package at the trimer's reduced disorder W' = s3 W / t3 (hop = 1),
to compare with the interpolation of the book's table. K=2 (3,0), U = 1.7t, W = 0.3t -> W' = 1.83*0.3/0.197 = 2.79."""
import sys, time
sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/statmech/probe"); sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/statmech/src"); sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/percolation/src")
from anderson import Ensemble
import numpy as np
Wp = float(sys.argv[1]) if len(sys.argv) > 1 else 2.79
t0 = time.time(); e = Ensemble(3, 0, Wp, P=200000, seed=1)
Emin = e.E_min(); bb = e.band_bottom()
Ec = e.mobility_edge(Emin + 0.02, bb + 0.5, iters=8, betas=(0.5,), sweeps=300, transient=150)
print(f"(3,0) W'={Wp}: E_min={Emin:.4f} band bottom={bb:.4f} E_c^loc={Ec:.4f}  tail width E_c-bottom={Ec-bb:+.4f}  ({time.time()-t0:.0f}s)", flush=True)

"""Mobility edge at the band bottom beyond the book's table, with the statmech package: (s,0) ensembles at reduced disorder W'."""
import sys, time, json
sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/statmech/probe"); sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/statmech/src"); sys.path.insert(0, "/Users/avazquez/av2atg/chygraph/percolation/src")
from anderson import Ensemble
out = {}
for s, Wps in ((3, (14.0, 16.0)), (4, (16.0, 20.0, 24.0, 28.0))):
    for Wp in Wps:
        t0 = time.time(); e = Ensemble(s, 0, Wp, P=200000, seed=1)
        Emin = e.E_min(); bb = e.band_bottom()
        Ec = e.mobility_edge(Emin + 0.02, 0.0, iters=9, betas=(0.5,), sweeps=300, transient=150)
        out[f"{s},{Wp}"] = dict(Emin=Emin, band_bottom=bb, Ec=Ec)
        print(f"({s},0) W'={Wp}: E_min={Emin:.4f} bottom={bb:.4f} E_c={Ec:.4f} ({time.time()-t0:.0f}s)", flush=True)
        json.dump(out, open("mobility_ext.json", "w"), indent=1)

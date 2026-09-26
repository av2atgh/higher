"""rrg_q2 on random 3-regular graphs conditioned on girth >= g (no short cycles), by rejection."""
import numpy as np, networkx as nx, sys, json, time
from rrg_q2 import run
def girth(G):
    best = 10 ** 9
    for s in G.nodes:
        dist = {s: 0}; parent = {s: None}; q = [s]
        while q:
            nq = []
            for v in q:
                for w in G[v]:
                    if w not in dist:
                        dist[w] = dist[v] + 1; parent[w] = v; nq.append(w)
                    elif parent[v] != w:
                        best = min(best, dist[v] + dist[w] + 1)
            q = nq
    return best
if __name__ == "__main__":
    N = int(sys.argv[1]); g = int(sys.argv[2]); nseeds = int(sys.argv[3])
    found = []; seed = 0; t0 = time.time()
    while len(found) < nseeds:
        seed += 1
        G = nx.random_regular_graph(3, N, seed=seed)
        if girth(G) >= g: found.append(seed)
    print(f"N={N}: girth>={g} seeds {found} after {seed} samples ({time.time()-t0:.0f}s)", flush=True)
    out = []
    for s in found:
        out += run(N, s, (0.0, 1.0, 1.4, 1.7, 2.0, 2.3)); json.dump(out, open(f"rrg_girth{g}_N{N}.json", "w"))

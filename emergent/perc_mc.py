"""Monte Carlo check of perc_cavity on random regular graphs: centres iid(p), activated = within r of a centre,
largest component of the activated subgraph; fraction of centres in it."""
import numpy as np, networkx as nx, sys
def run(K, r, p, N=100000, seed=0):
    rng = np.random.default_rng(seed); G = nx.random_regular_graph(K + 1, N, seed=seed)
    centre = rng.random(N) < p
    act = centre.copy(); frontier = set(np.flatnonzero(centre))
    for _ in range(r):
        nxt = set()
        for v in frontier:
            for w in G[v]:
                if not act[w]: act[w] = True; nxt.add(w)
        frontier = nxt
    H = G.subgraph(np.flatnonzero(act))
    big = max(nx.connected_components(H), key=len) if H.number_of_nodes() else set()
    return sum(centre[v] for v in big) / max(centre.sum(), 1)
if __name__ == "__main__":
    for K, r in ((2, 1), (3, 1), (2, 2)):
        print(f"K={K} r={r}: S_MC(p):", " ".join(f"{p:.3f}->{np.mean([run(K, r, p, seed=s) for s in range(2)]):.4f}" for p in (0.05, 0.1, 0.15, 0.2, 0.3, 0.4)), flush=True)

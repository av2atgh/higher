import numpy as np, networkx as nx
from tree3 import build, ground
from checks import ed_graph
# Cayley tree z=4 (K=3), depth 4, N=161: no loops, boundary at distance 4 from root
def cayley(z, depth):
    G = nx.Graph(); G.add_node(0); frontier = [0]; nxt = 1
    for d in range(depth):
        new = []
        for v in frontier:
            k = z if v == 0 else z - 1
            for _ in range(k):
                G.add_edge(v, nxt); new.append(nxt); nxt += 1
        frontier = new
    return G
for U in (4.0, 6.0, 8.0):
    H, S = build(3, U, 60); e_tree = ground(H)[0]
    for depth in (3, 4):
        G = cayley(4, depth); A = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
        print(f"K=3 U={U} Cayley depth {depth} N={G.number_of_nodes()}: ED {ed_graph(A, U):.6f}  tree-reduced {e_tree:.6f}", flush=True)
for U in (6.0, 8.0):
    H, S = build(3, U, 60); e_tree = ground(H)[0]
    es = [ed_graph(nx.to_scipy_sparse_array(nx.random_regular_graph(4, 80, seed=s), format="csr", dtype=float), U) for s in (1, 2, 3)]
    print(f"K=3 U={U} RRG(4,80): {np.round(es,5)}  tree-reduced {e_tree:.6f}", flush=True)
for K in (3,):
    for U in (30.0, 60.0):
        H, S = build(K, U, 40); e = ground(H)[0]
        print(f"K={K} U={U}: residual beyond 2nd order {e + 3*U + 3*(K+1)/(2*U):.3e}")

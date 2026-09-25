# Notes

## 2026-09-25 — odd fermion clusters / Borromean binding on random networks

### Quantum side

1. **Efimov tower needs 2.3 < d_s < 3.8.** Random graphs (ER, configuration,
   random regular) are tree-like with exponentially decaying return
   probability, i.e. infinite spectral dimension. No Efimov tower there. A
   finite, non-log-periodic Borromean trimer is still possible, as on the 3D
   cubic lattice (Mattis–Rudin 1984): the band edge sets a two-body binding
   threshold U_2 > 0 (as it does whenever d_s > 2), and the question is
   whether the trimer binds first. Computable: two-body on a z-regular tree is
   a one-body problem on the product graph, solved by cavity recursion; the
   three-body needs ED on finite RRGs.
2. **Where an Efimov window could exist:** networks with tunable spectral
   dimension (diamond hierarchical lattices, Song–Havlin–Makse fractal nets,
   Sierpinski-type). Most known ones sit below 2; need constructions landing
   in (2.3, 3.8).
3. **Odd clusters and statistics on graphs.** Quantum statistics of n
   particles on a graph is set by H_1 of the discretized configuration space
   (Ko–Park 2012; Harrison, Keating, Robbins, Sawicki, CMP 2014; Maciążek &
   Sawicki 2019). For 3-connected graphs only bosons/fermions; 2-connected
   allow one anyon phase; graphs with **cut vertices** have H_1 terms that
   grow with n, so three particles carry statistical phases that two do not.
   The giant component of a sparse random graph is full of cut vertices
   (dangling trees), so this n-dependence is generic there. Statistics, not
   binding, but it is the honest "three differ from two" statement for random
   graphs.

### Classical / network-science side

4. **Borromean hyperedges.** A 3-hyperedge with none of its three 2-faces
   present is the combinatorial Borromean object. It is a hypergraph that is
   not a simplicial complex (a 2-simplex would force its edges). In the
   chygraph language this is the gap between a complex and its 1-skeleton's
   clique complex. Real cases: obligate ternary protein complexes (missed by
   pairwise assays such as Y2H), PROTAC ternary complexes with positive
   cooperativity, triadic regulation (Sun, Radicchi, Kurths, Bianconi, Nat.
   Comm. 2023).
5. **Information-theoretic Borromean.** Three variables pairwise independent
   but jointly dependent (parity: Z = X xor Y). Pairwise MI = 0, interaction
   information nonzero, pure synergy in PID. Any pairwise network inference
   (correlation networks, graphical lasso) is blind to it. This is the exact
   statistical analogue of "three bound, two not".
6. **Percolation with three-body activation.** Hyperedge active only if all
   three nodes are (triadic / higher-order percolation, k-core-like
   discontinuous transitions). The "cooperative" regime where no pairwise
   process percolates but the three-body one does is the Borromean regime.
7. **Parity / odd-cycle analogue of "three fermions make a fermion".** Signs
   multiply along cycles: a triangle of three negative edges is frustrated
   (structural balance). Odd cycles are what make a random graph
   non-bipartite; the giant component of G(n,p) acquires odd cycles at the
   same threshold as the giant component, so the "odd cluster" sector is
   never empty above p_c.

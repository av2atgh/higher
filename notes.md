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

## 2026-09-25 — real random structures for odd Borromean fermion clusters

Constraint: two-spin equal-mass fermions (2+1) have no Borromean trimer in
free space; Efimov needs mass ratio > 13.6, Kartavtsev–Malykh trimers > 8.2.
Electron cases need three flavours, heavy-light imbalance, or graph
structure.

1. **Warm dilute nuclear matter** (supernova matter, neutron-star crust
   pasta, heavy-ion fireballs). Mott densities rise with binding: deuteron
   (2.2 MeV) dissolves before t / 3He (~8 MeV) (Röpke; Typel, Röpke, Klähn,
   Blaschke, Wolter, PRC 81, 015803, 2010). Between the two Mott densities the
   three-nucleon clusters are Borromean in the medium. Established.
2. **Amorphous semiconductors / chalcogenide glasses** as continuous random
   networks (Zachariasen–Polk; Weaire–Thorpe tight binding on a random
   4-regular graph). Attraction = Anderson negative-U centres (pairs, no ESR).
   Odd cluster = third carrier on a pair, ESR-active. Borromean regime on
   RRG(4) = the Mattis–Rudin computation in README next steps. Open.
3. **Disordered attractive Fermi systems** (dirty superconductors near
   BCS–BEC crossover, Fermi gases in speckle). Randomness supplies the
   two-body resonance tuning locally and effective mass contrast; Borromean
   trimers as Griffiths rare-region objects, Lifshitz tail below the two-body
   continuum. Own extrapolation; needs three flavours or mass imbalance to be
   clean.

Excluded: 3D Anderson-transition wavefunctions (effective d_s < 2, outside
the Efimov window); trions in disordered TMDs (odd cluster, not Borromean,
exciton always binds).

## 2026-09-25 — where a tree-like mean field pays off

Pick: three-flavour attractive Hubbard on RRG(z) (z = 4 for the continuous
random network) with site disorder. Nuclear matter has no graph (Beth–
Uhlenbeck territory). Two-flavour case has no trion phase.

Already done on the Bethe lattice (do not redo):
* pair sector, clean: DMFT, Keller–Metzner–Schollwöck PRL 2001; Garg,
  Krishnamurthy, Randeria PRL 2005; Toschi, Capone, Castellani PRB 2005.
* pair sector, disordered: cavity, Ioffe & Mézard PRL 2010 (SIT, glassy pair
  phase); Feigel'man, Ioffe, Kravtsov, Cuevas Ann. Phys. 2010 (fractal SC).
* three flavours, clean: DMFT, Inaba & Suga PRA 2009, PRL 2010 (colour SF vs
  trions). Cubic-lattice reference: Rapp, Zaránd, Honerkamp, Hofstetter PRL
  2007 (first-order SF → trion).

Open, and reachable by a three-colour cavity on the tree:
1. colour-SF Tc vs disorder W (three pairing channels in Ioffe–Mézard);
2. SF → trion Fermi-liquid line vs W; expectation: trions localize first, so
   disorder favours the trion side;
3. trion glass at strong W;
4. zero-density anchor: trimer vs pair binding threshold on the tree
   (Mattis–Rudin on RRG); if U_3 < U_2 the trion phase extends to zero
   density.

Cheap first pass: strong-coupling limit → hard-core dimer (t²/U) + trimer
(t³/U²) mixture on a random graph, Bethe–Peierls in the chygraph style;
gives trion-crystal / phase-separation lines at commensurate filling.

Realizations: cold atoms with speckle (clean SU(3)); amorphous chalcogenides
only approximately (three lone-pair p orbitals split by the random
environment).

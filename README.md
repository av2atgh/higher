# higher — when three are needed

Line of thinking started 2026-09-25. Question: Cooper pairs are the key to
BCS; is there any setting where **three** particles are needed to do the job,
and does the idea have a network analogue?

## 1. Physics inventory (2026-09-25)

A literal three-electron condensate is out: three fermions make a fermion,
which cannot Bose-condense. What exists instead:

* **Efimov / Borromean states — three bound, two not.** Efimov (1970): resonant
  two-body interaction gives an infinite tower of trimers with no dimer.
  Realized in cold Cs (Kraemer et al., Nature 2006), the He trimer (Kunitski
  et al., Science 2015), Borromean halo nuclei 6He, 11Li, and three magnons in a
  spin chain (Nishida, Kato, Batista, Nat. Phys. 2013). Requires an effective
  dimension 2.3 < d < 3.8 (Nielsen, Fedorov, Jensen, Garrido, Phys. Rep. 2001).
  On lattices a Borromean trimer can appear below the dimer threshold from band
  structure alone (Mattis & Rudin, PRL 1984; Mattis, RMP 1986).
* **Trions in SU(3) Fermi gases (cold-atom QCD).** Colour superfluid of pairs
  vs Fermi liquid of three-body singlets (Rapp, Zaránd, Honerkamp, Hofstetter,
  PRL 2007). "Cooper triples" on a Fermi sea: Niemann & Hammer PRA 2012;
  Akagami, Inotani, Ohashi PRA 2021; Tajima et al. PRR 2022. The triple wins by
  destroying the superfluid.
* **Z3 Read–Rezayi state.** Pfaffian = p-wave BCS of composite fermions;
  Read–Rezayi k=3 = three-body clustering, candidate for nu = 12/5, Fibonacci
  anyons.
* **Charge-6e superconductivity.** Three Cooper pairs bound into a 6e boson;
  h/6e Little–Parks in CsV3Sb5 rings (Ge et al., PRX 2024); vestigial order of
  a three-component pair density wave.

Accounting: odd fermion clusters give Borromean bound states, trionic Fermi
liquids, or parafermion Hall states. The superconducting analogue needs an
even multiple; 6e is the realized "triple".

## 2. Network analogues (see notes.md)

## Next steps

* DONE 2026-09-25 (trimer/results.md): on the tree of branching K the pair
  threshold is U_2 = 2t(K-1)/sqrt(K) exactly and the three-flavour trimer
  binds at U_3 < U_2 for every K >= 2 (U_3/U_2 from 0.93 at K=2 to 0.84 at
  K=10). Borromean window exists; 2+1 equal-mass fermions have no trimer.
* DONE 2026-09-25 (trimer/results.md, one deep site): a site too shallow
  to bind one particle (|eps0| < eps_c = (K-1)t/sqrt(K)) binds three before
  two at every depth, U_3/U_2 = 0.58-0.70 (K=2,3). Window three to four
  times wider than in the uniform sector. Birman-Schwinger thresholds
  (trimer/bs.py). L=60 runs may still be finishing (scan_bs_K*_L60.log).
* paper/main.tex: after a PRL-style referee report (~/Downloads/prl_review.md)
  the manuscript was revised (all six major points) and retargeted to
  Physical Review E. 2026-09-26: new section "Finite density" (trimer Fermi
  liquid vs colour superfluid, trimer/phase.py, fig3); after a colleague's
  review the manuscript was reframed and then refocused as "Emergent
  hyperedges from pairwise dynamics" (definition + operational test of an
  emergent hyperedge; the Bethe-lattice trimer as the worked case; two
  thresholds on expanders with statistics selecting the sector). Earlier
  versions: paper/main_v1.tex (PRL), main_v2.tex (PRE, physics framing).
* emergent/: percolation of the emergent hyperedges (exact cavity for the
  correlated layer, MC check) and localisation of the composites (book's
  statmech package tables + direct runs); paper Sec. "The emergent
  hyperedges connect and localise", fig4.
* Next: impurity scan at K = 4, 6, 10; finite disorder W on a random
  regular graph (fraction of sites in the Borromean window as a function
  of W and U); three-colour cavity at finite density (notes.md).
* Efimov window on networks: build graphs with spectral dimension in
  (2.3, 3.8) (hierarchical / diamond lattices) and look for the log-periodic
  trimer tower.
* Borromean hyperedges in interaction data: three-way dependence with zero
  pairwise dependence (parity/XOR, interaction information); ternary
  complexes missed by pairwise assays.

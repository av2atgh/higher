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

* Borromean threshold on a random regular graph: does a three-particle bound
  state appear at an attraction U_3 < U_2 (dimer threshold) as on the cubic
  lattice (Mattis–Rudin)? Two-body problem is exact by cavity on the tree;
  three-body by exact diagonalization on RRG(N, z) with N ~ 200.
* Efimov window on networks: build graphs with spectral dimension in
  (2.3, 3.8) (hierarchical / diamond lattices) and look for the log-periodic
  trimer tower.
* Borromean hyperedges in interaction data: three-way dependence with zero
  pairwise dependence (parity/XOR, interaction information); ternary
  complexes missed by pairwise assays.

# Trimer vs pair binding on the Bethe lattice (2026-09-25)

Model: hopping -t on the tree of branching K (degree z = K+1), on-site
attraction -U per pair on the same site, three distinguishable particles
(three flavours). Code: tree3.py (reduced model), checks*.py (validation),
scan.py / scan2.py (results), logs scan.log, scan2.log, scan.json.

## Sector

Automorphism-invariant sector of the infinite tree: wavefunction depends
only on the orbit of the configuration. For three points the orbit is the
triple of depths (a,b,c) below the median, positive-depth particles in
distinct branches; orbit multiplicity n = z^(k) (z-1)^(S-k). Relative motion
is square-summable with that measure, centre of mass uniform. This is the
tree analogue of the zero-total-momentum sector of a lattice, and it is what
a random regular graph gives in the thermodynamic limit once rare short
loops are ignored. Continuum edge per free particle: -2 sqrt(K) t, the band
edge of the localization chapter, Eq. (cleantree) in the book.

Caveat: for a tightly bound composite, this sector puts the centre of mass
in the Perron (uniform) mode, energy E_int - z t_eff, which lies below the
L2 band bottom E_int - 2 sqrt(K) t_eff of the composite. A finite Cayley
tree confines the centre of mass and sits above; a finite random regular
graph adds short loops and sits below. Both were checked (checks2.py).

## Pair: exact

Relative coordinate d = distance. Symmetrised: half line with hopping
J = 2t sqrt(K) for d >= 1, J0 = 2t sqrt(K+1) between d = 0 and 1, site
energy -U at d = 0. Cavity element of the bulk half line g = 1/J at the
edge, and G_00^{-1} = -U - E - J0^2 g = 0 gives

    U_2 = 2J - J0^2/J = 2 t (K-1)/sqrt(K),      E_2(U) = -J(l + 1/l),
    U = J(l + 1/l) - J0^2 l/J,   0 < l < 1.

K = 1 (line) gives U_2 = 0, as it must. Large U: E_2 = -U - 4 z t^2/U,
i.e. on-site shift -2z t^2/U plus the Perron energy of the pair hopping
2t^2/U on the tree.

## Trimer: numerical, invariant sector, truncation a+b+c <= L

Threshold U_3 = smallest U with E_3 < -6 sqrt(K) t (bisection, tol 2e-3;
the truncated continuum lies above the edge so this is an upper bound that
was already converged at L = 60 and unchanged at L = 120, 180).

| K  | z  | U_3    | U_2    | U_3/U_2 | window (U_2-U_3)/U_2 |
|----|----|--------|--------|---------|----------------------|
| 2  | 3  | 1.310  | 1.4142 | 0.926   | 7.4 %                |
| 3  | 4  | 2.073  | 2.3094 | 0.898   | 10.2 %               |
| 4  | 5  | 2.637  | 3.0000 | 0.879   | 12.1 %               |
| 6  | 7  | 3.508  | 4.0825 | 0.859   | 14.1 %               |
| 10 | 11 | 4.800  | 5.6921 | 0.843   | 15.7 %               |

**The trimer binds before the pair on every tree with K >= 2: a Borromean
window U_3 < U < U_2, widening with branching.** Above U_2 the trimer stays
below the dimer + free particle continuum E_2(U) - 2 sqrt(K) t at every U
scanned (scan.log, column "bound?"), so the trion is the ground composite
at all U > U_3 in the three-flavour model, on the tree of degree 4 as on the
cubic lattice (Mattis and Rudin, PRL 1984).

## 2+1 sector (two identical fermions + one distinguishable, equal mass)

Antisymmetric under exchange of the two identical particles, attraction only
between unlike particles. No trimer: E_3 stays above the dimer + particle
continuum for all U up to 14 (K = 2) and 6 (K = 3), with the gap shrinking
slowly (scan2.log). Consistent with free space, where equal-mass 2+1 needs
a mass ratio above 8.2 (Kartavtsev and Malykh) for a trimer. So the
Borromean effect on the tree requires three flavours, as stated in notes.md.

## Validation (checks.py, checks2.py, checks3.py)

* K = 1: reduced model vs exact diagonalisation of three particles on a
  41-site ring, agreement to 1e-8 at U = 1 and 3.
* Pair: numerical half line vs closed form, 1e-8.
* Large U: E_3 = -3U - 3z t^2/(2U) + O(t^4/U^3); residual at U = 30, 60
  scales as 1/U^3.
* Move rules: total out-weight 3z on every orbit and detailed balance
  n(s) M(s,s') = n(s') M(s',s) on every pair, K = 1, 2, 3, 5.
* Finite graphs bracket the reduced model as the sector argument predicts:
  Cayley tree of degree 4 (depth 3 and 4) above, random 4-regular graphs
  (N = 40, 80) below, both by a few 1e-2 at U = 4-8.

## Consequences for the phase diagram (notes.md, 2026-09-25 entry)

At zero density on the random regular graph the trion, not the pair, is
the first composite to form. In the three-flavour attractive Hubbard model
on the tree the trion Fermi liquid therefore extends down to zero density
for U > U_3 and the colour superfluid needs a finite density to compete.
The Borromean window (7 to 16 %) is the natural place to look for a
disorder-driven transition: site disorder that pushes local coupling
across U_3 but not U_2 creates trions and no pairs.

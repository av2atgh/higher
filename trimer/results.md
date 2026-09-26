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
continuum for all U up to 14 (K = 2 and K = 3), with the gap shrinking
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

## One deep site (dilute disorder), 2026-09-25

Code tree4.py (orbits under the root stabiliser, impurity eps0 at the root,
validated to 1e-5 against Cayley-tree ED at m=2,3), bs.py (Birman–Schwinger
threshold at the continuum edge, U_c = 1/lambda_max of V^{1/2}(H0-E)^{-1}V^{1/2};
reproduces bisection to its 2e-3 grid), scan_bs.py, logs scan_bs_K{K}_L{L}.log.

One particle binds to the site for |eps0| > eps_c = (K-1)t/sqrt(K) (exact,
cavity resolvent at the edge). For 0 < |eps0| < eps_c (L=40 trimer, 80 pair):

| |eps0|/eps_c | K=2 U_2 | K=2 U_3 | ratio | K=3 U_2 | K=3 U_3 | ratio |
|---|---|---|---|---|---|---|
| 0    | 1.985 | 1.391 | 0.70 | 3.264 | 2.172 | 0.67 |
| 0.3  | 1.954 | 1.285 | 0.66 | 3.197 | 1.961 | 0.61 |
| 0.5  | 1.790 | 1.118 | 0.62 | 2.840 | 1.665 | 0.59 |
| 0.8  | 1.261 | 0.759 | 0.60 | 1.908 | 1.099 | 0.58 |
| 0.9  | 0.979 | 0.598 | 0.61 | 1.453 | 0.849 | 0.58 |
| 0.99 | 0.596 | 0.415 | 0.70 | 0.819 | 0.545 | 0.67 |

U_3 < U_2 at every depth; the window is 30-42 %, three to four times the
uniform-sector window. eps0 -> 0 gives the square-summable thresholds of the
clean tree (L=60: K=2 1.974/1.379, K=3 3.257/2.165; still drifting down by
~0.01 from L=40), above the uniform-sector values 1.414/1.310 and
2.309/2.073, as the Perron-vs-band-bottom argument predicts.

## 2+1 sector at strong coupling (analytic), 2026-09-26

Dimer D = (up, down) on site i, energy -U; extra up fermion on j != i.
First order in t: the down of the dimer hops i -> j, which is a degenerate
state (dimer at j, free up at i). Operator algebra:
  -t c†_{j,down} c_{i,down} c†_{i,up} c†_{i,down} c†_{j,up} |0>
    = +t c†_{j,up} c†_{j,down} c†_{i,up} |0>,
i.e. the exchange amplitude is +t, opposite in sign to a hop. In the
relative coordinate d = d(i,j) >= 1 (invariant sector) the exchange is a
diagonal term +t at d = 1; the up cannot hop onto i (Pauli), so d = 0 is
absent. The relative Hamiltonian at order t^0 (t/U) is the tree adjacency
with the dimer site removed plus a positive potential +t on its neighbours,
whose spectrum lies above -2 sqrt(K) t (subgraph of the tree plus a positive
operator). The dimer itself hops only at order 2t^2/U. Hence no 2+1 bound
state below the dimer + fermion continuum at strong coupling on the tree;
the exchange, which binds bosonic trimers in 1D (Valiente, Petrosyan,
Saenz 2010), is repulsive for the fermionic 2+1 case. Consistent with the
finite-L scan (gap to the continuum positive, decreasing as the hard-wall
energy of a free relative motion) and with Mattis-Rudin's fermion result and
Abdullaev, Khalkhuzhaev, Kholmatov (arXiv:2502.01099): equal masses on the
3D lattice, no discrete spectrum below the continuum at large coupling;
trimers require a mass ratio above a threshold.

## Revision after the PRL report (2026-09-26)

Square-summable (L2) sector on the clean tree (tree4 with eps0 = 0):
* pair, exact from the radial kernel (l2_pair.py): 1/U_2 = sum_d b_d P_d(2 sqrt K),
  b_d = <ii|(H0 + 4 sqrt K)^{-1}|jj>; uniform sector = sum_d b_d n_d, reproduces
  2(K-1)/sqrt K after 1/D Richardson. Values: K=2 1.965677, 3 3.251185, 4 4.268903,
  5 5.136053, 6 5.903985, 10 8.402470.
* trimer, Birman-Schwinger on L = 30..60 extrapolated with fitted power (1.7-1.9;
  same fit on the pair reproduces the exact values to 1e-3) (l2_trimer.py, extrap.py):
  K=2 1.369, 3 2.159, 4 2.741, 5 3.218, 6 3.633, 10 4.949; ratios 0.696, 0.664,
  0.642, 0.627, 0.615, 0.589; windows 30-41 %.
* uniform K=5: U3 = 3.103, U2 = 3.578, ratio 0.867.

Cubic lattice, zero total momentum, box in relative coordinates (cubic.py):
pair 8.06 at R=14 -> 7.92 by 1/R extrapolation (Watson 7.9136); trimer 5.186,
5.168, 5.163, 5.160, 5.159 at R = 2..6 -> 5.158. Window 35 %, U3/U2 = 0.652.

Efimov count (efimov_count.py): uniform sector at U = U2^u: one state below the
edge for K = 2, 3 at L = 60, 120, 180; next state at +0.028, +0.0097, +0.0048
(K=2), i.e. L^-2, continuum. L2 sector: 7-8 states below the edge at U = U2 are the
box-quantised centre-of-mass band of the single bound trimer; not a count.

Finite RRG, Q-sector (orthogonal to Perron on every particle), K=2 (rrg_q2.py,
rrg_girth.py, rrg_pair_scaling.py): see logs; pair coincidence does not vanish
between N = 100 and 200 even at girth >= 6; trimer coincidence ~0.2 at U = 1.7,
300x its U = 0 value; pair scaling to N = 1600 in rrg_pair_scaling.log.

Disorder average, isolated-site approximation (disorder_avg.py, L=60 curves):
W = 2 eps_c: U = 0.5 U2 -> Borromean 0.14/0.17 (K=2/3), pair 0.05/0.07;
U = 0.7 U2 -> 0.37/0.35 and 0.13/0.15. W = 4 eps_c: 0.25 hold one particle,
Borromean 0.07/0.09 at U = 0.5 U2. Finite disordered Cayley ED (disorder_ed.py)
abandoned: coincidence probabilities measure the IPR of the lowest orbital.

Finite density (bcs.py, l2_energies.py, phase.py; 2026-09-26): BCS mean field of
two colours on the KM DOS has 2 mu -> E_2^u as n -> 0 (the condensate is the
Perron mode); a Fermi liquid of trimers fills the L2 band with hopping
t_3 = (E3_L2 - E3_u)/(K+1-2 sqrt K). Hartree -3U n^2 common, dropped. Mixed
phase allowed. Trimer liquid wins at n -> 0 for all U > U_3 (L2); superfluid
enters at n_1(U), see phase_K{2,3}.json and paper fig3. Runs at 1.4 GB peak
each; never run two in parallel with an RRG N=200 job (16 GB machine).
Onset of the superfluid (phase_onset_K{2,3}.json; trimer Fermi level = binding):
K=2: n1 = 0.003 (U=1.38), 0.015 (1.40), 0.030 (1.42), 0.052 (1.45), 0.090 (1.50),
0.128 (1.55), 0.167 (1.60), 0.245 (1.70); crosses 0.1 at 1.51. K=3: 0.004 (2.17),
0.030 (2.20), 0.085 (2.25), 0.144 (2.30), 0.207 (2.35), 0.271 (2.40); crosses 0.1
at 2.26. Above n1 the superfluid fraction stays < 3% at n = 0.1 and < 5% at 0.3
(trimer band narrow: t3 = 0.17 / 0.14 near threshold). Paper Sec. "Finite density", fig3.

Colleague's review (2026-09-26) applied: manuscript reframed around the one new
physics point (two thresholds on a non-amenable graph, Alon-Boppana gap; statistics
selects the sector: bosonic composites condense in the Perron mode -> uniform
threshold, fermionic fill a band -> L2 threshold). Cross-sector ratio U3(L2)/U2^u
= 0.968, 0.935, 0.914, 0.899, 0.890, 0.869 (K = 2..10), window 3-13 %: the one a
trion liquid enjoys. Large-K pair limits (semicircle_limit.py): U2^u -> 2 t*,
U2 -> 3.3075 t* (1/b0, b0 = 0.302347), ratio 1.654; a + b/sqrt K fits of the
trimer ratios -> 0.50 (L2), 0.78 (u), 0.79 (cross), indicative only (same fit on
the pair gives 1.54 for 1.654). Known-physics recapitulated: Mattis-Rudin
mechanism, Mattis 2+1 no-go (exchange argument lattice-independent by Cauchy
interlacing), no Efimov tower, trionic phase (Rapp). 6Li in optical lattices
tests the cubic window, not the tree; hyperbolic circuit-QED lattices named as
the home of the two-sector physics.

Refocus (2026-09-26): manuscript now "Emergent hyperedges from pairwise
dynamics". Sec. II defines an emergent hyperedge (irreducible k-body bound state;
without faces = Borromean) and the operational test on finite graphs (P_k finite
as N grows, P_{k-1} -> 0; free: N^{1-k}; bound (k-1)-cluster + free: 1/N). The
RRG girth-6 section is the test. Nine network-science/information references
added from memory (battiston2021, bick2023, bianconi2021, lambiotte2019,
iacopini2019, neuhauser2020, sun2023, rosas2019, williams2010): verify.

Emergent hyperedges: percolation and localisation (emergent/, 2026-09-26).
* perc_cavity.py: exact cavity for correlated site percolation on the tree (site active if
  within r of a Borromean centre; message = distance to nearest centre in the subtree +
  connection probability conditioned on the outside distance). p_c(r): K=2: 0.499 (r=0),
  0.0729, 0.0137, 0.0030; K=3: 0.333, 0.0239, 0.0022, 0.00023. MC on RRG N=1e5 agrees
  (perc_mc.py). Heterogeneous radii version checked against single-radius.
* trimer_density.py: r90 of site-bound trimers 2..10 (site_density_K*.json); s3, s2
  (uniform_density.json). perc_map.py: S = 1.000 wherever f_B > 0, also compact-only.
* trimer_localization.py + tail_fraction.py + mobility_ext.py + glass_line.py (statmech
  package: anderson_wc/lines tables, Ensemble.mobility_edge at W' = 14,16 (3,0) and
  16..28 (4,0); regular_instance ED for the tail fraction): W_c^(3) = 2.1 -> 0.55 t (K=2,
  U = 1.38 -> 4), 2.6 -> 0.65 t (K=3); W_c^(2) ~ 4 t (K=2), 6 t (K=3) once bound;
  W_glass(U): K=2 0.6 (1.38) -> 1.4 (1.5) -> 1.7 (1.7); K=3 1.1 -> 2.1.
  Direct mobility-edge check at W'=2.79 agrees with the table interpolation to 2%.

Review round 2 (2026-09-26): t3 checks: strong coupling 3t^3/(2U^2); trimer band width in the
root-fixed sector (emergent/trimer_bandwidth.py; band = box ladder before the gap, the isolated
L-independent state above it is an internal excitation): t3 = 0.089/0.149/0.201 (K=2, U=4/3/2.5)
vs sector 0.087/0.135/0.166; K=3: 0.058/0.089/0.114 (U=5/4/3.5) vs 0.057/0.083/0.102. Composite-
impurity route (composite_impurity.py) abandoned: box artefact. Localisation quoted for compact
trimers only (U >= 1.7 K=2, >= 2.5 K=3), 20% uncertainty near the low end. Extra points U = 1.8,
1.9 (K=2), 2.6 (K=3). Single-particle glass line at equal density (glass_line_single.json):
15-16 t (K=2), 29-31 t (K=3). Alon-Boppana -> Friedman keeps the gap open, AB says the tree
saturates it. Percolation section cut to a paragraph in Sec. X.

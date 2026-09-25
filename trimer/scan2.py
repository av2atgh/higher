import numpy as np
from tree3 import pair_threshold, pair_energy
from scan import E3, bisect_threshold
print("branching trend (L=60, tol 2e-3):")
for K in (2, 3, 4, 6, 10):
    U2 = pair_threshold(K); U3 = bisect_threshold(K, 60, tol=1e-3)
    print(f"  K={K:2d} z={K+1:2d}: U3={U3:.4f}  U2={U2:.4f}  U3/U2={U3/U2:.4f}  (U2-U3)/U2={(U2-U3)/U2:.4f}", flush=True)
print("2+1 fermion sector at larger U (L=90): E3_anti vs dimer+particle continuum")
for K in (2, 3):
    for U in (6.0, 8.0, 10.0, 14.0):
        ea = E3(K, U, 90, sector="anti21"); e2 = pair_energy(K, U) - 2*np.sqrt(K)
        print(f"  K={K} U={U:5.1f}: E3_anti={ea:10.5f}  E2+edge1={e2:10.5f}  diff={ea-e2:+.5f}", flush=True)

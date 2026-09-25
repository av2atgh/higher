import numpy as np, json
from scan import E3
from tree3 import pair_energy
K = 3
rows = []
for U in np.arange(1.9, 2.81, 0.05):
    rows.append([U, E3(K, U, 90), pair_energy(K, U)])
json.dump(rows, open("fine_K3.json", "w"))
print(np.array(rows))

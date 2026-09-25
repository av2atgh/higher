import numpy as np
from tree3 import states, moves
for K in (1, 2, 3, 5):
    z = K + 1
    S = states(30, z); idx = {s: i for i, s in enumerate(S)}
    def n(s):
        k = sum(1 for d in s if d > 0); ff = 1
        for j in range(k): ff *= (z - j)
        return ff * (z - 1) ** (sum(s) - k)
    bad = 0; badsum = 0
    for s in S:
        tot = 0
        fwd = {}
        for s2, w in moves(s, z):
            tot += w
            fwd[s2] = fwd.get(s2, 0) + w
        if tot != 3 * z: badsum += 1
        for s2, w in fwd.items():
            if sum(s2) > 30: continue
            back = sum(w2 for s3, w2 in moves(s2, z) if s3 == s)
            if n(s) * w != n(s2) * back: bad += 1
    print(f"K={K}: states {len(S)}, out-weight != 3z: {badsum}, detailed-balance violations: {bad}")

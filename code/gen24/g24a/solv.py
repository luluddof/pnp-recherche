import johnson
bad = []
for k in range(1, 26):
    for j1 in range(0, k):
        j2 = k - 1 - j1
        if johnson.solve(k, j1, j2) is None: bad.append((k, j1, j2))
print("unsolvable (k, j1, j2) with j1+j2 = k-1, k <= 25:", bad)

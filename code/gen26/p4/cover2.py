import sys, glob, cover
N, p, q = sys.argv[1:4]
for r, n in [(3, 46), (4, 73), (5, 36)]:
    fs = glob.glob(f"cubechain_{N}_{p}_{q}_r{r}*.txt")
    print(f"r={r}:", cover.cov(fs, n, skip=(0,)))

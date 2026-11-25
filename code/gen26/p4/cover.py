import re, sys
def cov(files, n, skip=()):
    st = {}
    for f in files:
        try:
            for line in open(f):
                m = re.match(r"\s+cube (\d+) .*: (\w[\w-]*) ", line)
                if m:
                    i, s = int(m.group(1)), m.group(2)
                    if st.get(i) != 'UNSAT': st[i] = s
        except FileNotFoundError: pass
    uns = sorted(i for i in st if st[i] == 'UNSAT'); other = {i: st[i] for i in st if st[i] != 'UNSAT'}
    missing = [i for i in range(n) if i not in st and i not in skip]
    return len(uns), other, missing
#print("r=4:", cov(["cube_593_r4_a.txt", "cubechain_5_9_3_r4.txt", "cubechain_5_9_3_r4_a.txt", "cubechain_5_9_3_r4_b.txt"], 73, skip=(0,)))
#print("r=5:", cov(["cubechain_5_9_3_r5.txt"], 36, skip=(0,)))

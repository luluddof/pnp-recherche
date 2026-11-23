import subprocess, glob, os, sys, json
from concurrent.futures import ThreadPoolExecutor
H = os.path.dirname(os.path.abspath(__file__))
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 2000000
pat = sys.argv[2] if len(sys.argv) > 2 else 'r*'
heurs = sys.argv[3].split(',') if len(sys.argv) > 3 else ['var','pair','clause','linla','varla']
jobs = [(h, f) for f in sorted(glob.glob(f'{H}/inst/{pat}.cnf')) for h in heurs]
def run(j):
    h, f = j
    out = subprocess.run([f'{H}/xdpll', h, str(LIMIT)], stdin=open(f), capture_output=True, text=True, timeout=3000).stdout.split()
    return h, os.path.basename(f), out[0], int(out[1]), int(out[2]), float(out[3])
res = list(ThreadPoolExecutor(4).map(run, jobs))
with open(f'{H}/res_{pat.replace("*","")}.json', 'w') as fo: json.dump(res, fo)
bad = [r for r in res if r[2] not in ('UNSAT', 'LIMIT')]
print('non-UNSAT:', bad)

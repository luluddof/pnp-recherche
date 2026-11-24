#!/bin/bash
# usage: runlist.sh listfile R maxiter ; each line "N p q"
cd /tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen24/g24c
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cat $1 | xargs -P ${P:-3} -L 1 bash -c 'timeout 580 python3 vp.py $0 $1 $2 '"$2"' 7 '"$3"' > logs/vp_$0_$1_$2.txt 2>&1'

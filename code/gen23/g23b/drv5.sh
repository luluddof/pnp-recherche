for m in 0 1 2 3 4; do timeout 300 python3 mysat.py 5 $((4-m)) $((7-m)) $m; done
for m in 0 1 2 3 4 5 6; do timeout 300 python3 mysat.py 5 $((6-m)) $((6-m)) $m; done
timeout 300 python3 mysat.py 5 7 7 0

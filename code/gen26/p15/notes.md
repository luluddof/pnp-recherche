# G26-P15 independent verification of G26-P11 exact points
Script: myverify.py (pure-Python Fraction Gaussian elimination + integer Bareiss, no numpy/flint, verify_exact.py not used).
Convention (from JSON): C[(x,y),(u,v)] = sum_l E_l[x,u] C_l[y,v]; X[(x,y),(u,v)] = [x=v][y=u]. All entries integers.
Sanity (sanity.py): C=0 gives rank 36 at N=6; perturbing one entry of the (8,6)@6 solution gives rank 10.

| point | rank_Q(X-C) | #terms | osr(C) | verdict | Pareto vs known |
|---|---|---|---|---|---|
| (8,6)@6 | 8 | 6 | 6 | CONFIRMED | yes (known: 15@q<=2, 12@3, 7@7) |
| (20,3)@7 | 20 | 3 | 3 | CONFIRMED | yes (known: 21@q<=2 only) |
| (17,4)@7 | 17 | 4 | 4 | CONFIRMED | yes |
| (14,5)@7 | 14 | 5 | 5 | CONFIRMED | yes |
| (12,6)@7 | 12 | 6 | 6 | CONFIRMED | yes |
| (11,7)@7 | 11 | 7 | 7 | CONFIRMED | yes |
| (22,5)@8 | 22 | 5 | 5 | CONFIRMED | yes (known: 28@q<=2, 8@8) |
| (20,6)@8 | 20 | 6 | 6 | CONFIRMED | yes |
| (22,6)@9 | 22 | 6 | 6 | CONFIRMED | yes (known: 36@q<=2) |

Also beats padding bound f(N+1,q) <= f(N,q)+2N+1 from known/claimed smaller-N points
(e.g. (7,7)@6 -> (20,7)@7; (8,8)@8 -> (25,8)@9).
p-1 impossibility: not attempted (no cheap argument found; trivial bounds are vacuous).

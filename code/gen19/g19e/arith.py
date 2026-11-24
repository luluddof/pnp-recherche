# T8: the arithmetic corollary, brute force: integers r >= 0, 1 <= b < M^2 with (a),(b) for all k <= M and (c) r >= Q  ==> 9 r b >= 2 Q M^2.
# Also WITHOUT (c) (border case): the corollary needs b <= M^2/2 (as stated in G18f) -- show the failure otherwise.
worst = 9; fail_border = None; n = 0
for M in range(1, 13):
    for Q in (1, 2, 3, 4, 7, 16, 30):
        for b in range(1, M*M):
            # least r allowed by (a),(b)
            r = max(k*Q*(M - k*b) for k in range(M+1)); r = max(r, 0)
            while any(b*Q < k*(M*Q - k*r) for k in range(M+1)): r += 1
            if 9*r*b < 2*Q*M*M and fail_border is None and 2*b > M*M: fail_border = (M, Q, r, b)
            if 2*b <= M*M: assert 9*r*b >= 2*Q*M*M, (M, Q, r, b)      # G18f statement (border-valid)
            r2 = max(r, Q); assert 9*r2*b >= 2*Q*M*M, (M, Q, r2, b); n += 1
            worst = min(worst, 9*r2*b/(2*Q*M*M))
print(f"T8 ok on {n} (M,Q,b): with (c), 9rb >= 2QM^2 for all 1 <= b < M^2; worst ratio rb/(2QM^2/9) = {worst:.4f}")
print("   without (c) and with b > M^2/2 the inequality system alone allows e.g. (M,Q,r,b) =", fail_border)

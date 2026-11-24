"""Accounting floors for D = 2 (max_x L(x)+L(M-x))/(4 lambda), lambda = 1, h1 = 1; and the D-direction table."""
def rho(L, M, steps=4000):
    return max(L(M*i/steps) + L(M*(steps-i)/steps) for i in range(steps+1))/4
print("partitions-like sqrt profile, (V): M=4h1      ", rho(lambda E: E**.5, 4))
print("linear profile, no (V): E0=h1, H=0, M=2h1     ", rho(lambda E: E, 2))
print("linear profile with (V): M=4h1 (cap 2 lambda) ", rho(lambda E: min(2, E), 4))
print("staircase with (V): M=4h1                     ", rho(lambda E: 0 if E < 1 else (1 if E <= 3 else 2), 4))
print("G9a concave extremal profile, (V)             ", rho(lambda E: min(E, 2/3 + E/3), 4))
print("D  1/D  D^-1/2  D^-(D-1)/D")
for D in (2, 3, 4, 6, 8): print(D, round(1/D, 4), round(D**-.5, 4), round(D**(-(D-1)/D), 4))

# G22b — anatomy of the exact small solutions; theory of the monomial shape   (working notes, saved as I go)
Tags: [P] proved here (short proof, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [L] literature (from memory).
Convention (G19b): X(x,y,u,v) = [x=v][y=u], C = sum_l G_l[x,u] G'_l[y,v], R = X - C, rr(R) <= p (flattening (x,y)|(u,v)).
Operator form: C(W) = sum_l G_l W G'_l^T, X(W) = W^T.

## 0. First findings (to be detailed below)
 * [E] In the four "symmetric" exact solutions ((5,3)@4 twice, (4,4)@4, (7,7)@6): G'_l = G_l^T, the supports M_l are disjoint,
   and R is a DIRECT SUM OF p RANK-ONE +-1 BLOCKS: R(x,y,u,v) = sum_j H_j[x,y] H_j[v,u], H_j signed partial permutations with
   disjoint supports covering exactly the complement T of S = union M_l.  So
        W^T = sum_l G_l W G_l + sum_j tr(H_j W) H_j      (all W),
   and the N x N array "cell (a,b) -> its class (G_l or H_j) and sign" is a CONSISTENTLY SIGNED INTERCALATE MATRIX
   (= integer sum-of-squares formula [N,N,p+q]) whose intercalations always pair one G-colour with one H-colour.

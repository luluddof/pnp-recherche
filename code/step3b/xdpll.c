// DPLL over affine subspaces of GF(2)^n (tree-like Res(+) style), n <= 64.
// State = RREF linear system (fully reduced, rows indexed by pivot var).
// Propagation modes: clause is evaluated by projecting the affine subspace onto its
// undetermined variables (exact for <=3 vars): forbidden point outside -> satisfied;
// projection = 2 points incl. forbidden -> force other point; 1 point -> conflict.
// Heuristics: var, pair, clause, linla (lookahead over var/pair/clause-parity forms),
//             varla (lookahead over variables only; fair baseline for linla).
// Usage: xdpll HEUR [leaflimit] < dimacs   -> prints "UNSAT leaves nodes secs" | "SAT ..." | "LIMIT ..."
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
typedef uint64_t u64;
#define MAXM 4096
#define K 3
int n, m;
int clen[MAXM], cvar[MAXM][K], cneg[MAXM][K];
typedef struct { u64 piv; u64 row[64]; u64 rhs; } St;  // rhs bit v = rhs of row with pivot v
long long npar = 0, nbr = 0, leaves = 0, nodes = 0, limit = 2000000000LL;
int heur = 0;
// per-propagation outputs (active clause info)
int nact, act[MAXM], actk[MAXM]; u64 actU[MAXM]; u64 actF[MAXM]; // U mask, forbidden value bits on U

static inline u64 reduce(const St *s, u64 x, int *r) {
  u64 t = x & s->piv; int rr = 0;
  while (t) { int v = __builtin_ctzll(t); x ^= s->row[v]; rr ^= (s->rhs >> v) & 1; t &= t - 1; }
  *r = rr; return x;
}
static inline int inspan(const St *s, u64 x) { int r; return reduce(s, x, &r) == 0; }
// returns -1 conflict, 0 redundant, 1 added
static int add_eq(St *s, u64 x, int r) {
  int rr; x = reduce(s, x, &rr); r ^= rr;
  if (!x) return r ? -1 : 0;
  int p = __builtin_ctzll(x);
  u64 t = s->piv;
  while (t) { int v = __builtin_ctzll(t); t &= t - 1;
    if ((s->row[v] >> p) & 1) { s->row[v] ^= x; if (r) s->rhs ^= 1ULL << v; } }
  s->piv |= 1ULL << p; s->row[p] = x; if (r) s->rhs |= 1ULL << p; else s->rhs &= ~(1ULL << p);
  return 1;
}
static inline int det(const St *s, int v, int *val) {
  if (((s->piv >> v) & 1) && s->row[v] == (1ULL << v)) { *val = (s->rhs >> v) & 1; return 1; }
  return 0;
}
// propagate to fixpoint; returns 0 conflict, 1 ok. fills act arrays.
static int propagate(St *s) {
  for (;;) {
    int changed = 0; nact = 0;
    for (int c = 0; c < m; c++) {
      int uv[K], uf[K], ku = 0, sat = 0;
      for (int i = 0; i < clen[c]; i++) { int val, v = cvar[c][i];
        if (det(s, v, &val)) { if (val != cneg[c][i]) { sat = 1; break; } }
        else { uv[ku] = v; uf[ku] = cneg[c][i]; ku++; } }  // forbidden value of var = cneg (lit false)
      if (sat) continue;
      if (ku == 0) return 0;
      if (ku == 1) { int r = add_eq(s, 1ULL << uv[0], !uf[0]); if (r < 0) return 0; changed = 1; continue; }
      // projection onto U
      u64 forms[8]; int fv[8], nf = 0;
      for (int sub = 1; sub < (1 << ku); sub++) {
        u64 x = 0; for (int i = 0; i < ku; i++) if ((sub >> i) & 1) x |= 1ULL << uv[i];
        int r; if (reduce(s, x, &r) == 0) { forms[nf] = sub; fv[nf] = r; nf++; } }
      int cnt = 0, fok = 0, other = -1, fpt = 0;
      for (int i = 0; i < ku; i++) if (uf[i]) fpt |= 1 << i;
      for (int a = 0; a < (1 << ku); a++) { int ok = 1;
        for (int j = 0; j < nf; j++) if ((__builtin_popcount(a & (int)forms[j]) & 1) != fv[j]) { ok = 0; break; }
        if (ok) { cnt++; if (a == fpt) fok = 1; else other = a; } }
      if (!fok) continue;               // satisfied on whole subspace
      if (cnt == 1) return 0;           // cannot happen (U undetermined) but safe
      if (cnt == 2) {                   // force other point
        for (int i = 0; i < ku; i++) { int r = add_eq(s, 1ULL << uv[i], (other >> i) & 1); if (r < 0) return 0; }
        changed = 1; continue; }
      u64 U = 0, F = 0; for (int i = 0; i < ku; i++) { U |= 1ULL << uv[i]; if (uf[i]) F |= 1ULL << uv[i]; }
      act[nact] = c; actk[nact] = ku; actU[nact] = U; actF[nact] = F; nact++;
    }
    if (!changed) return 1;
  }
}
static int rank_of(const St *s) { return __builtin_popcountll(s->piv); }

// ---------- heuristics: produce candidate forms ----------
static void var_scores(const St *s, double *sc) {
  double p[64] = {0}, q[64] = {0};
  for (int i = 0; i < nact; i++) { double w = actk[i] == 2 ? 5 : 1; u64 U = actU[i];
    while (U) { int v = __builtin_ctzll(U); U &= U - 1; if ((actF[i] >> v) & 1) q[v] += w; else p[v] += w; } }
  for (int v = 0; v < n; v++) sc[v] = p[v] * q[v] * 1024 + p[v] + q[v];
}
static u64 best_var(const St *s) {
  double sc[64]; var_scores(s, sc); int b = -1;
  for (int v = 0; v < n; v++) { int val; if (det(s, v, &val)) continue; if (b < 0 || sc[v] > sc[b]) b = v; }
  return 1ULL << b;
}
static int pc[64][64];
static u64 best_pair(const St *s) {
  memset(pc, 0, sizeof pc);
  for (int i = 0; i < nact; i++) { int w = actk[i] == 2 ? 5 : 1; int vs[K], k = 0; u64 U = actU[i];
    while (U) { vs[k++] = __builtin_ctzll(U); U &= U - 1; }
    for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) pc[vs[a]][vs[b]] += w; }
  int bi = -1, bj = -1, bs = -1;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (pc[i][j] > bs) {
    u64 x = (1ULL << i) | (1ULL << j); if (inspan(s, x)) continue; bs = pc[i][j]; bi = i; bj = j; }
  if (bi < 0 || bs <= 0) return best_var(s);
  return (1ULL << bi) | (1ULL << bj);
}
// form to branch on for active clause i: parity of U if unknown, else a pair inside U (closes the
// "wrong parity" side quickly). *dim = dimension of projection of subspace onto U.
static u64 clause_form(const St *s, int i, int *dim) {
  u64 U = actU[i]; int vs[K], k = 0; u64 T = U; while (T) { vs[k++] = __builtin_ctzll(T); T &= T - 1; }
  int rk = 0; // rank of span ∩ <e_U>
  for (int sub = 1; sub < (1 << k); sub++) { u64 x = 0; for (int j = 0; j < k; j++) if ((sub >> j) & 1) x |= 1ULL << vs[j];
    if (inspan(s, x)) rk++; }
  int r = 0; while ((1 << r) - 1 < rk) r++;  // rk = 2^r - 1 nonzero forms
  *dim = k - r;
  if (!inspan(s, U)) return U;
  for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) { u64 x = (1ULL << vs[a]) | (1ULL << vs[b]); if (!inspan(s, x)) return x; }
  return 0;
}
static u64 best_clause(const St *s) {
  double sc[64]; var_scores(s, sc); u64 bx = 0; double bsc = -1; int bd = 99;
  for (int i = 0; i < nact; i++) { int d; u64 x = clause_form(s, i, &d); if (!x) continue;
    double t = 0; u64 U = actU[i]; while (U) { t += sc[__builtin_ctzll(U)]; U &= U - 1; }
    if (d < bd || (d == bd && t > bsc)) { bd = d; bsc = t; bx = x; } }
  if (!bx) return best_var(s);
  return bx;
}
// lookahead over candidate list
#define TOPK 15
static u64 lookahead(const St *s, int linear) {
  u64 cand[200]; int nc = 0;
  double sc[64]; var_scores(s, sc);
  // top vars
  int used[64] = {0}; int nv = linear ? TOPK : 3 * TOPK;
  for (int t = 0; t < nv; t++) { int b = -1; for (int v = 0; v < n; v++) { int val;
      if (used[v] || det(s, v, &val)) continue; if (b < 0 || sc[v] > sc[b]) b = v; }
    if (b < 0) break; used[b] = 1; cand[nc++] = 1ULL << b; }
  if (linear) {
    // top pairs
    memset(pc, 0, sizeof pc);
    for (int i = 0; i < nact; i++) { int w = actk[i] == 2 ? 5 : 1; int vs[K], k = 0; u64 U = actU[i];
      while (U) { vs[k++] = __builtin_ctzll(U); U &= U - 1; }
      for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) pc[vs[a]][vs[b]] += w; }
    for (int t = 0; t < TOPK; t++) { int bi = -1, bj = -1, bs = 0;
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (pc[i][j] > bs) { bs = pc[i][j]; bi = i; bj = j; }
      if (bi < 0) break; pc[bi][bj] = 0; u64 x = (1ULL << bi) | (1ULL << bj);
      if (!inspan(s, x)) cand[nc++] = x; }
    // top clause parities (3-var active clauses by activity)
    double csc[MAXM]; int taken[MAXM]; u64 cf[MAXM];
    for (int i = 0; i < nact; i++) { int d; taken[i] = 0; csc[i] = 0; u64 U = actU[i];
      while (U) { csc[i] += sc[__builtin_ctzll(U)]; U &= U - 1; } cf[i] = clause_form(s, i, &d);
      csc[i] += 1e12 * (3 - d); if (!cf[i] || __builtin_popcountll(cf[i]) < 2) taken[i] = 1; }
    for (int t = 0; t < TOPK; t++) { int b = -1;
      for (int i = 0; i < nact; i++) if (!taken[i] && (b < 0 || csc[i] > csc[b])) b = i;
      if (b < 0) break; taken[b] = 1; cand[nc++] = cf[b]; }
  }
  int r0 = rank_of(s); double best = -1; u64 bx = cand[0];
  for (int c = 0; c < nc; c++) {
    int d[2];
    for (int b = 0; b < 2; b++) { St t = *s; int ok = add_eq(&t, cand[c], b) >= 0 && propagate(&t);
      if (!ok) { return cand[c]; } d[b] = rank_of(&t) - r0; }
    double scr = (double)d[0] * d[1] * 1024 + d[0] + d[1];
    if (scr > best) { best = scr; bx = cand[c]; }
  }
  return bx;
}
static u64 choose(const St *s) {
  switch (heur) {
    case 0: return best_var(s);
    case 1: return best_pair(s);
    case 2: return best_clause(s);
    case 3: return lookahead(s, 1);
    case 4: return lookahead(s, 0);
  }
  return best_var(s);
}
St model;
// returns 1 SAT, 0 UNSAT, -1 limit
static int search(St *s) {
  nodes++;
  if (!propagate(s)) { leaves++; return 0; }
  if (nact == 0) { model = *s; return 1; }
  if (leaves >= limit) return -1;
  u64 x = choose(s); nbr++; if (__builtin_popcountll(x) > 1) npar++;
  // need to re-propagate state info? choose used act arrays; fine (lookahead clobbers them but we're done)
  for (int b = 0; b < 2; b++) {
    St t = *s; int r = add_eq(&t, x, b);
    if (r < 0) { leaves++; nodes++; continue; }
    int res = search(&t); if (res != 0) return res;
  }
  return 0;
}
int main(int argc, char **argv) {
  const char *H[] = {"var", "pair", "clause", "linla", "varla"};
  for (int i = 0; i < 5; i++) if (argc > 1 && !strcmp(argv[1], H[i])) heur = i;
  if (argc > 2) limit = atoll(argv[2]);
  char line[4096]; int c = 0;
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] == 'c') continue;
    if (line[0] == 'p') { sscanf(line, "p cnf %d %d", &n, &m); continue; }
    char *p = line; int lit, off, k = 0;
    while (sscanf(p, "%d%n", &lit, &off) == 1) { p += off; if (lit == 0) break;
      if (k >= K) { fprintf(stderr, "clause too long\n"); return 2; }
      cvar[c][k] = abs(lit) - 1; cneg[c][k] = lit < 0; k++; }
    if (k) { clen[c] = k; c++; }
  }
  m = c; if (n > 64) { fprintf(stderr, "n>64\n"); return 2; }
  St s; memset(&s, 0, sizeof s);
  clock_t t0 = clock();
  int r = search(&s);
  double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
  if (r == 1) {  // extract point: free vars 0
    int a[64]; for (int v = 0; v < n; v++) a[v] = 0;
    for (int v = 0; v < n; v++) if ((model.piv >> v) & 1)
      a[v] = ((model.rhs >> v) & 1) ^ (__builtin_popcountll(model.row[v] & ~(1ULL << v)) & 0) ;
    // free vars are 0, so pivot var = rhs (other bits are non-pivots = 0)
    for (int cc = 0; cc < m; cc++) { int ok = 0; for (int i = 0; i < clen[cc]; i++) if (a[cvar[cc][i]] != cneg[cc][i]) ok = 1;
      if (!ok) { printf("BADMODEL\n"); return 3; } }
    printf("SAT %lld %lld %.3f\nv", leaves, nodes, secs);
    for (int v = 0; v < n; v++) printf(" %d", a[v] ? v + 1 : -(v + 1)); printf(" 0\n");
  } else printf("%s %lld %lld %.3f\n", r == 0 ? "UNSAT" : "LIMIT", leaves, nodes, secs);
  fprintf(stderr, "parity_branches %lld / %lld\n", npar, nbr);
  return 0;
}

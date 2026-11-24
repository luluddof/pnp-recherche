# G17a — Adversarial review of G14b (OpenAI Appendix D, tensor-contraction route)

Agent report (AI, same model family as the prover; not an independent human review).
Scripts: `rev1.py`–`rev5.py` (`rev2.py` timed out and has no saved output; its toy checks rerun in
`rev4.py`/`rev5.py`). Outputs: `rev1_out.txt`, `rev3_out.txt`, `rev4_out.txt`, `rev5_out.txt`.

## (a) Reconstruction: faithful
Checked item by item against the source; Dias–Lima App. B reads the same constants.
- Both inferred superscripts are right: rendered page 58 shows L = 2^a and D_arr = r^{2a}.
- D.3 read in full: it is only a reduction to Theorem D.1 (2⌈n^{1/4}⌉ group vertices, wires index
  monomials, s = nν, P0 = ℓ = n+1) and does not touch the D.4 mechanics.
- Minor discrepancies, none changing an exponent:
  1. The recovery loss has a second term, evaluation at −t/α (≤ 4M(a+1)), which model VT omits.
  2. Mode coefficients are divisible by 2^{w_Z + a − i}, not just 2^{w_Z}.
  3. "ν bounded below" is automatic (≈ 1.97 ≤ ν ≤ 6 by Lemma D.3), not an extra hypothesis.
  4. Theorem B's counting half overlaps Dias–Lima Prop. 4.1; the taper optimisation is the new part.

## Models
- DF is a faithful one-level abstraction. On toy packets, verified exactly: K_0 = I, (D.9) holds
  entrywise, K_t has rank 1 across groups.
- VT contains App. D, but its axiom "p_t ≥ −v2(E_t(0))" is **refuted** as a necessity: it is a property
  of formula (D.9), not of the data. Exact 2-adic computation: at a=1, M=2, q=40 the formula loss is 40
  but 10 bits suffice; at a=2, M=3, q=(9,5) loss_0 = 26 but (p_0,p_1) = (13,13) suffices.
- Outside VT: other linear recovery rules, non-Cauchy or multivariate deformations, odd-prime grids,
  sparse packings, and treating modes as wires or re-summing β instead of enumerating.

## (b) Theorems
- **Theorem A: confirmed.** Rank is over the field Frac(O); det ≡ ±1 mod 2 gives det ≠ 0, so
  "unimodular" is not needed. N ≥ 2^{(2−o(1))M/4^u} is correct but very weak: at the root
  4^u ≈ s^{1/65536}, an exponential gap in u from 2^{O(M/√u)}. It says nothing about √log.
- **Theorem B: partly.**
  - The inequality max Λ_i ≥ (2−1/L)M/(4u), the (D.10) check, the count 2^{√(np)−n} and the LP all hold
    (re-derived and re-run).
  - It is a tightness statement about App. D's own analysis and enumeration, not a lower bound for a
    family.
  - Link (i): the prover's sketch fails, since the paper's taper is about 3× over-determined and the
    witness degree is too high. Replacement supplied by the reviewer: F_m = ∏((T−q)/(−q))^{m_q} with
    tapered m, h′ ≈ M/u. It has degree ≤ M(2L−1), F_m(0) = 1, is integral on all discs, and all its
    Taylor coefficients at t are divisible by 2^{g_t}, with g_t ≥ (0.5–2)·M·D_t/u. So any recovery rule
    on Lemma D.5's hypothesis class needs p_t > g_t at some used node. This matches the exact minimum in
    toys (10, 19, 28 bits).
  - For the actual kernel span, the precision equals the black-box value in full-packet toys but is far
    smaller for a sparse packet (4 vs 14 bits), so it is unproved in general.
  - Link (ii) is open. Obvious skips (degree cap, redundant modes b ≥ D/2, true weights) do not change
    the exponent, by the reviewer's estimate, with a small numeric spot-check.
  - Link (iii) is harmless: a one-level count already lower-bounds the work.

## (c) Honest assessment
- Claimable today:
  - Theorem A, a standard flattening bound.
  - Lemma D.5's precision Θ(M·D_t/u) is optimal up to a constant of about 12 for black-box recovery.
  - App. D as written, with any valuation taper, enumerates 2^{Θ(M/√u)} tuples.
- Missing before "√log is optimal within family X":
  - a rank-type bound in DF, where 2^{M/4^u} is proved against a 2^{M/√u} target;
  - realisability of the witness by network kernels;
  - non-vanishing or incompressibility of the β-sum.
- The potential is real but modest: a tight analysis plus one clean black-box lemma, not a barrier.

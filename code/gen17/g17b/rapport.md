# G17b — Article 1 (pre-v3) checked against the OpenAI preprint itself

Agent report (AI, same model family as the authors of the paper; not an independent human review).
Source checked: OpenAI preprint "The Exponential Time Hypothesis Is False", revision of 25 October 2026
(not copied into this repository). Paper checked: `resultats/brouillons/pre-v3/pre_v3.tex`.
The agent's work folder was named `g18a` (numbering mistake at launch).

**Summary.** No assumption of pre-v3 is refuted by the preprint. Main defects: the object of Section 6
(the preprint copies sparsification leaves, not F), and a missing statement that our bound is slack on
the preprint's own design.

## Q1 (H*) — copy rule MATCHES; where it is applied DIFFERS
§2.8, p.18: "If a variable x occurs s ⩾ 2 times in this 3-CNF, replace its occurrences by distinct copies
x1,…,xs, preserving their signs, and add Tovey's implication cycle (¬x1∨x2)∧···∧(¬xs−1∨xs)∧(¬xs∨x1). …
Each copy appears in one clause inherited from the 3-CNF and two cycle clauses, hence in exactly three
clauses. Variables with at most one occurrence are kept unchanged, including unused declared variables."
- All variables with ≥2 occurrences are copied, on a cycle; fresh variables arise only from OR gates
  "for each longer clause", so none for a 3-CNF.
- This is a full copy expansion, so (E3), hence (E3′), holds, with ν₀ = N₀.
- But §2.8 copies "after disjoint sparsification", giving Ĝ for each leaf G. "F⋆" is Dias–Lima's
  notation; it equals Ĝ only when the tree is trivial.

## Q2 (leaves) — MATCHES Appendix A
(A.1)–(A.2): branch when "deg_c(H) ⩾ θ_{c−|H|}", θ_j = M^j, children "F_A = G∧H, F_B = …"; "do not
propagate unit clauses". The remark's leaf is the all-A-children leaf. The claim that few clauses are
deleted w.h.p. is NOT FOUND: the preprint has no random-formula analysis.

## Q3 (Definition 5.1) — faithful generalisation of §2.3–2.5 and Definition 2.1
"⟨P_{e,v}(0),X⟩ = X_v, deg P_{e,v} ⩽ R"; condition (4): "P_{e,v}(t)|S_h lies in the row space of C_h for
every h ≠ g". Definition 5.1 drops integrality, the row bound D and |S_g| ⩽ 2n/r.
- Parameters (2.10)/(2.19): r = ⌊√λ/ν²⌋ groups, Q = ⌊λ/(128ν)⌋ points 1..Q, R = 4⌈√λ⌉, q = ⌊3nR/Q⌋
  (Dias–Lima write ⌈·⌉), at most 2τn = 2n/(rλ) rows, D = r^{2Q} ⩽ (2n)^{1/128}.
- (2.1) has "n ⩽ M ⩽ 3n" via dummy clauses, which is our Lean hypothesis.
- The design is not variable-uniform: "Different clauses may use different replacement linear forms for
  the same variable."

## Q4 (cost) — product-type; DP-type NOT FOUND
§2.4: "Try every measurement tuple y ∈ Y := ∏{−D,…,D}^{c_g}", time "(2D+1)^{2τn} 2^{2n/r}"; (2.11)
multiplies the three costs. So 2^{ρ_t} is paid. §3.2 warns: "A rank bound alone would not control the
number of possible measurements".

## Q5 (routes)
"Appendix D proves both bounds in Theorem 1.1 by tensor contraction"; "The rational interpolation
argument in Section 6.4 gives the weaker exponent O_k(ρ(n))" (parity form, original variables).
Clause-product gives (2.12), 2^{O(n(log log n)²/√log n)}, on 3-occurrence 3-CNF, and
2^{O_k(n(log log n)^{k+1}/(log n)^{1/(2k)})} in general (§2.8, "A weaker algorithm"; "Taking k = 3 already
refutes ETH").

## Q6 / corrections for v3, by importance
1. Abstract, "that this is the preprint's formula … was not checked" → "the copy rule is that of §2.8
   (checked); the preprint applies it to sparsification leaves, so F⋆ is Dias–Lima's formula".
2. §6.1, "We did not read…" and (H*) → quote §2.8; "kept or dropped" → "kept". Remark on variants: (a)
   and (b) do not occur for 3-CNF; (c) → "sparsification always precedes". Promote the leaf remark to the
   main statement; "first-child" → "A-child"; "consists of" → "contains".
3. Add to the introduction and scope: the preprint's design has s_t ≈ n/r and
   q = Θ(n log log n/√log n), so Theorem 2 holds through q alone and constrains only small-group
   redesigns.
4. Scope row for Cor. 1(3), "which we did not examine" → "it is not variable-uniform (§2.3)".
5. Scope question 3 and §3.2 → state the answer with the §2.4 quote and the rank caveat.
6. Related work, "According to [DL]…" → cite Theorem 1.1, Appendix D, §6.4 and §2.8 directly; at bounded
   density the headline bound is D.2 applied to F.
7. Lean: note that n ⩽ M holds in (2.1). Footnote the ⌊·⌋/⌈·⌉ discrepancy.

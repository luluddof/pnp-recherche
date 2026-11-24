# Referee report (G8R)

**Manuscript.** "Measurements cannot absorb the cut: a determinant bound for clause-product interpolation schemes" (`resultats/note_limits_s9.md`, about 290 lines of Markdown with LaTeX, about 6 pages).
**Answers the question in.** F. Dias, Y. Lima, "Limits and Constants of OpenAI's Claimed ETH Refutation", arXiv:2611.00517v1, §9 (cited below as [DL]).
**Recommendation.** **Major revision.** The mathematics is correct as far as I could check it. The revision is needed for two scope claims (one of them wrong) and for a presentation that is far from arXiv standard (details in §7).

---

## 1. Summary of the claimed contribution

[DL, Def. 5.1] defines a "clause-product interpolation scheme". Each occurrence $(e,v)$ is deformed into a degree-$R$ vector $P_{e,v}(T)$ with $P_{e,v}(0)=e_v$. $H(0)=\#F$ is recovered from Hermite data at nonzero points $t$ (order $q_t$). At each point a *design* partitions the variables into groups, and a clause may see a group other than its own only through a subspace $\mathcal R_{t,h}$ of measurement rows. [DL] prove $2^{\Omega(\sqrt n)}$ on expanders. Their technique cannot go further, because a measured group of size $s$ can "absorb" all of its $\approx bs$ foreign occurrences. They leave open (§9) whether measurement rows can absorb the cut, and they give only a pairwise $2\times2$-minor remark for one-dimensional groups.

The note's contributions, in my words:

1. **Witness determinant (Lemma 1).** Pick one occurrence per variable. The resulting $n\times n$ matrix $\Phi^\omega(T)$ has $\Phi^\omega(0)=I$ and $\deg\det\Phi^\omega\le nR$. At each point its kernel has dimension at least a combinatorial *deficiency* minus $\rho_t$. Since $\det\Phi^\omega$ vanishes at $t$ to order at least $\dim\ker\Phi^\omega(t)$, the deficiencies summed over points (minus $\rho_t$) total at most $nR$.
2. **Averaging and anticoncentration (Lemma 2).** For a uniformly random witness, the expected deficiency at $t$ is at least $c_0/(2\sqrt{2s_t})$ times a weighted count of *split* variables (variables whose clauses go to at least two groups), with $c_0=(1-e^{-2})/4$.
3. **Split variables form a separator (Lemma 3)**, unless one group's clauses already hold more than $2n/3$ of the variables.
4. **Theorem B.** Some point satisfies $\rho_t+q_t+1>\min\bigl(c\,\mathrm{sep}(F)/\sqrt{s_t},\ 2n/3-s_t\bigr)$. So rows can absorb at most a $\sqrt{s_t}$ factor of the cut.
5. **Theorem A.** For variable-uniform deformations there is no $\sqrt s$ loss: $\rho_t+q_t+1>\min(\mathrm{sep}(F),2n/3-s_t)$.
6. **Lemma 4.** $\mathrm{sep}(F)=\Omega(n)$ w.h.p. for random 3-CNF with $\Delta>3$, adapted from [DL, Lemma 8.2].
7. **Corollary C.**
   - The §9 target design (groups of size $\tilde O(\sqrt n)$, $\tilde O(\sqrt n)$ rows and Taylor order) does not exist on such formulas. Some point needs $\rho_t+q_t=\Omega(n^{3/4}/\mathrm{polylog})$.
   - In the product-type cost model every scheme costs $2^{\Omega(n^{2/3})}$.
   - Variable-uniform schemes cost $2^{\Omega(n)}$.

**Do the claims match the proofs?** For the theorems, lemmas and Corollary C(1–3) as stated in §2: yes. Three framing claims do not match:

- **The abstract's "improving the paper's $2^{\Omega(\sqrt n)}$".** [DL, Thm 8.3] is about the scheme applied to the preprint's bounded-occurrence formula $F^\star$ and to sparsification leaves. The note's bound is for the scheme applied to $F$ itself, or to bounded-occurrence formulas with linear separators. The note concedes in §5 that $\mathrm{sep}(F^\star)\ge cn^\star$ is unverified, but the abstract does not say so. The fair comparison is with [DL, Thm 8.1 + Lemma 8.2] on $F$, which also give $2^{\Omega(\sqrt n)}$. Against that comparison, the improvement is real.
- **"DP-type model: open" (abstract and §5).** This is wrong as stated; see §2.9 below. In [DL]'s literal DP-type model the cost on expanders is $2^{\Theta(\sqrt n)}$ by a trivial scheme. So no lower-bound improvement is possible in that model, and the gap is not open there.
- **The title.** "Measurements cannot absorb the cut" overstates Theorem B, whose own title is "Measurements absorb at most a $\sqrt s$ factor of the cut". Whether a $\sqrt s$ fraction can be absorbed is left open (§5, tightness).

---

## 2. Correctness

I re-derived every step listed below. All constants were recomputed:

- $c_0=0.216166$;
- $c_3=c_0\cdot 2/(9\sqrt2)=0.0339672$, which equals $2(1-e^{-2})/(36\sqrt2)$;
- $\sqrt{4/7}/2=0.37796$.

I also re-ran `g6a/anticonc_check.py` (worst ratio $0.7397$, as reported) and `g6r/lemma1_check.py` (see §2.10).

### 2.1 Fact R (§3.1): correct
$f=\prod_t(T-t)^{q_t+1}$ has zero Hermite data and $f(0)\ne0$. If $\deg f\le MR$, this contradicts the exactness of $\Lambda$.

*Attribution.* This is the argument of [DL, Lemma 6.2] (and of their (R1)), applied to a different polynomial. The note should say so.

### 2.2 Lemma 1 (§3.2): correct
- $\Phi(0)=I$ holds, since row $v$ is $P_{\omega(v)}$ with $\omega(v)$ an occurrence of $v$, and $P_\bot=e_v$.
- The column-block rank bound $\operatorname{rank}\Phi(t)\le\sum_h\operatorname{rank}B_h$ holds.
- Rows outside $I^\omega_{t,h}$ restricted to $S_h$ lie in $\mathcal R_{t,h}$, by (3) or because $e_v|_{S_h}=0$ in the $\bot$ case.
- $(y-r)_+\ge y_+-r$ holds.
- The divisibility step (change of basis, $(T-t)^k\mid\det$, distinct points) is standard and correct.

**Observation (simplification, not an error).** The hypothesis $M\ge n$, used in Theorems A and B, Corollary C and the "Note" after Lemma 4, is superfluous. Rows of non-occurring variables are constant, so $\deg\det\Phi^\omega\le R\cdot\#\{v:d_v\ge1\}\le MR$, which is all Theorem B needs. With this change, Lemma 1 can state $\le n_{\rm occ}R$, the hypothesis disappears, and the paragraph defending $M\ge n$ for random 3-CNF is no longer needed.

**Missing citation.** "Order of vanishing of $\det\ge\dim\ker$" is classical (local Smith form, or the multiplicity lemmas of the polynomial method). Cite a standard source.

### 2.3 Lemma 2 (§3.3): correct, constant verified
- **Deficiency.** $\operatorname{def}=\frac12\sum_h|X_h-a_h|$ follows from $\sum_hX_h=n$, since each variable has exactly one home group.
- **Symmetrisation.** $\mathbb E|X-a|\ge\frac12\mathbb E|X-X'|$ follows from the triangle inequality. $\mathbb EZ^2=2V$ is correct.
- **Step (i).** Conditioned on $K\neq\emptyset$, $Z$ is a sum of $|K|$ independent fair signs, so $\Pr(Z=0\mid K)\le\binom{2j}{j}4^{-j}\le\frac12$. Also $\Pr(K\ne\emptyset)\ge1-e^{-2V}$. For $V<1$, $(1-e^{-2V})/V$ is decreasing. Hence $\mathbb E|X-a|\ge c_0V$. Correct.
- **Step (ii).** $\mathbb EZ^4=\sum\tau_v+3\sum_{u\ne v}\tau_u\tau_v\le2V+12V^2$. The Hölder form $\mathbb EZ^2\le(\mathbb E|Z|)^{2/3}(\mathbb EZ^4)^{1/3}$ is correct (exponents $3/2$ and $3$). For $V\ge1$, $(2V)^{3/2}/(2V+12V^2)^{1/2}\ge2\sqrt2V/\sqrt{14V}=\sqrt{4/7}\sqrt V$. Correct.
- **Passage to $V_h/\sqrt{2s}$.** Both cases are correct. The case $\mu_h>2s$ needs $\frac12\ge c_0/\sqrt{2s}$, which holds for $s\ge1$; say so.
- **Weights.** $1-\sum_h\pi_{v,h}^2\ge w(d_v)$ for split $v$ is correct: the extremal split is $(d-1,1)$.
- **Remarks.**
  - The binding constraint on $c_0$ is step (i) at $V=1$. Step (ii) gives $0.378$, and the experiments suggest about $0.74$. A sharper constant is easy to get but not needed.
  - The fourth-moment / Hölder lower bound $\mathbb E|Z|\ge(\mathbb EZ^2)^{3/2}/(\mathbb EZ^4)^{1/2}$ is classical (Littlewood; Berger 1997). Anticoncentration of Bernoulli sums is classical too (Erdős–Littlewood–Offord, Kolmogorov–Rogozin). Cite these.
  - Notation clash: $X$ is the evaluation set, and $X$, $X'$, $X_h$ are random variables in this proof.

### 2.4 Lemma 3 (§3.4): correct
- An edge between two non-split variables forces equal home groups, so components of $\Pr(F)-\Sigma_t$ stay inside one $J_{t,h}$.
- Case 2: $\frac12\sum_h|\mu_h-a_h|=\sum_h(\mu_h-a_h)_+\ge\mu_{h_0}-a_{h_0}>2n/3-s_t$. Correct.
- Minor gap: "a non-split variable has a well-defined home group" relies, for $d_v=0$, on the home-group convention for $\bot$. Say so in the definition of $J_{t,h}$.

### 2.5 Theorem B (§3.5): correct
- Jensen applies to the convex map $(\cdot)_+$.
- The chain $\le nR\le MR<\sum(q_t+1)$ is valid.
- Pigeonhole gives a point with $\mathbb E\operatorname{def}_t<\rho_t+q_t+1$.
- The case split covers all cases.
- $\Sigma_t\subseteq\{d_v\ge2\}$, so $W(F)$ applies.
- $c_b=c_0(b-1)/(\sqrt2b^2)$ is correct.

### 2.6 Theorem A (§3.6): correct
- The unconstrained rows of block $h$ are exactly $J_{t,h}$: a split variable is constrained in every block, and a non-split one in every block except its home group.
- $\sum_h|J_{t,h}|=n-\sigma'_t$ holds.
- In Case 2, $\sum_h\min(a_h,|J_{t,h}|)<s_t+n/3$.
- Nice and clean. The phrase "global linear change of variables $x\mapsto\Phi(T)x$" reuses the symbol $\Phi$ with a different meaning; rename it.

### 2.7 Lemma 4 (§3.7): correct, given [DL, App. A]
- The adaptation from a treewidth separator (components $\le n/2$) to a $2/3$-separator is sound. The construction of $A$ with $n/3\le|A|\le2n/3$ works because each added component is $<n/3$ in the second case.
- I checked [DL, Lemma A.1]'s proof: $\alpha^3+(1-\alpha)^3\le1/3$ on $[1/3,2/3]$, and the union bound $2^{n-\Delta n/3}\to 0$. It supports the reading "at least $m/3$ clauses meet both $S$ and its complement". The source's wording is "shares at least $m/3$ clauses with its complement", and the note's paraphrase is the right interpretation.
- Edge case: $|Y|=0$ makes $|Y|\log(en/|Y|)$ undefined. It is excluded directly by (A.1), but this should be said.
- The bound on $W(F)$ (Markov on degrees, $D=\lceil6\Delta/c_\Delta\rceil$) is correct.
- "$M=3\Delta n$ deterministically" is correct in [DL]'s model of three distinct variables per clause. State the model explicitly.

### 2.8 Corollary C (§3.8): correct, with small gaps
- **C.1.** Correct.
- **C.2.**
  - $\binom nj\ge2^j$ for $j\le n/2$ is correct.
  - The split at $s_t\ge n^{2/3}$ is correct.
  - The final $\Omega(n^{2/3})$ is correct.
- **C.3.** The case $s_t\ge n/6$, where the cost is at least $2^{n/6}$ directly, is left implicit in "Same argument". Write it out.
- **Balanced reformulation (end of §2 and §3.8).** Two corrections:
  - It holds "at some point $t$", not for the design in general.
  - It needs the first term of the min to be the active one, i.e. $\gamma n/\sqrt{s_t}\le 2n/3-s_t$. For large $n$ that means roughly $s_t\le n/2$, i.e. $g_t\ge3$.

### 2.9 Scope error: the DP-type model is degenerate, not "open"
§5 says: "In the DP-type model, the gap between $2^{\Omega(\sqrt n)}$ [DL, Thm 8.3] and the preprint's $2^{O(\rho(n))}$ remains open." That is false in [DL]'s model as written. The following scheme shows it:

- constant deformations $P_o\equiv e_v$, so $R=0$ and $H\equiv\#F$;
- $X=\{1\}$, $q=0$, $\Lambda(f)=f(1)$, which is exact up to degree $MR=0$;
- at the single point, $\sqrt n$ groups of size $\sqrt n$, all fully measured ($\mathcal R_{t,h}=\mathbb Q^{S_h}$, so (3) holds trivially).

If $R\ge1$ is required, take $MR+1$ points with $q_t=0$ and the same design. This scheme satisfies Def. 5.1:

- $s_t=\sqrt n$;
- $w_t\le\#\text{groups}-1<\sqrt n$;
- the Taylor charge is $\binom n0=1$.

So the DP-type charge is $2^{\sqrt n}$. Hence [DL]'s $2^{\Omega(\sqrt n)}$ is *tight in their DP-type model*, on every formula, and nothing above $2^{\sqrt n}$ can be proved in that model. This observation also appears, unreviewed, in the project journal (entry G7b, "Prop. O1"). I verified it independently against the text of Def. 5.1 and §5.2. The model under-charges because the $2^{w_t}$ term counts one bit per measured group, whatever its number of rows.

Consequences for the note:

- §5 and the abstract must say that the DP-type question is degenerate as posed, not open. A meaningful DP version needs a corrected charge, for example a width weighted by rows.
- The cost significance of Corollary C(1) lies in the product-type model. In the DP-type model the §9 *cost* target is achieved trivially, by a design with $\sqrt n$ rows per group instead of polylog many.
- The note could usefully point out this flaw in [DL]'s §9 framing; it strengthens the case for the product-type result.

### 2.10 Computational section (§4): factual inaccuracies
- **Item 4** says the reviewer's check covered "6 seeds times 6 witness maps". Re-running `g6r/lemma1_check.py` shows that seeds 0, 4 and 5 return `None`: the linear system for the deformations was infeasible and no scheme was built. Only 3 seeds × 6 witness maps were checked. For those, the maximum of $\sum_t(\operatorname{def}_t-\rho_t)_+$ is 1 or 3, against $nR=18$. The test is weak, and the text overstates it.
- **Item 1.** In every run of `theorem_checks.txt`, $\sum_t(\mathbb E\operatorname{def}_t-\rho_t)_+=0.00$. So the averaging step of Theorem B is tested only trivially. The caveat in the note mentions large $\rho_t$ but should state this consequence explicitly.
- **Item 1.** `run_checks.py` also configures $(n,g,|X|,R)=(36,6,9,3)$, but the output file contains only $n=24$ and $n=30$. Either the run did not finish or the output was truncated. Explain this, or remove the configuration.
- The note gives no commands, dependencies (sympy, Python version) or run times, and it points to repository-internal paths (`code/gen6/...`) that a reader cannot reach.

---

## 3. Faithfulness to the source

- **Definition 5.1.** It is reproduced accurately in (a), (b), (c) and (3), and condition (3) matches the source exactly. Omissions:
  - **The edges of the interaction graph are not defined.** The note writes only "The interaction graph $\mathcal I_t$ has the groups as vertices, and $w_t=\operatorname{tw}(\mathcal I_t)$". The source continues: "for each clause $e$ the groups $\{\alpha_t(e)\}\cup\{h\ne\alpha_t(e):P_{e,v}(t)|_{S_{t,h}}\ne0$ for some $v\in I_e\}$ form a clique". Since $w_t$ enters the cost model, this must be restored.
  - The sentence that the preprint's clause-product form is an instance, "with $\mathcal R_{t,h}$ the row space of its integer measurement matrix $C_{t,h}$ and $\Lambda$ Hermite interpolation with $q_t=\lceil 3nR/|X|\rceil$", is omitted. It is useful context, because it shows the theorem applies to the preprint's own construction.
  - "clauses of at most $k$ literals" is omitted. This is harmless.
- **Cost model §5.2.** The formula is accurate, and so is the "additionally $\kappa(t)\ge2^{\rho_t}$". Two omissions:
  - [DL]'s caveat: "they are statements about running time only to the extent that the model's charges are paid".
  - The explicit exclusion: "An evaluator that also counted inside groups by dynamic programming would not pay $2^{s_t}$, and our model does not cover it."

  Both matter, because the note's claim (§1.2, §5) that $s_t=\tilde O(\sqrt n)$ is *necessary* for the §9 cost depends entirely on the $2^{s_t}$ charge.
- **The §9 quotation.** It is verbatim and correct. Points to fix:
  - The elided sentence "(Occurrences of the same variable can be absorbed together at no cost.)" is precisely why the note uses one witness occurrence per variable. Restore it.
  - $p_{e,v,v'}$ is undefined in both papers. Define it as $(P_{e,v}(T))_{v'}$.
  - It would help to also quote the preceding sentence ("Our lower-bound technique cannot exceed $2^{\tilde O(\sqrt n)}$: its only constraint is $(w+1)s+q\gtrsim\mathrm{tw}$"). It explains why a new idea was needed.
- **Paraphrase in the abstract.** "They resolve the case of one row per group by a $2\times2$-minor argument" overstates [DL]. They observe that a pair is absorbed at at most $2R$ points. They do not turn this into a bound even for one row per group. The sentence "Turning this into a bound for groups with several measurement rows is open" is ambiguous on the one-row case. Write instead: "they give a pairwise $2\times2$-minor observation for one-row groups".
- **"In §9 they ask whether…".** [DL] state a direction ("one must show that measurements cannot absorb the cut"), not a formal question. The note's interpretation (all points, $s_t,q_t,\rho_t=\tilde O(\sqrt n)$) is reasonable, and it is stated explicitly, which is good.
- **Theorem 8.3.** It is described correctly ($F^\star$ and sparsification leaves).

---

## 4. Novelty and significance

**New relative to [DL]:**

- The $n\times n$ witness determinant. It generalises the $2\times2$ principal minor of $\Phi$ that [DL] use in §9. It also improves the degree bookkeeping relative to [DL, Lemma 6.2], since $\deg\le nR$ where $\prod_oc_o$ has degree up to $MR$.
- The deficiency / rank inequality per column block.
- The random-witness averaging with an explicit anticoncentration constant.
- The "split variables form a separator" lemma.
- Theorem A, which shows that variable-uniform deformations gain nothing from measurements.

Together these give the first bound in which *rows* are charged rather than *measured groups*. Compare [DL, Thm 7.1(ii)]: deleting $q_t+b\mu_ts_t$ edges isolates the groups, which gives $\mu_ts_t\gtrsim n$ on expanders. The note gives $\rho_t\sqrt{s_t}\gtrsim n$. This is the right one-line comparison, and the note should state it explicitly; at the moment it never mentions Thm 7.1(ii).

**Significance (honest assessment):**

- **Positive.**
  - It settles, negatively and with clean proofs, the specific "from below" question of [DL, §9] in the product-type model, the model that charges what the preprint's evaluation actually pays.
  - It moves the product-type barrier on expanders from $2^{\Omega(\sqrt n)}$ to $2^{\Omega(n^{2/3})}$ for the scheme applied to $F$.
  - The determinant-plus-anticoncentration technique is short and reusable.
- **Limiting.**
  1. It is a barrier for one abstraction of one route. It says nothing about the preprint's headline bound (tensor contraction) or its parity form, both outside Def. 5.1, as the note correctly says.
  2. It is a lower bound in a *charging* model, not a running-time lower bound.
  3. No upper bound is given in the product-type model, so whether $n^{2/3}$ is the truth is open.
  4. It does not yet apply to the preprint's actual route ($F^\star$ or sparsification leaves).
  5. In the DP-type model it adds nothing, and cannot, because the model is degenerate (§2.9).

  Overall this is a solid, modest negative result: appropriate for an arXiv note or a workshop paper, not a major venue on its own.

**Prior work.** I am not aware of this result in the literature. The note's novelty rests on an AI search; it should say so with the same caveat [DL] use ("rests on searches, not a systematic survey").

---

## 5. Structure and presentation compared with [DL]

| [DL] component | Note | Gap |
|---|---|---|
| Title, authors, affiliations, date | Title only; status line; no author | **No human author takes responsibility.** arXiv requires human authors; AI tools cannot be authors. |
| Abstract | Present, about 250 words | Too long and formula-heavy. $c_b$ is used before $b$ is defined, and only $c_3$ is given. Overclaims (see §1 and §3). |
| 1 Introduction with numbered contributions | **Missing.** §1 "Setting" goes straight to definitions | Add context (the preprint, [DL]'s barrier, why §9 matters) and a numbered list of contributions (Lemma 1, Thm B, Thm A, Cor C, scope remark on the DP model). |
| 2 Related work | **Missing** | Cover: [DL] Thm 7.1 and 8.1–8.3; polynomial-method multiplicity and determinant lemmas; anticoncentration (Littlewood–Offord, Kolmogorov–Rogozin, moment method); balanced separators (Lipton–Tarjan) and their relation to treewidth; expansion of random CNF (Ben-Sasson–Wigderson, Ben-Sasson–Impagliazzo); limits-of-methods literature (as in [DL] §2). |
| 3 Plain-terms explanation | Missing | A half-page intuitive overview of "one row per cut variable vs $\sqrt s$" with a toy example would help non-specialists. |
| Technical sections with numbered results | Present (§1–3) | Numbering is mixed: Fact R, Lemmas 1–4, Theorems **B then A**, Corollary C. Use sequential numbering (Lemma 2.1, Theorem 2.5, …) and state A before B or merge them. Statements are in §2 and proofs in §3; acceptable, but give each proof a back-reference. |
| Notation | §1.3 list | No notation table. Clashes: $X$ (point set) vs random variables $X,X',X_h$; $\Pr(F)$ (primal graph) vs $\Pr(\cdot)$ (probability) in §3.3; $w(d)$ vs $w_t$; $\Phi$ (witness matrix) vs $\Phi$ (Thm A change of variables); $\rho_t$ vs the preprint's $\rho(n)$ (acknowledged); $\sigma'_t$ has a prime with no unprimed $\sigma$; about eight constants ($c_0,c_b,c_3,c_\Delta,c'_\Delta,c^*_\Delta,\gamma,c$). |
| 10 Experiments with reproducible commands | §4 "Computational checks" | No commands, environment, versions or run times. Paths are internal. One factual error (item 4) and one vacuous test (item 1); see §2.10. No table of results. |
| Code availability / listing appendix | Missing | Give a public repository or an appendix listing, as [DL] App. C does, plus one command per check. |
| 11 Methodology and AI use | §6, before an internal "Issues found" section | Content is fine and candid. Move it near the end (before the references), add who takes responsibility, and state the novelty-search caveat. |
| 12 Conclusion | Missing (§5 "Scope and open problems" partly fills the role) | Add a short conclusion, and fix the DP-model statement there. |
| References with full citations | **One reference** ([DL]) | Needs a real bibliography: the preprint itself (OpenAI, arXiv id), [DL], and the standard tools listed above. |
| Appendices (proofs, code) | None | Optional: move Lemma 2's anticoncentration and Lemma 4 to an appendix to shorten the main line. |
| Figures and tables | None | A table comparing the bounds in [DL] and here (model × quantity charged × bound × formula family) would help. A schematic of $\Phi$ with column blocks would also help. |
| "Issues found (relative to the source notes and the journal)" | Present | Internal changelog that refers to unpublished notes and a journal. **Remove it from any public version**, and also the in-text pointers ("source notes round this down", "as observed by the adversarial review"). |

**Writing quality.** Sentences are clear, short and precise. Proofs are well organised and easy to check, which is a real strength. Length (about 6 pages) is appropriate for the content. The main presentation problems are the missing framing (introduction, related work, conclusion), internal project references, mixed numbering and notation clashes.

---

## 6. Required and suggested changes

### Required (must fix)
1. **DP-type scope (§2.9).** Replace "DP-type model: open" (abstract, §5) with a correct statement. In [DL]'s literal DP-type model a trivial scheme costs $2^{\sqrt n}$, so their $2^{\Omega(\sqrt n)}$ is tight in that model, and any improvement needs a corrected charge. Include the one-paragraph construction.
2. **Comparison with [DL, Thm 8.3].** In the abstract and Corollary C(2), say that $2^{\Omega(n^{2/3})}$ is for the scheme applied to $F$ (or to bounded-occurrence formulas with $\mathrm{sep}=\Omega(n)$), not to $F^\star$ or the sparsification leaves. The fair comparison is [DL, Thm 8.1 + Lemma 8.2]. Also compare explicitly with [DL, Thm 7.1(ii)] ($\mu_ts_t\gtrsim n$ vs $\rho_t\sqrt{s_t}\gtrsim n$).
3. **Overstatements.** Fix three:
   - the title, e.g. "…cannot absorb more than a $\sqrt s$ fraction of the cut";
   - the abstract's "they resolve the case of one row per group";
   - the abstract's "improving the paper's $2^{\Omega(\sqrt n)}$".
4. **Faithfulness.** Restore the definition of the edges of $\mathcal I_t$, and [DL]'s two cost-model caveats (charges must be paid; DP inside groups is not covered). Restore the elided §9 sentence on same-variable occurrences, and define $p_{e,v,v'}$.
5. **Computational section.**
   - Correct item 4 (3 of 6 seeds produced schemes, so 18 checks, not 36).
   - Say explicitly that $\sum_t(\mathbb E\operatorname{def}_t-\rho_t)_+=0$ in all runs, so Theorem B's averaging step is untested.
   - Explain the missing $n=36$ run.
   - Give commands, dependencies and a public code location.
6. **Authorship and responsibility.** Add human author(s) who take responsibility, or do not post to arXiv.
7. **Missing structure.** Add an Introduction with a numbered contributions list, a Related work section, a Conclusion, and a full bibliography. The bibliography must include the OpenAI preprint and the standard results used (determinant multiplicity, anticoncentration and moment bounds, separators, random-CNF expansion). Attribute Fact R's argument to [DL, Lemma 6.2].
8. **Internal material.** Remove the "Issues found" section and all references to "source notes", "journal" and "adversarial review", except as a neutral sentence in the AI-use statement.
9. **Notation clashes.** Fix $X$ (points vs random variables), $\Pr$ (primal graph vs probability), $\Phi$ (two meanings) and $w$ vs $w_t$.
10. **Small proof gaps.**
    - Corollary C.3: the case $s_t\ge n/6$.
    - Balanced reformulation: "at some $t$", with the condition $s_t\lesssim n/2$.
    - Lemma 4: the edge case $|Y|=0$, and the random-CNF model (three distinct variables) stated.
    - Lemma 3: the home group of non-occurring variables in $J_{t,h}$.
    - Lemma 2: $1/2\ge c_0/\sqrt{2s}$ stated.

### Suggested (nice to have)
1. Drop the hypothesis $M\ge n$ everywhere: $\deg\det\Phi^\omega\le n_{\rm occ}R\le MR$ always. This also removes the "Note" after Lemma 4.
2. Prove $\mathrm{sep}(F^\star)\ge cn^\star$ (or the weighted analogue) for the Tovey-cycle formula. It reduces to a degree-weighted version of Lemma 4, because components of $\Pr(F)-Y$ map into components of $\Pr(F^\star)-Y^\star$. That would make Corollary C apply to [DL, Thm 8.3(a)]'s setting and justify the comparison claim fully. Consider also the sparsification leaves (Thm 8.3(b)), which [DL] handle by deleting $O(\eta n)$ clauses; separators are similarly robust.
3. Add a comparison table of bounds (model, charged quantity, bound, family) and a figure of the block structure of $\Phi^\omega(t)$.
4. Add a short plain-language section (rows vs cut variables; why the loss is $\sqrt s$, namely fluctuations of a random witness).
5. Discuss tightness: give the best known product-type upper bound on expanders, or a candidate design with $\rho_t\approx n/\sqrt{s_t}$, or explain why the random-witness loss may be an artefact (a deterministic witness choice or a different minor could remove it).
6. Comment on what a corrected DP-type model would need (rows-weighted width), and point to it as the natural next question.
7. Improve $c_0$ (for example by using (ii) down to smaller $V$, or a sharper Littlewood–Offord bound), or simply state that the constant is not optimised.
8. Shorten the abstract to about 150 words, move the displayed inequality into the introduction, and simplify $c_3=(1-e^{-2})/(18\sqrt2)$.
9. Strengthen the reviewer test (`lemma1_check.py`): retry infeasible seeds, use larger $n$, and include instances with $\operatorname{def}_t>\rho_t$ at several points, so that the inequalities are tested away from trivial cases.
10. Add the novelty-search caveat to the methodology section, matching [DL] §11.

---

## 7. Recommendation: major revision

The core mathematics is correct and nicely compact. I checked Fact R, Lemmas 1–4, Theorems A and B and Corollary C line by line and found no error affecting the stated theorems, only small gaps that are easy to close. The witness-determinant plus anticoncentration argument is a genuine, if modest, advance on the question [DL] left open in §9: it charges rows rather than measured groups and gives $2^{\Omega(n^{2/3})}$ in the product-type model.

As a self-contained arXiv note or workshop submission, however, it needs a major revision:

- One scope statement is wrong ("DP-type model: open"). [DL]'s DP-type model is degenerate, and a trivial scheme already achieves $2^{\sqrt n}$ in it.
- The headline comparison with [DL, Thm 8.3] mixes $F$ and $F^\star$.
- The title and abstract overstate both the result and [DL]'s one-row remark.
- The computational section contains a factual overstatement and an untested step.
- Measured against the reference paper, the manuscript lacks an introduction with contributions, related work, a bibliography, a conclusion, reproducible code and a responsible author, and it still contains internal project material.

None of this requires new mathematics, except optionally the $F^\star$ separator bound. With these changes I would expect to recommend acceptance.

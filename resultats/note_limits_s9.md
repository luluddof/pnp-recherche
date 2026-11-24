# Measurements cannot absorb the cut: a determinant bound for clause-product interpolation schemes

*A note on the open question of §9 of Dias–Lima, "Limits and Constants of OpenAI's Claimed ETH Refutation" (arXiv:2611.00517).*

**Status.** This result was produced by AI agents and has **not yet been reviewed by human experts**. See §6.

---

## Abstract

Dias and Lima (arXiv:2611.00517, Definition 5.1) abstract the "clause-product" form of an interpolation algorithm for #SAT as a *scheme*. A scheme has three parts: a polynomial deformation of every clause occurrence; Hermite recovery of $H(0)=\#F$ from data at nonzero points; and, at each point $t$, a *design*. A design partitions the variables into groups, and a clause may see the groups other than its own only through a few linear *measurement rows*. On expanding formulas they show that such a scheme costs $2^{\Omega(\sqrt n)}$. In §9 they ask whether measurements can "absorb the cut". Concretely, they ask whether a design with about $\sqrt n$ groups, Taylor order $O(\sqrt n)$ and polylogarithmically many rows per measured group could exist on an expander; it would give cost $2^{\tilde O(\sqrt n)}$. They resolve the case of one row per group by a $2\times 2$-minor argument and leave the case of several rows open.

We answer the question negatively for groups of size $s_t=\tilde O(\sqrt n)$. The main tool is a *witness determinant*: choose one occurrence of each variable, and let $\Phi(T)$ be the $n\times n$ matrix whose rows are the chosen deformation vectors. Then $\det\Phi$ is a nonzero polynomial of degree at most $nR$, and its order of vanishing at $t$ is bounded below by a combinatorial *deficiency* of the design at $t$ minus the number $\rho_t$ of measurement rows. We average over random witnesses and use an anticoncentration bound with an explicit constant. The result is that at some point $t$,
$$\rho_t+q_t+1>\min\Bigl(c_b\,\frac{\operatorname{sep}(F)}{\sqrt{s_t}},\ \frac{2n}{3}-s_t\Bigr),\qquad c_3=\tfrac{2(1-e^{-2})}{36\sqrt2}\approx 0.0340,$$
where $\operatorname{sep}(F)$ is the $2/3$-vertex-separator number of the primal graph and $b$ bounds the number of occurrences of a variable. If the deformation depends only on the variable (not on the clause), no absorption at all is possible. Consequences on formulas with $\operatorname{sep}(F)\ge cn$, including random 3-CNF of density $\Delta>3$:

- the §9 target design does not exist;
- every scheme costs $2^{\Omega(n^{2/3})}$ in the paper's *product-type* cost model, improving the paper's $2^{\Omega(\sqrt n)}$ in that model.

The DP-type cost model is left open.

---

## 1. Setting

### 1.1 Definition 5.1 of Dias–Lima

We follow [DL, §5] (Dias–Lima, arXiv:2611.00517v1). Let $F$ be a CNF over $x_1,\dots,x_n$, and let $I_e$ be the variable set of clause $e$. The *occurrences* are $O=\{(e,v):v\in I_e\}$, and $M=|O|$. Write $d_v$ for the number of occurrences of $v$, and $b$ for an upper bound on all $d_v$. The *primal graph* $\Pr(F)$ joins two variables that share a clause. Let $\varphi_e$ be the satisfaction polynomial of $e$, so that $\#F=\sum_{x\in\{0,1\}^n}\prod_e\varphi_e(x)$.

> **Definition 5.1 (Clause-product interpolation scheme) [DL].** A scheme for $F$ consists of:
>
> **(a) Deformation.** For each occurrence $o=(e,v)$, a vector $P_o(T)\in\mathbb Q[T]^n$ of degree at most $R$ with $P_o(0)=e_v$, and
> $$H(T)=\sum_{x\in\{0,1\}^n}\prod_e\varphi_e\bigl(\langle P_{e,v}(T),x\rangle\bigr)_{v\in e}\qquad\text{([DL, (2)])}.$$
> The *own coefficient* of $o$ is $c_o(T)=(P_o(T))_v$. Since $\deg\varphi_e=|I_e|$, $\deg H\le MR$.
>
> **(b) Recovery.** A finite set $X\subset\mathbb Q\setminus\{0\}$, orders $q_t\ge 0$ for $t\in X$, and a linear functional $\Lambda$ on Hermite data $(f^{(j)}(t))_{t\in X,\,j\le q_t}$ such that $\Lambda(f)=f(0)$ for every $f\in\mathbb Q[T]$ of degree at most $MR$.
>
> **(c) Designs.** For each $t\in X$:
> - a partition of $[n]$ into groups $S_{t,1},S_{t,2},\dots$;
> - a group $\alpha_t(e)$ for each clause $e$;
> - a subspace $\mathcal R_{t,h}\subseteq\mathbb Q^{S_{t,h}}$ (the span of the measurement rows) for each group $h$,
>
> such that for every clause $e$, every $v\in I_e$ and every group $h\neq\alpha_t(e)$,
> $$P_{e,v}(t)\big|_{S_{t,h}}\in\mathcal R_{t,h}.\tag{3}$$
> Write $s_t=\max_h|S_{t,h}|$, $\rho_t=\sum_h\dim\mathcal R_{t,h}$, and $\mu_t$ for the number of measured groups ($\mathcal R_{t,h}\ne0$). The interaction graph $\mathcal I_t$ has the groups as vertices, and $w_t=\operatorname{tw}(\mathcal I_t)$.

**Cost model [DL, §5.2].** Each point $t$ is charged
$$\kappa(t)\ \ge\ \max\Bigl(2^{s_t},\,2^{w_t},\,\tbinom{n}{\min(q_t,\lfloor n/2\rfloor)}\Bigr)\qquad\text{(DP-type)}.$$
The preprint's *product-type* evaluation is additionally charged $\kappa(t)\ge 2^{\rho_t}$. A scheme costs at least $\max_{t\in X}\kappa(t)$. In [DL], Theorem 8.3 gives cost $2^{\Omega_\Delta(\sqrt n)}$ w.h.p. on random 3-CNF, in the DP-type (and hence also the product-type) model. In that theorem the scheme is applied to the preprint's bounded-occurrence formula $F^\star$ and to the leaves of its sparsification.

### 1.2 The question of §9

[DL, §9, "The gap"] reads, in part:

> "From above, a design on an expander with $\sqrt n$ groups, $O(\sqrt n)$ interaction width, Taylor order $O(\sqrt n)$ and $O(1)$ or polylogarithmically many rows per measured group would bring the scheme to $2^{\tilde O(\sqrt n)}$. From below, one must show that measurements cannot absorb the cut. A measured group of dimension one can absorb two foreign occurrences $(e,v),(e',v')$ with $v\ne v'$ at $t$ only if the $2\times2$ minor $c_{e,v}(t)c_{e',v'}(t)-p_{e,v,v'}(t)p_{e',v',v}(t)$ vanishes; that polynomial equals $1$ at $T=0$, so each such pair is absorbed at no more than $2R$ points. […] Turning this into a bound for groups with several measurement rows is open."

Such a design would have $s_t=\tilde O(\sqrt n)$, $q_t=\tilde O(\sqrt n)$ and $\rho_t\le(\#\text{groups})\cdot\mathrm{polylog}\,n=\tilde O(\sqrt n)$ at every point. In both cost models, $s_t=\tilde O(\sqrt n)$ is also necessary for cost $2^{\tilde O(\sqrt n)}$, because $2^{s_t}$ is charged. Our determinant replaces the $2\times2$ minor by an $n\times n$ minor, and this handles any number of rows.

### 1.3 Notation

Fix a point $t$, and write $S_h=S_{t,h}$, $a_h=|S_h|$ and $r_h=\dim\mathcal R_{t,h}$. Then $\sum_h a_h=n$ and $\sum_h r_h=\rho_t$. Write $(y)_+=\max(y,0)$.

- **Witness maps.** A witness map $\omega$ assigns to each variable $v$ with $d_v\ge1$ one of its occurrences $\omega(v)=(e,v)$. For a variable that never occurs we set $\omega(v)=\bot$ and $P_\bot(T):=e_v$ (constant).
- **Home group.** The *home group* of $v$ at $t$ under $\omega$ is $h^\omega_t(v)=\alpha_t(e)$ if $\omega(v)=(e,v)$. If $\omega(v)=\bot$, it is the group containing $v$.
- **Deficiency.** Put $I^\omega_{t,h}=\{v:h^\omega_t(v)=h\}$ and define the deficiency
$$\operatorname{def}^\omega_t=\sum_h\bigl(a_h-|I^\omega_{t,h}|\bigr)_+ .$$
- **Split variables.** A variable is *split at $t$* if $d_v\ge1$ and its occurrences lie in clauses assigned by $\alpha_t$ to at least two different groups. Let $\Sigma_t$ be the set of split variables and $\sigma'_t=|\Sigma_t|$. A non-split variable has a well-defined home group $h_t(v)$, independent of $\omega$. Put $J_{t,h}=\{v\notin\Sigma_t:h_t(v)=h\}$.
- **Separators.** A set $Y\subseteq[n]$ is a *$2/3$-separator* if every connected component of $\Pr(F)-Y$ has at most $2n/3$ vertices. $\operatorname{sep}(F)$ is the minimum size of a $2/3$-separator.
- **Weighted separators.** With $w(d)=2(d-1)/d^2$, let
$$W(F)=\min\Bigl\{\sum_{v\in Y}w(d_v)\ :\ Y\text{ a }2/3\text{-separator with }d_v\ge2\ \forall v\in Y\Bigr\},$$
with $W(F)=+\infty$ if no such $Y$ exists. Since $w$ is decreasing on $d\ge2$, $W(F)\ge\frac{2(b-1)}{b^2}\operatorname{sep}(F)$ when all $d_v\le b$.
- **Anticoncentration constant.** $c_0:=(1-e^{-2})/4=0.21616\ldots$. The source notes round this down to $0.216$.

---

## 2. Main results

**Fact R (Hermite counting).** In every scheme, $\sum_{t\in X}(q_t+1)\ge MR+1$.

**Lemma 1 (witness determinant).** For every scheme and every witness map $\omega$,
$$\sum_{t\in X}\bigl(\operatorname{def}^\omega_t-\rho_t\bigr)_+\ \le\ \sum_{t\in X}\dim\ker\Phi^\omega(t)\ \le\ nR,$$
where $\Phi^\omega(T)$ is the $n\times n$ matrix with row $v$ equal to $P_{\omega(v)}(T)$.

**Lemma 2 (expected deficiency).** Choose $\omega$ at random: independently for each variable, uniformly among its occurrences. For every $t$, with $s=s_t$ and $\mu_h=\mathbb E|I^\omega_{t,h}|$,
$$\mathbb E\operatorname{def}^\omega_t\ \ge\ \max\Bigl(\tfrac12\sum_h|\mu_h-a_h|,\ \ \frac{c_0}{2\sqrt{2s}}\sum_{v\in\Sigma_t}w(d_v)\Bigr).$$
The proof rests on the following inequality: for a sum $X$ of independent Bernoulli variables with variance $V$ and any $a\in\mathbb R$, $\mathbb E|X-a|\ge c_0\min(V,\sqrt V)$.

**Lemma 3 (split variables form a separator).** For every $t$, either $\Sigma_t$ is a $2/3$-separator (so $\sigma'_t\ge\operatorname{sep}(F)$), or
$$\tfrac12\sum_h|\mu_h-a_h|\ >\ \tfrac{2n}{3}-s_t .$$

**Theorem B (measurements absorb at most a $\sqrt{s}$ factor of the cut).** Let $F$ satisfy $M\ge n$ (for example, every variable occurs). Then every scheme for $F$ in the sense of Definition 5.1 has a point $t\in X$ with
$$\rho_t+q_t+1\ >\ \min\Bigl(\frac{c_0}{2\sqrt2}\cdot\frac{W(F)}{\sqrt{s_t}},\ \ \frac{2n}{3}-s_t\Bigr).$$
This holds for arbitrary clause-dependent deformations, partitions, assignments $\alpha_t$ and subspaces $\mathcal R_{t,h}$. If every variable occurs at most $b$ times, then
$$\rho_t+q_t+1\ >\ \min\Bigl(c_b\,\frac{\operatorname{sep}(F)}{\sqrt{s_t}},\ \frac{2n}{3}-s_t\Bigr),\qquad c_b=\frac{c_0\,(b-1)}{\sqrt2\,b^2},\quad c_3=0.03397\ldots$$

**Theorem A (variable-uniform deformations: no absorption).** Let $M\ge n$, and suppose $P_{e,v}=P_v$ for all clauses $e\ni v$; this is a global linear change of variables $x\mapsto\Phi(T)x$. Then some $t\in X$ satisfies
$$\rho_t+q_t+1>\min\bigl(\operatorname{sep}(F),\ \tfrac{2n}{3}-s_t\bigr).$$

**Lemma 4 (random 3-CNF; adapting [DL, Lemma 8.2]).** Let $F$ be a uniformly random 3-CNF with $m=\Delta n$ clauses, $\Delta>3$. With high probability, $\operatorname{sep}(F)\ge c_\Delta n$ and $W(F)\ge c'_\Delta n$ for constants $c_\Delta,c'_\Delta>0$. Moreover $M=3\Delta n\ge n$ deterministically.

**Corollary C.** Let $\mathcal F$ be one of the following families:

- **(i)** formulas with all $d_v\le b$ ($b$ fixed), $M\ge n$ and $\operatorname{sep}(F)\ge cn$ for a fixed $c>0$;
- **(ii)** uniformly random 3-CNF with $m=\Delta n$, $\Delta>3$ (w.h.p.).

Then, for $n$ large:

1. **(The §9 design does not exist.)** For every fixed $k$, no scheme on $F\in\mathcal F$ has $s_t,\ q_t,\ \rho_t\le\sqrt n\log^k n$ for all $t\in X$. More precisely, if $s_t\le\sqrt n\log^kn$ for all $t$, then some $t$ has $\rho_t+q_t=\Omega(n^{3/4}/\log^{k/2}n)$. This holds in particular for the design described in §9 ($\sqrt n$ groups of size $\tilde O(\sqrt n)$, Taylor order $O(\sqrt n)$, polylogarithmically many rows per group).
2. **(Product-type cost.)** Every scheme on $F\in\mathcal F$ costs at least $2^{\Omega(n^{2/3})}$ in the product-type model of [DL, §5.2]. The paper has $2^{\Omega(\sqrt n)}$.
3. **(Variable-uniform schemes.)** A variable-uniform scheme on $F\in\mathcal F$ costs $2^{\Omega(n)}$ in the product-type model.

For balanced designs ($g_t$ groups with $s_t=\Theta(n/g_t)$), Theorem B says: either $q_t=\Omega(\sqrt{n g_t})$, or the average number of rows per group, $\rho_t/g_t$, is $\Omega(\sqrt{s_t})$.

---

## 3. Proofs

### 3.1 Fact R

Suppose $\sum_t(q_t+1)\le MR$. Then $f(T)=\prod_{t\in X}(T-t)^{q_t+1}$ has degree at most $MR$, and all its Hermite data vanish. By linearity, $\Lambda(f)=\Lambda(0)=0$. But $f(0)=\prod_t(-t)^{q_t+1}\neq0$ because $0\notin X$, which contradicts (b). $\square$

### 3.2 Proof of Lemma 1

**The determinant is nonzero of degree at most $nR$.** Every entry of $\Phi(T)=\Phi^\omega(T)$ has degree at most $R$. Also $\Phi(0)=I$, because $P_o(0)=e_v$ and $P_\bot=e_v$. So $\det\Phi$ is a nonzero polynomial of degree at most $nR$.

**Rank bound at a point.** Fix $t\in X$ and a group $h$, and let $B_h=\Phi(t)[:,S_h]$ be the $n\times a_h$ column block. Take a row $v\notin I^\omega_{t,h}$, so $h^\omega_t(v)\neq h$. There are two cases:

- If $\omega(v)=(e,v)$, then $h\ne\alpha_t(e)$, and (3) gives $\Phi_v(t)|_{S_h}\in\mathcal R_{t,h}$.
- If $\omega(v)=\bot$, then $v\notin S_h$, so $\Phi_v(t)|_{S_h}=e_v|_{S_h}=0\in\mathcal R_{t,h}$.

So all rows of $B_h$ outside $I^\omega_{t,h}$ lie in a space of dimension $r_h$. Hence $\operatorname{rank}B_h\le\min(a_h,\ r_h+|I^\omega_{t,h}|)$. The blocks partition the columns, so $\operatorname{rank}\Phi(t)\le\sum_h\operatorname{rank}B_h$. This gives
$$\dim\ker\Phi(t)\ \ge\ \sum_h\bigl(a_h-r_h-|I^\omega_{t,h}|\bigr)_+\ \ge\ \sum_h\Bigl(\bigl(a_h-|I^\omega_{t,h}|\bigr)_+-r_h\Bigr)=\operatorname{def}^\omega_t-\rho_t,$$
using $(y-r)_+\ge y_+-r$ for $r\ge0$. Since the kernel dimension is nonnegative, it is at least $(\operatorname{def}^\omega_t-\rho_t)_+$.

**Order of vanishing.** Let $k=\dim\ker\Phi(t)$. Choose $U\in GL_n(\mathbb Q)$ whose first $k$ columns span $\ker\Phi(t)$. The first $k$ columns of $\Phi(T)U$ vanish at $T=t$, so each is divisible by $T-t$. Hence $(T-t)^k$ divides $\det(\Phi U)=\det U\cdot\det\Phi$. The points of $X$ are distinct, so $\prod_{t\in X}(T-t)^{\dim\ker\Phi(t)}$ divides the nonzero polynomial $\det\Phi$, and therefore $\sum_t\dim\ker\Phi(t)\le\deg\det\Phi\le nR$. $\square$

*Remark.* The bound is tight. For the "cut" designs in §4, $\sum_t\dim\ker\Phi(t)=nR$ exactly.

### 3.3 Proof of Lemma 2

Write $X_h=|I^\omega_{t,h}|$. Then $\sum_hX_h=n=\sum_ha_h$, so $\sum_h(a_h-X_h)=0$, and therefore $\operatorname{def}^\omega_t=\sum_h(a_h-X_h)_+=\frac12\sum_h|X_h-a_h|$. Taking expectations, $\mathbb E\operatorname{def}_t=\frac12\sum_h\mathbb E|X_h-a_h|$.

**Structure of $X_h$.** $X_h=\sum_v\mathbf 1[h^\omega_t(v)=h]$ is a sum of independent Bernoulli variables with parameters $\pi_{v,h}$, where $\pi_{v,h}$ is the fraction of $v$'s occurrences in clauses assigned to $h$. For a variable that never occurs, $\pi_{v,h}=1$ for its own group. Hence $\mu_h=\sum_v\pi_{v,h}$ and $V_h:=\operatorname{Var}X_h=\sum_v\pi_{v,h}(1-\pi_{v,h})\le\mu_h$. Jensen gives $\mathbb E|X_h-a_h|\ge|\mu_h-a_h|$, which is the first bound.

**Anticoncentration.** Let $X$ be a sum of independent Bernoulli($p_v$) variables with variance $V$, and let $a\in\mathbb R$. Let $X'$ be an independent copy. By the triangle inequality, $\mathbb E|X-a|\ge\frac12\mathbb E|Z|$ with $Z=X-X'=\sum_v\xi_v$. The $\xi_v$ are independent and symmetric on $\{-1,0,1\}$, with $\Pr(\xi_v=\pm1)=p_v(1-p_v)$ each, so $\mathbb EZ^2=2V$.

- **(i) Bound for all $V$.** Condition on the set $K=\{v:\xi_v\ne0\}$. If $K\ne\emptyset$, $Z$ is a sum of $|K|$ independent signs, so $\Pr(Z=0\mid K)\le\frac12$. Since $\Pr(K\neq\emptyset)=1-\prod_v(1-2p_v(1-p_v))\ge1-e^{-2V}$, and $|Z|\ge1$ whenever $Z\ne0$, we get $\mathbb E|Z|\ge\Pr(Z\ne0)\ge\frac12(1-e^{-2V})$.
- **(ii) Bound for $V\ge1$.** Write $\tau_v=\mathbb E\xi_v^2=\mathbb E\xi_v^4$. Then $\mathbb EZ^4=\sum_v\tau_v+3\sum_{u\ne v}\tau_u\tau_v\le 2V+12V^2$. Hölder's inequality ($\mathbb EZ^2\le(\mathbb E|Z|)^{2/3}(\mathbb EZ^4)^{1/3}$) gives $\mathbb E|Z|\ge(2V)^{3/2}/(2V+12V^2)^{1/2}\ge\sqrt{4/7}\,\sqrt V\ge0.7559\sqrt V$ for $V\ge1$.

For $V<1$, (i) and the monotonicity of $(1-e^{-2V})/V$ give $\mathbb E|X-a|\ge\frac14(1-e^{-2V})\ge\frac{1-e^{-2}}4V=c_0V$. For $V\ge1$, (ii) gives $\mathbb E|X-a|\ge0.378\sqrt V\ge c_0\sqrt V$. In both cases, $\mathbb E|X-a|\ge c_0\min(V,\sqrt V)$.

**From anticoncentration to a bound in terms of $V_h/\sqrt{2s}$.** We show $\mathbb E|X_h-a_h|\ge c_0V_h/\sqrt{2s}$ in two cases:

- **If $\mu_h\le2s$:** then $V_h\le2s$, and since $s\ge1$, $\min(V_h,\sqrt{V_h})\ge V_h/\sqrt{2s}$.
- **If $\mu_h>2s$:** then $a_h\le s$, so $\mathbb E|X_h-a_h|\ge\mu_h-s\ge\mu_h/2\ge V_h/2\ge c_0V_h/\sqrt{2s}$.

**Summing over groups.** $\sum_hV_h=\sum_v\bigl(1-\sum_h\pi_{v,h}^2\bigr)$. For $v\in\Sigma_t$, the $\pi_{v,h}$ are multiples of $1/d_v$ with at least two of them nonzero, so $\sum_h\pi_{v,h}^2\le\bigl((d_v-1)^2+1\bigr)/d_v^2$ and $1-\sum_h\pi_{v,h}^2\ge w(d_v)$. All other terms are nonnegative, so $\mathbb E\operatorname{def}_t\ge\frac{c_0}{2\sqrt{2s}}\sum_{v\in\Sigma_t}w(d_v)$. $\square$

### 3.4 Proof of Lemma 3

**Components of $\Pr(F)-\Sigma_t$ stay inside one $J_{t,h}$.** Take an edge $uv$ of $\Pr(F)$ with $u,v\notin\Sigma_t$. Both variables lie in a common clause $e$, so both occur. Neither is split, so $h_t(u)=\alpha_t(e)=h_t(v)$. Variables that never occur are isolated vertices of $\Pr(F)$. Hence every connected component of $\Pr(F)-\Sigma_t$ lies inside a single $J_{t,h}$.

**Case 1: $|J_{t,h}|\le2n/3$ for every $h$.** Then $\Sigma_t$ is a $2/3$-separator, so $\sigma'_t\ge\operatorname{sep}(F)$.

**Case 2: $|J_{t,h_0}|>2n/3$ for some $h_0$.** Every $v\in J_{t,h_0}$ has $\pi_{v,h_0}=1$, so $\mu_{h_0}\ge|J_{t,h_0}|>2n/3$, while $a_{h_0}\le s_t$. Since $\sum_h\mu_h=n=\sum_ha_h$,
$$\tfrac12\sum_h|\mu_h-a_h|=\sum_h(\mu_h-a_h)_+\ \ge\ \mu_{h_0}-a_{h_0}\ >\ \tfrac{2n}3-s_t.\qquad\square$$

This sharpens the bound $(2n/3-s_t)/2$ stated in the source notes by a factor of 2, as observed by the adversarial review (see "Issues found").

### 3.5 Proof of Theorem B

**Averaging Lemma 1.** Average Lemma 1 over the random $\omega$ of Lemma 2. Since $(\cdot)_+$ is convex, Jensen's inequality gives
$$\sum_{t\in X}\bigl(\mathbb E\operatorname{def}_t-\rho_t\bigr)_+\ \le\ \mathbb E\sum_t(\operatorname{def}_t-\rho_t)_+\ \le\ nR\ \le\ MR\ <\ \sum_t(q_t+1).$$
The last two steps use $M\ge n$ and Fact R. So some $t$ has $(\mathbb E\operatorname{def}_t-\rho_t)_+<q_t+1$, that is, $\mathbb E\operatorname{def}_t<\rho_t+q_t+1$.

**Lower bound at that point.** Apply Lemma 3:

- **In Case 1,** $\Sigma_t$ is a $2/3$-separator consisting of variables with $d_v\ge2$. So $\sum_{v\in\Sigma_t}w(d_v)\ge W(F)$, and Lemma 2 gives $\mathbb E\operatorname{def}_t\ge\frac{c_0}{2\sqrt2}W(F)/\sqrt{s_t}$.
- **In Case 2,** Lemma 2 (first bound) and Lemma 3 give $\mathbb E\operatorname{def}_t>2n/3-s_t$.

**Bounded occurrences.** If $d_v\le b$ for all $v$, then $W(F)\ge\frac{2(b-1)}{b^2}\operatorname{sep}(F)$ and $\frac{c_0}{2\sqrt2}\cdot\frac{2(b-1)}{b^2}=c_b$. $\square$

### 3.6 Proof of Theorem A

Use the same determinant with $\Phi_v=P_v$, and $\Phi_v=e_v$ for variables that never occur. Row $v$ is constrained by (3) in block $h$ whenever some clause of $v$ is assigned to a group other than $h$, or, for a variable that never occurs, whenever $h$ is not its own group. So the unconstrained rows of block $h$ are exactly $J_{t,h}$. As in §3.2,
$$\dim\ker\Phi(t)\ \ge\ \sum_h(a_h-|J_{t,h}|)_+-\rho_t\ =\ n-\sum_h\min(a_h,|J_{t,h}|)-\rho_t .$$
Lemma 1's divisibility argument and Fact R then give a point $t$ with $n-\sum_h\min(a_h,|J_{t,h}|)<\rho_t+q_t+1$. Now apply Lemma 3's case split:

- **In Case 1,** $\sum_h\min(a_h,|J_{t,h}|)\le\sum_h|J_{t,h}|=n-\sigma'_t\le n-\operatorname{sep}(F)$.
- **In Case 2,** $\sum_h\min(a_h,|J_{t,h}|)\le a_{h_0}+\sum_{h\ne h_0}|J_{t,h}|<s_t+n/3$.

$\square$

### 3.7 Proof of Lemma 4 (the hypothesis $\operatorname{sep}(F)\ge cn$ for random 3-CNF)

We use the two facts that [DL] prove in its Appendix A, for $\Delta>3$, w.h.p.:

- **(A.1)** Every set $S$ of variables with $n/3\le|S|\le2n/3$ shares at least $m/3$ clauses with its complement, that is, at least $m/3$ clauses meet both $S$ and $[n]\setminus S$.
- **(A.2)** For every $j\ge1$, every set of $j$ variables occurs in at most $C_\Delta j\log(en/j)$ clauses.

**Bound on $\operatorname{sep}(F)$.** Let $Y$ be a $2/3$-separator with $|Y|\le n/9$. Then $\Pr(F)-Y$ has at least $8n/9$ vertices, in components of size at most $2n/3$. Choose $A$ as follows. If some component has at least $n/3$ vertices, let $A$ be that component. Otherwise, add components one at a time until the union reaches $n/3$. In both cases $n/3\le|A|\le2n/3$.

A clause meeting both $A$ and $[n]\setminus A$ must contain a variable of $Y$, since otherwise it would join $A$ to another component of $\Pr(F)-Y$. By (A.1) and (A.2),
$$m/3\le C_\Delta|Y|\log(en/|Y|).$$
Write $y=|Y|/n$. With $m=\Delta n$ this reads $\Delta/3\le C_\Delta\,y\log(e/y)$. The right side tends to $0$ as $y\to0$, so $y\ge c^*_\Delta>0$. Hence $\operatorname{sep}(F)\ge c_\Delta n$ with $c_\Delta=\min(1/9,c^*_\Delta)$. This is the argument of [DL, Lemma 8.2], with the treewidth separator (components of size at most $n/2$) replaced by a $2/3$-separator.

**Bound on $W(F)$.** $\sum_vd_v=M=3\Delta n$, so at most $3\Delta n/D$ variables have $d_v>D$. Take $D=\lceil6\Delta/c_\Delta\rceil$. Any $2/3$-separator $Y$ with all $d_v\ge2$ has $|Y|\ge c_\Delta n$, so it contains at least $c_\Delta n/2$ variables with $2\le d_v\le D$, each of weight at least $w(D)$. Hence $W(F)\ge\frac{c_\Delta}{2}w(D)\,n=:c'_\Delta n$. $\square$

*Note.* Random 3-CNF has, w.h.p., about $e^{-3\Delta}n$ variables that never occur. This is why Lemmas 1–3 are stated with the convention $P_\bot=e_v$ and the hypothesis $M\ge n$, rather than "every variable occurs".

### 3.8 Proof of Corollary C

Under (i) let $\gamma=c_bc$, and under (ii) let $\gamma=\frac{c_0}{2\sqrt2}c'_\Delta$. Theorem B, together with Lemma 4 for (ii), gives a point $t$ with
$$\rho_t+q_t+1>\min\bigl(\gamma n/\sqrt{s_t},\ 2n/3-s_t\bigr).$$

**1. The §9 design.** If $s_t\le\sqrt n\log^kn$, the right side is at least $\min(\gamma n^{3/4}/\log^{k/2}n,\ 2n/3-o(n))=\gamma n^{3/4}/\log^{k/2}n$ for large $n$. This exceeds $2\sqrt n\log^kn+1$.

**2. Product-type cost.** At this $t$, using $\binom nj\ge(n/j)^j\ge2^j$ for $j\le n/2$,
$$\kappa(t)\ge2^{\max(s_t,\ \rho_t,\ \min(q_t,\lfloor n/2\rfloor))}.$$
If $s_t\ge n^{2/3}$ we are done. Otherwise $2n/3-s_t\ge n/2$ for large $n$, and $\gamma n/\sqrt{s_t}\ge\gamma n^{2/3}$, so $\max(\rho_t,q_t)\ge(\gamma n^{2/3}-1)/2$. Hence $\max(\rho_t,\min(q_t,\lfloor n/2\rfloor))=\Omega(n^{2/3})$. (The balance point is $\max(s,n/\sqrt s)\ge n^{2/3}$.)

**3. Variable-uniform schemes.** Same argument with Theorem A. If $s_t<n/6$, then $\rho_t+q_t+1>\min(\operatorname{sep}(F),n/2)=\Omega(n)$. Under (ii), $\operatorname{sep}(F)\ge c_\Delta n$ by Lemma 4.

**Balanced reformulation.** If $s_t=\Theta(n/g_t)$, then $n/\sqrt{s_t}=\Theta(\sqrt{ng_t})$ and $\sqrt{ng_t}/g_t=\Theta(\sqrt{s_t})$. $\square$

---

## 4. Computational checks

All checks are in `code/gen6/g6a/` (by the original agent) and `code/gen6/g6r/` (by the reviewer) of the project repository. They are sanity checks of the algebra and the inequalities; they are not evidence of the asymptotic statements.

1. **Lemma 1 and Theorem B inequality chain** (`run_checks.py`, output `theorem_checks.txt`).
   - **Instances.** Random 3-regular formulas (every variable in exactly 3 clauses) with:
     - $n=24$: $|X|=8$, $R=2$, uniform $q=18$ (so $\sum(q+1)=152>MR=144$);
     - $n=30$: $|X|=10$, $R=2$, $q=18$ ($190>180$).

     Groups of size $s=6$ (4 or 5 groups), random partitions per point, majority assignment $\alpha_t$.
   - **Design families.** Four valid families: variable-uniform cut+measure; clause-dependent cut+measure; "global-sum" (one all-ones row per group plus leftovers); random dense deformation with minimal $\mathcal R_{t,h}$.
   - **What was checked.**
     - Condition (3) exactly, over $\mathbb Q$ (Fractions).
     - For 4 random witness maps per design: $\dim\ker\Phi^\omega(t)\ge\operatorname{def}^\omega_t-\rho_t$ at every $t$, and $\sum_t\dim\ker\Phi^\omega(t)\le nR$. Equality $nR$ (48 and 60) is attained by the cut designs, so Lemma 1 is tight.
     - $\mathbb E\operatorname{def}_t$, computed from the exact distribution of each $X_h$ (in floating point), is at least the Lemma 2 bound $\max(c_3\sigma'_t/\sqrt{s_t},\ \frac12\sum_h|\mu_h-a_h|)$.
   - **Caveat.** These naive designs have $\rho_t\approx(0.4$–$1.0)\,n$, far from the regime of Theorem B, so they are not a tightness test of the $\sqrt s$ loss.
2. **End-to-end recovery** (`e2e_checks.txt`).
   - **Instance.** $n=12$, 3 groups of 4, $R=1$, $|X|=6$, $q=6$ ($42>MR=36$), two design families.
   - **Exact computations.** $H(T)$ was computed exactly modulo $2^{61}-1$; $H(0)=\#F=977$, and Hermite recovery from the data at $X$ returns 977.
   - **Product-type evaluation.** At every $t$, the product-type design evaluation (enumerating tuples of measurement readings, 64 to 3072 tuples, and counting inside groups) equals the brute-force $H(t)$ exactly in $\mathbb Q$.
3. **Anticoncentration** (`anticonc_check.py`). 4000 random sums of up to 40 Bernoulli variables, minimising over integer $a$ (a median, where the minimum of $\mathbb E|X-a|$ is attained). The worst ratio $\mathbb E|X-a|/\min(V,\sqrt V)$ is $0.74$, against the proved $0.216$.
4. **Independent check by the reviewer** (`g6r/lemma1_check.py`, sympy).
   - **Instances.** Random schemes satisfying (3) by construction: $n=9$, 7 clauses, $R=2$, 4 points, 3 groups, one measurement row per group. The deformations are solved from the linear constraints, with random values for the free parameters.
   - **What was checked,** for 6 seeds times 6 witness maps: $P_o(0)=e_v$; $\det\Phi(0)=1$ and $\deg\det\Phi\le nR$; $\dim\ker\Phi(t)\ge\operatorname{def}_t-\rho_t$; the multiplicity of the root $t$ of $\det\Phi$ is at least $\dim\ker\Phi(t)$; and $\sum_t(\operatorname{def}_t-\rho_t)_+\le nR$.

---

## 5. Scope and open problems

- **What is proved.** Theorems A and B use only (a), (b), (c) of Definition 5.1, the degree bound $\deg H\le MR$ and $M\ge n$. They hold for arbitrary clause-dependent deformations, arbitrary partitions and assignments, and arbitrary measurement subspaces, with no cost model. The cost model enters only in Corollary C(2,3).
- **DP-type model: open.** The DP-type model charges $2^{w_t}$, not $2^{\rho_t}$. Theorem B forces many rows, but rows spread over many groups need not enter a single DP bag ($w_t\le\mu_t\le\rho_t$ only bounds $w_t$ from above). In the DP-type model, the gap between $2^{\Omega(\sqrt n)}$ [DL, Thm 8.3] and the preprint's $2^{O(\rho(n))}$ remains open. Here $\rho(n)$ is the preprint's exponent function, unrelated to $\rho_t$. Corollary C(1) does exclude the specific §9 target design in both models, since that design is described by its parameters.
- **Group size.** The refutation of the §9 design assumes $s_t=\tilde O(\sqrt n)$. This assumption is implied by the target cost $2^{\tilde O(\sqrt n)}$ in both cost models, because $2^{s_t}$ is charged.
- **Hypothesis $\operatorname{sep}(F)\ge cn$.** It is proved here for uniformly random 3-CNF ($\Delta>3$) by Lemma 4. For random 3-regular variable–clause incidences it is assumed, not proved. We have not checked whether it holds for the preprint's bounded-occurrence formula $F^\star$ (Tovey cycles) or for the leaves of its sparsification, which [DL, Thm 8.3] treats via treewidth. Our random-3-CNF corollary is for $F$ itself.
- **Tightness of the $\sqrt s$ loss: open.** The loss comes only from fluctuations of a random witness. It is open whether some design achieves $\rho_t\approx n/\sqrt{s_t}$, which would give product-type cost $2^{\tilde\Theta(n^{2/3})}$. By Theorem A, such a design must be clause-dependent.
- **Not covered.** The preprint's parity form and its tensor contraction escape Definition 5.1: [DL, §5.1] explains that they break the degree bound or the exactness of recovery in (b). Our results say nothing about them.
- **Relation to [DL, §9].** The $2\times2$ minor of §9 is a $2\times2$ minor of the witness matrix $\Phi$. Lemma 1 uses the full determinant, which handles any number of rows per group at once.

---

## 6. How this result was obtained

This note, the results and the proofs were produced by AI agents (Anthropic's Claude) in an automated research workflow:

- One agent found the result and wrote the original notes, proofs and checks.
- A separate adversarial AI review checked each lemma line by line, wrote independent tests, and found no prior publication of the result. It confirmed the result with small fixes, which are incorporated here: the convention $P_\bot=e_v$ for variables that never occur; the sharper second case of Lemma 3; and the explicit hypothesis $\operatorname{sep}(F)\ge cn$, with its proof for random 3-CNF adapted from [DL, Lemma 8.2].
- This note was written by a third agent from those materials. It has **not yet been reviewed by human experts**, and errors may remain.

The quoted definition and question are taken from the arXiv HTML version of [DL]. The cited facts (A.1), (A.2) are from its Appendix A; their proofs were read but are not reproduced here.

**Reference.** [DL] F. Dias, Y. Lima, *Limits and Constants of OpenAI's Claimed ETH Refutation*, arXiv:2611.00517v1 [cs.CC], 31 Oct 2026.

---

## Issues found (relative to the source notes and the journal)

1. **Lemma 3 / Theorem B, second term.** The notes and the journal state $(2n/3-s_t)/2$. The proof gives $>2n/3-s_t$, because $\sum_h\mu_h=\sum_ha_h$ (this is the reviewer's fix). The weaker form remains true. This note states the sharper form, consistent with Theorem A.
2. **Constant $c_3$.** The journal writes the constant as "0.034". The proved value is $c_3=c_0(b-1)/(\sqrt2b^2)$ with $c_0=(1-e^{-2})/4$, so $c_3=0.03397\ldots$, slightly *below* $0.034$; with the notes' rounded $c_0=0.216$ it is $0.03394$. Writing "$\approx0.034$" is fine; writing "$0.034$" as a proved constant is very slightly too strong.
3. **Source of $c_0$.** The notes describe $0.216$ as "numerically checked". In fact it follows rigorously from step (i) at $V=1$: $(1-e^{-2})/4=0.2162$. The numerical check is only a sanity check.
4. **Hypothesis "every variable occurs".** The notes' Theorem B assumes that every variable occurs. Random 3-CNF violates this w.h.p., so the notes' random-3-CNF corollary needs the $e_v$ convention and the hypothesis $M\ge n$ used here. Random 3-CNF also has unbounded degrees, so it needs the weighted quantity $W(F)$ (the notes' parenthetical remark), not $c_b\operatorname{sep}(F)$.
5. **"Equivalently" (notes, Corollary 1).** The reformulation "average rows per group $\ge c\sqrt{s_t}$ unless $q_t\ge c\sqrt{ng_t}$" holds for balanced designs ($s_t=\Theta(n/g_t)$). For unbalanced designs only the $s_t$ form is proved.
6. **"Exact" $\mathbb E\operatorname{def}$ in the checks.** `run_checks.py` computes $\mathbb E\operatorname{def}_t$ from the exact distribution, but in floating point, not in Fractions.
7. **Notes, Corollary 3 ("the preprint's deformation must be clause-dependent").** This follows only for formulas to which both the preprint's subexponential bound and $\operatorname{sep}\ge cn$ apply. Whether $\operatorname{sep}(F^\star)\ge cn^\star$ holds is unverified, so this note does not state the remark unconditionally.

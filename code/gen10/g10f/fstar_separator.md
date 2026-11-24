# Linear separators for the bounded-occurrence formula $F^\star$ (agent G10f)

Status: complete paper proof, not machine-checked, not reviewed by a second agent.
Constants were recomputed numerically and the deterministic lemma was sanity-tested on
200 small random instances (`check.py`, same directory). Everything about the preprint's
construction is taken from Dias--Lima's description only ("variable copies on Tovey cycles",
"copies tied by Tovey's implication cycle $(\neg x_1\vee x_2)\wedge\dots\wedge(\neg x_s\vee x_1)$").

## 0. Where the naive reduction fails

The suggested route "a $2/3$-separator $Y^\star$ of $\Gr(F^\star)$ projects to a $2/3$-separator
$Y$ of $\Gr(F)$, then apply Lemma~\ref{lem:random}" does **not** work as stated.
Contracting copy-cycles does give a set $Y$ with $|Y|\le|Y^\star|$, but the components of
$\Gr(F)-Y$ are only balanced for the *occurrence weight* $d(K)=\sum_{v\in K}d_v\le 2n^\star/3$,
not for the number of vertices: a component may contain more than $2n/3$ variables (the
low-degree ones) and still carry fewer than $2/3$ of the occurrences. So $Y$ need not be a
$2/3$-separator of $F$, property (A1) (sets with $n/3\le|S|\le2n/3$) does not apply, and
$\sepw(F)$ of the paper is the wrong notion (it weights separator vertices by $\theta(d_v)$,
and still measures components by cardinality). What is needed is an expansion property for
cuts that are balanced *in clauses*. Lemma~\ref{lem:P} below provides it by a first-moment
argument; it replaces (A1), and holds for $\Delta\ge2$.

## 1. LaTeX-ready text

```latex
\subsection{The bounded-occurrence formula $F^\star$}\label{sec:fstar}

Throughout this subsection $F$ is a 3-CNF on the variable set $V=[n]$ with $m$ clauses, each
on three distinct variables; $d_v$ is the number of clauses containing $v$, so
$\sum_vd_v=3m$, and $d(S)=\sum_{v\in S}d_v$. A clause \emph{lies inside} $S\subseteq V$ if
its three variables are in $S$, and \emph{crosses} $S$ if it has a variable in $S$ and a
variable outside $S$.

\begin{definition}[Copy expansion]\label{def:copy}
A CNF $F^\star$ on a variable set $V^\star$ is a \emph{copy expansion} of $F$ if there is a
map $\pi:V^\star\to V$ such that
\begin{itemize}
\item[(E1)] for every $v\in V$ the fibre $\pi^{-1}(v)$ is empty or induces a connected
 subgraph of $\Gr(F^\star)$;
\item[(E2)] for every clause $e$ of $F$ there is a clause $e^\star$ of $F^\star$ such that
 $\pi$ maps the variable set of $e^\star$ bijectively onto $I_e$;
\item[(E3)] $|\pi^{-1}(v)|=d_v$ if $d_v\ge1$, and $|\pi^{-1}(v)|\le1$ if $d_v=0$.
\end{itemize}
Let $n_0$ be the number of $v$ with $d_v=0$ and $\pi^{-1}(v)\neq\emptyset$. Then
$n^\star=|V^\star|=3m+n_0$.
\end{definition}

The formula $F^\star$ of~\cite[Theorem~8.3(a)]{DL} is of this form: every variable $v$ with
$d_v\ge1$ is replaced by $d_v$ copies $v_1,\dots,v_{d_v}$, one per occurrence, and the copies
are tied by Tovey's implication cycle
$(\neg v_1\vee v_2)\wedge\dots\wedge(\neg v_{d_v}\vee v_1)$~\cite{Tov84}. The fibre of $v$ is
a cycle (an edge if $d_v=2$, a vertex if $d_v=1$), so (E1) holds; (E2) and (E3) hold by
construction; a variable with $d_v=0$ is kept or dropped. Every variable of this $F^\star$
occurs at most three times (once in a clause of $F$ and at most twice in its cycle). The
order of the copies on the cycles plays no role below.

\begin{lemma}[Transfer]\label{lem:T}
Let $F^\star$ be a copy expansion of $F$ and let $Y^\star$ be a $2/3$-separator of $F^\star$
with $|Y^\star|\le n^\star/9$. Put $Y=\pi(Y^\star)$ and $L=d(Y)$. Then $|Y|\le|Y^\star|$ and
there is a partition $V=A\cup B\cup Y$ such that
\begin{enumerate}
\item[\textup{(1)}] every clause of $F$ lies inside $A$, or lies inside $B$, or contains a
 variable of $Y$; in particular at most $L$ clauses cross $A$;
\item[\textup{(2)}] at least $\frac m3-L-\frac{n_0}3$ clauses lie inside $A$, and at least
 $\frac{2m}9-L-\frac{n_0}3$ clauses lie inside $B$.
\end{enumerate}
\end{lemma}

\begin{proof}
The components of $\Gr(F^\star)-Y^\star$ have at most $2n^\star/3$ vertices each and at
least $8n^\star/9$ vertices in total. As in the proof of Lemma~\ref{lem:random}, there is a
union $A^\star$ of components with $n^\star/3\le|A^\star|\le2n^\star/3$ (the largest
component if it has at least $n^\star/3$ vertices; otherwise add components one at a time
until the union reaches $n^\star/3$). Let $B^\star=V^\star\setminus(Y^\star\cup A^\star)$;
then $|B^\star|\ge n^\star-\frac23n^\star-\frac19n^\star=\frac29n^\star$, and no edge of
$\Gr(F^\star)$ joins $A^\star$ and $B^\star$.

Let $v\notin Y$. Its fibre is disjoint from $Y^\star$ and connected by (E1), so it lies
inside $A^\star$ or inside $B^\star$. Let $A$ be the set of $v\notin Y$ whose fibre is
nonempty and inside $A^\star$, and $B=V\setminus(Y\cup A)$.

(1) Let $e$ be a clause of $F$ with $I_e\cap Y=\emptyset$, and $e^\star$ as in (E2). The
variables of $e^\star$ are pairwise adjacent in $\Gr(F^\star)$ and none is in $Y^\star$, so
all are in $A^\star$ or all are in $B^\star$; hence $I_e\subseteq A$ or $I_e\subseteq B$. A
clause that crosses $A$ lies neither inside $A$ nor inside $B$, so it contains a variable of
$Y$, and at most $d(Y)=L$ clauses contain a variable of $Y$.

(2) If $x\in A^\star$ and $\pi(x)\notin Y$, the fibre of $\pi(x)$ meets $A^\star$, so
$\pi(x)\in A$. Hence $A^\star\subseteq\pi^{-1}(A)\cup\pi^{-1}(Y)$, and by (E3)
\[
\frac{n^\star}3\ \le\ |A^\star|\ \le\ d(A)+d(Y)+n_0,\qquad\text{so}\qquad
d(A)\ \ge\ m-L-n_0 ,
\]
using $n^\star\ge3m$. Let $m_A$ be the number of clauses inside $A$. An occurrence of a
variable of $A$ lies in a clause inside $A$ (three such occurrences per clause) or, by (1),
in a clause containing a variable of $Y$ (at most two such occurrences per clause, at most
$L$ clauses). So $d(A)\le3m_A+2L$ and $m_A\ge\frac m3-L-\frac{n_0}3$. In the same way
$B^\star\subseteq\pi^{-1}(B)\cup\pi^{-1}(Y)$ gives
$d(B)\ge\frac29n^\star-L-n_0\ge\frac{2m}3-L-n_0$ and $d(B)\le3m_B+2L$, so
$m_B\ge\frac{2m}9-L-\frac{n_0}3$.
\end{proof}

From now on $F$ is the random 3-CNF of Lemma~\ref{lem:random} with $m=\Delta n$, and
$H(x)=-x\ln x-(1-x)\ln(1-x)$.

\begin{lemma}[No clause-balanced sparse cut]\label{lem:P}
Let $\Delta\ge2$. With high probability there is no set $S\subseteq[n]$ such that at least
$m/6$ clauses lie inside $S$, at least $m/6$ clauses lie inside $[n]\setminus S$, and at most
$m/20$ clauses cross $S$.
\end{lemma}

\begin{proof}
Put $\beta=1/6$, $\varepsilon=1/20$. Fix $S$ and $\alpha=|S|/n$. A clause lies inside $S$
with probability $p\le\alpha^3$ and inside the complement with probability
$p'\le(1-\alpha)^3$ (as in Lemma~\ref{lem:A1}). If the event holds for $S$, let $I$ and $I'$
be the sets of indices of the clauses inside $S$ and inside its complement. They are
disjoint, $|I|,|I'|\ge\beta m$ and $|I|+|I'|\ge(1-\varepsilon)m$. The number of pairs
$(I,I')$ with these properties is at most
$\sum_{j\le\varepsilon m}\binom mj\,2^{m}\le e^{H(\varepsilon)m}2^m$ (choose the complement
of $I\cup I'$, then split the rest). For a fixed pair with $|I|=i$, $|I'|=i'$, the clauses
are independent, so the probability that all clauses of $I$ lie inside $S$ and all clauses of
$I'$ inside the complement is
\[
p^{i}p'^{\,i'}\ \le\ \bigl(\alpha^{i}(1-\alpha)^{i'}\bigr)^3\ \le\
\exp\bigl(-3(i+i')H(\tfrac{i}{i+i'})\bigr)\ \le\ \exp\bigl(-3(1-\varepsilon)H(\beta)m\bigr).
\]
The second inequality is the maximum of $\alpha\mapsto\alpha^i(1-\alpha)^{i'}$ on $[0,1]$,
attained at $\alpha=i/(i+i')$. The third holds because $i,i'\ge\beta m$ and $i+i'\le m$ give
$i/(i+i')\in[\beta,1-\beta]$, where $H\ge H(\beta)$ by concavity and symmetry. By the union
bound over the pairs $(I,I')$ and over the $2^n$ sets $S$, the probability that some $S$
exists is at most
\[
\exp\bigl(n\ln2-\kappa\,\Delta n\bigr),\qquad
\kappa=3(1-\varepsilon)H(\beta)-\ln2-H(\varepsilon)=1.28410-0.69315-0.19852=0.39243\ldots
\]
Since $\ln2/\kappa=1.766\ldots<2\le\Delta$, this tends to $0$.
\end{proof}

We also use two facts.
\begin{itemize}
\item[(A2$'$)] With high probability, for every $1\le j\le n$ and every set $S$ of $j$
 variables, $d(S)\le C_\Delta\,j\log(en/j)$ with $C_\Delta=3e^2\Delta$. This is what the
 proof of Lemma~\ref{lem:A2} shows: it bounds $L_S=d(S)$, and uses only
 $C_\Delta\ge3$, that is $\Delta\ge e^{-2}$.
\item[(N)] Let $N_0$ be the number of variables with $d_v=0$. If $\Delta\ge2$, then
 $N_0\le n/60$ with high probability. Indeed $\PP(d_v=0)=(1-3/n)^m\le e^{-3\Delta}\le e^{-6}$,
 so $\E N_0\le e^{-6}n<n/120$; and for $u\ne v$,
 $\PP(d_u=d_v=0)=\bigl(\tfrac{(n-3)(n-4)}{n(n-1)}\bigr)^m\le\bigl(\tfrac{n-3}{n}\bigr)^{2m}
 =\PP(d_u=0)\PP(d_v=0)$ because $n(n-4)<(n-1)(n-3)$. So $\mathrm{Var}(N_0)\le\E N_0$ and, by
 Chebyshev's inequality, $\PP(N_0\ge n/60)\le\E N_0/(n/120)^2\le14400\,e^{-6}/n\to0$.
\end{itemize}

\begin{theorem}[Linear separators for $F^\star$]\label{thm:fstar}
Fix $\Delta\ge2$ and let $F$ be a random 3-CNF as in Lemma~\ref{lem:random}. Let
$\eta^\star=0.000241\ldots$ be the solution in $(0,1)$ of $y\log(e/y)=1/(60e^2)$. With
probability tending to $1$ as $n\to\infty$, the following holds simultaneously for every copy
expansion $F^\star$ of $F$ (in particular for every choice of the orders on the Tovey
cycles): every $2/3$-separator $Y^\star$ of $F^\star$ satisfies $|\pi(Y^\star)|\ge\eta^\star n$
or $|Y^\star|>n^\star/9$. Consequently
\[
\sep(F^\star)\ \ge\ \eta^\star n\ \ge\ \frac{\eta^\star}{3\Delta+1}\,n^\star .
\]
\end{theorem}

\begin{proof}
Assume the events of Lemma~\ref{lem:P}, (A2$'$) and (N); they depend on $F$ only. Let
$Y^\star$ be a $2/3$-separator of $F^\star$ with $|Y^\star|\le n^\star/9$, and let $Y$, $L$,
$A$, $B$ be as in Lemma~\ref{lem:T}. Note $n_0\le N_0\le n/60\le m/120$. Suppose
$L\le m/20$. Then at most $m/20$ clauses cross $A$, at least
$\frac m3-\frac m{20}-\frac m{360}\ge\frac m6$ clauses lie inside $A$, and at least
$\frac{2m}9-\frac m{20}-\frac m{360}=\frac{61m}{360}\ge\frac m6$ clauses lie inside
$B\subseteq[n]\setminus A$. This contradicts Lemma~\ref{lem:P} with $S=A$. Hence
$d(Y)=L>m/20>0$, so $Y\ne\emptyset$, and by (A2$'$) with $j=|Y|$,
\[
\frac{\Delta n}{20}\ <\ 3e^2\Delta\,|Y|\log\frac{en}{|Y|},\qquad\text{that is}\qquad
y\log(e/y)>\frac1{60e^2}\quad\text{with }y=|Y|/n\in(0,1].
\]
As $y\mapsto y\log(e/y)$ is increasing on $(0,1]$, $y\ge\eta^\star$. So
$|Y^\star|\ge|Y|\ge\eta^\star n$. If $|Y^\star|>n^\star/9$, then
$|Y^\star|>3\Delta n/9\ge\eta^\star n$ as well. Finally $n^\star=3\Delta n+n_0\le(3\Delta+1)n$.
\end{proof}

\begin{corollary}[The scheme applied to $F^\star$]\label{cor:fstar}
Fix $\Delta>3$ (or only $\Delta\ge2$), let $F$ be a random 3-CNF as in
Lemma~\ref{lem:random}, and let $F^\star$ be its bounded-occurrence formula built with
Tovey's cycles. On an event of probability $1-o(1)$, $F^\star$ belongs to
family~\textup{(i)} of Corollary~\ref{cor:main} with $n^\star$ in place of $n$, $b=3$ and
$c=\eta^\star/(3\Delta+1)$. Hence, with $\gamma^\star=c_3\,\eta^\star/(3\Delta+1)$ and
$n^\star\ge n_0(\gamma^\star)$, every clause-product interpolation scheme for $F^\star$ has
product-type charge at least
\[
2^{(\gamma^\star (n^\star)^{2/3}-1)/2}\ \ge\ 2^{(\gamma^\star(3\Delta n)^{2/3}-1)/2}
=2^{\Omega_\Delta(n^{2/3})},
\]
no scheme for $F^\star$ has $s_t,q_t,\rho_t\le\sqrt{n^\star}\log^\ell n^\star$ at all points,
and variable-uniform schemes for $F^\star$ have product-type charge $2^{\Omega_\Delta(n)}$.
\end{corollary}

\begin{proof}
Every variable of $F^\star$ occurs at most $3$ times, and Theorem~\ref{thm:fstar} gives
$\sep(F^\star)\ge cn^\star$. Corollary~\ref{cor:main}, family~(i), applies to $F^\star$:
Theorems~\ref{thm:uniform} and~\ref{thm:general} hold for every CNF, and the clause width
plays no role.
\end{proof}

\begin{remark}[Sparsification leaves]\label{rem:leaves}
Let $G$ be obtained from $F$ by deleting at most $r$ clauses and adding at most $u\le r$
unit clauses, and let $G^\star$ be a copy expansion of $G$ (with $d_v$ in (E3) counting the
occurrences in $G$, unit clauses included). The proof of Lemma~\ref{lem:T} gives, with
$L=d_F(Y)$: every clause of $F$ lies inside $A$, inside $B$, contains a variable of $Y$, or
is deleted; at least $\frac m3-L-\frac{r+u+n_0}3$ clauses of $F$ lie inside $A$ and at least
$\frac{2m}9-L-\frac{r+u+n_0}3$ inside $B$; here $n^\star\ge3(m-r)$, the fibres over $A\cup Y$
have at most $3m_A+2L+L+u+n_0$ elements, and $n_0\le N_0+3r$. If $r\le m/600$ and
$N_0\le m/600$, the assumption $L\le m/21$ contradicts Lemma~\ref{lem:P} (at most
$L+r\le m/20$ crossing clauses; at least $\frac{2m}9-\frac m{21}-\frac m{300}\ge\frac m6$
clauses inside $B$). So, with high probability, simultaneously for all such $G$,
$\sep(G^\star)\ge\eta^{\star\star}n$ with $\eta^{\star\star}=0.000228\ldots$ the solution
of $y\log(e/y)=1/(63e^2)$, and $n^\star(G)\le(3\Delta+1)n$. By the description of the
preprint's sparsification in the proof of \cite[Theorem~8.3(b)]{DL}, the leaf on the
first-child path is such a $G$ with $r\le\eta(\theta_2)n$, $u\le r/\theta_2$ and
$\eta(\theta_2)\to0$; under that description, Corollary~\ref{cor:fstar} extends to the
leaves once $\eta(\theta_2)\le\Delta/600$. We rely on~\cite{DL} for this description and
did not check it against the preprint.
\end{remark}
```

## 2. Notes for the editor

1. **Hypotheses actually used.** $F$ exactly-3-uniform with distinct variables per clause
   (true in the model); $F^\star$ satisfies (E1)-(E3); $b=3$ only enters Corollary. The density
   needed is $\Delta\ge2$ (Lemma P needs $\Delta>\ln 2/\kappa=1.766$; (N) needs $e^{-3\Delta}<1/120$,
   i.e. $\Delta>1.6$). $\Delta>3$ is therefore *not* needed for $F^\star$; it is needed only by (A1).
2. **Constants.** $H(1/6)=0.450561$, $H(1/20)=0.198515$, $\kappa=0.392437$,
   $\eta^\star=2.418\times10^{-4}$, $\eta^{\star\star}=2.290\times10^{-4}$;
   $\gamma^\star=c_3\eta^\star/(3\Delta+1)\approx 8.2\times10^{-6}/(3\Delta+1)$. Very small but explicit.
   (For comparison the paper's own $\eta_\Delta$ solves $y\log(e/y)=1/(9e^2)$: $\eta_\Delta=0.002098$,
   independent of $\Delta$ since $\Delta/(3C_\Delta)=1/(9e^2)$; the paper could say so.)
3. **Arithmetic in Theorem.** $n_0/3\le m/360$; $1/3-1/20-1/360=101/360\ge1/6$;
   $2/9-1/20-1/360=(80-18-1)/360=61/360\ge60/360$.
4. **Variant constructions.** If the preprint copies only variables with more than 3
   occurrences (fibre of size 1 otherwise), (E3) becomes $|\pi^{-1}(v)|\le\max(d_v,1)$ and
   $n^\star\ge3m-3n\ge 3\xi m$ with $\xi=1-1/\Delta>2/3$. Lemma T then gives
   $m_A\ge\xi m/3-L-n_0/3$, $m_B\ge2\xi m/9-L-n_0/3\ge 4m/27-\dots$, and Lemma P must be used
   with $(\beta,\varepsilon)=(1/8,1/50)$: $\kappa=0.31652$, valid for $\Delta>2.19$. Same
   conclusion with smaller constants (not written out in full). If equality gadgets use
   auxiliary variables, map them to $v$ by $\pi$; (E3) then needs $|\pi^{-1}(v)|\le K d_v$ and the
   constants change with $K$ (not done).
5. **What is not proved.** That the preprint's $F^\star$ is exactly the Tovey expansion with one
   copy per occurrence (taken from DL's wording); the leaf description in Remark (taken from DL's
   proof of 8.3(b)); $\eta(\theta_2)\to0$ w.h.p. (DL's claim).
6. **Comparison with DL Theorem 8.3.** With the Corollary the comparison becomes direct:
   same formula $F^\star$, same $n^\star=\Theta_\Delta(n)$; DL: $2^{\Omega_\Delta(\sqrt n)}$ in the
   DP-type (hence product-type) model; here $2^{\Omega_\Delta(n^{2/3})}$ in the product-type model
   only. Table~\ref{tab:compare}, the "What this paper does not show" item on $F^\star$, the
   sentence after Corollary~\ref{cor:main}, and Open problem 2 would have to be updated.

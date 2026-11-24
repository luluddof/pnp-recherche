# Outils réutilisables (`PnpLean.Outils`)

Ce fichier liste les résultats de ce projet qui ne parlent **ni de SAT ni des articles**, et que d'autres chercheurs peuvent réutiliser. Un seul fichier les rassemble : `PnpLean/Outils.lean` (espace de noms `PnpOutils`).

```bash
cd lean && lake exe cache get && lake build PnpLean.Outils
```

## Les deux conditions d'entrée

1. **Juste.** L'outil est prouvé en Lean 4 + Mathlib, sans `sorry`, `admit`, `axiom` ni `native_decide`. Il dépend seulement de `propext`, `Classical.choice`, `Quot.sound`. Les programmes Python du projet ne respectent pas cette condition : ils sont testés, pas prouvés. Ils ne sont donc **pas** dans cette liste.
2. **Absent de Mathlib.** Chaque énoncé a été cherché dans Mathlib (novembre 2026) et n'y est pas.

**Limite de la condition 2.** La recherche porte sur Mathlib seulement. Les autres bibliothèques (Isabelle, Rocq, autres projets Lean) n'ont pas été examinées. Les mathématiques sont classiques : ce qui est nouveau, c'est la preuve vérifiée par machine, pas le théorème.

## Liste des outils

### 1. Matrices de polynômes

| Nom | Énoncé | Limite |
|---|---|---|
| `pow_dvd_det_of_nullity` | Pour Φ carrée sur ℚ[X] et un point t : (X − t)^(dim ker Φ(t)) divise det Φ. Cas particulier connu : multiplicité géométrique ≤ multiplicité algébrique (Φ = X·1 − A, dans Mathlib sous la forme `finrank_eigenspace_le`). | Sur ℚ seulement |
| `sum_nullity_le` | Si les coefficients de Φ ont un degré ≤ R et Φ(0) = 1, alors Σ_{t∈T} dim ker Φ(t) ≤ n·R. | Sur ℚ |
| `sum_nullity_le_of_det_ne_zero` | Le même énoncé avec l'hypothèse plus faible det Φ(0) ≠ 0. | Sur ℚ |
| `deficiency_le_nullity_add` | Borne de déficit par blocs : Σ_γ (colonnes de γ − lignes de γ) ≤ dim ker A + Σ_γ dim R_γ. | Sur ℚ |

### 2. Probabilités sur un produit fini, sans théorie de la mesure

| Nom | Énoncé | Limite |
|---|---|---|
| `bernoulli_anticonc`, `c0` | Pour une somme X d'indicatrices indépendantes, de variance V, et tout réel a : c₀·min(V, √V) ≤ E\|X − a\|, avec c₀ = (1 − e⁻²)/4. | Espaces finis, mesure uniforme |
| `uavg_pi_mul` | Indépendance : si g ne dépend pas de la coordonnée a, alors E[f(ω_a)·g(ω)] = E[f]·E[g]. | Voir « Non inclus » |

### 3. Chemins polynomiaux de matrices sur le cube booléen

| Nom | Énoncé | Limite |
|---|---|---|
| `sum_rank_ge_multipoint` | J de degré ≤ d avec J(0) = X^{⊗ℓ} ; points distincts non nuls τ_k de directions c(k) ; J(τ_k) = K_k + E_k avec K_k[x,y] = 0 si x et y diffèrent en c(k). Alors Σ_x max(0, N − d − 2μ(x)) ≤ Σ_k rang E_k. | Tout corps |
| `rank_eval_zero_le` | A matrice de polynômes de forme quelconque ; préfacteurs f_u, g_v de valeur 1 en 0 ; deg(f_u·A_uv·g_v) < N. Si en N points distincts les préfacteurs annulent A(τ_k) − E_k, alors rang A(0) ≤ Σ_k rang E_k. | Tout corps |
| `rank_submatrix_le_sum_rank` | Version avec étiquettes dans Q^ℓ et budgets : le rang d'une sous-matrice de A(0) est ≤ Σ_k rang E_k. | Tout corps |
| `sum_rank_ge_N2` | Forme à deux paramètres, un point par direction. | Tout corps |
| `sum_rank_ge_block` | Forme par blocs : (\|T\| − d)·2^(ℓ−\|S\|) ≤ Σ_t rang E_t. | Sur ℚ |
| `four_bad_edges` | Sur le cube Q_{d+2}, au moins 4 mauvaises arêtes. | Sur ℚ |
| `card_subcube` | Un sous-cube a 2^(ℓ−\|S\|) points. | — |

## Non inclus, et pourquoi

- **`UAvg` (moyennes uniformes)** : Mathlib a déjà `Finset.expect` (notation 𝔼). Notre définition `uavg` fait double emploi. Seuls les deux résultats qui s'appuient dessus sont listés.
- **`rank_add_le'`, `rank_sum_le`** : déjà dans Mathlib (`Matrix.rank_add_le`).
- **Les théorèmes des articles** (`theoremA`, `theoremB`, `Scheme`, …) : ils dépendent de la définition 5.1 de Dias et Lima. Ils ne sont pas génériques.
- **Le code Python** (`code/bench/`, scripts des agents) : testé, pas prouvé.

## Pour en faire une vraie bibliothèque (pas fait)

1. Généraliser la section 1 de ℚ à un corps quelconque.
2. Réécrire la section 2 avec `Finset.expect` au lieu de `uavg`.
3. Choisir une licence. **Sans licence, personne n'a le droit de réutiliser le code.** C'est une décision de l'auteur.
4. Proposer `pow_dvd_det_of_nullity` et `bernoulli_anticonc` à Mathlib.

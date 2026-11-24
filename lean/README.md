# Vérifications Lean, classées par article

Lean 4.37.0-rc1 + Mathlib. Chaque dossier correspond à un article.

```bash
cd lean && lake exe cache get && lake build
```

Cette commande compile tout ce qui est **vérifié** (les dossiers `Article1/` et `Article2/`). Les dossiers `EnCours/` contiennent du travail non terminé : ils ne font pas partie de la cible par défaut.

**Règle :** un fichier hors de `EnCours/` ne contient ni `sorry`, ni `admit`, ni `axiom`, ni `native_decide`. Ses théorèmes dépendent seulement des axiomes standards de Lean (`propext`, `Classical.choice`, `Quot.sound`).

---

## Article 1 — « Measurement Rows Absorb at Most a √s Factor of the Cut »

Article : `resultats/article/main.tex` (anglais) et `main_fr.tex` (français).

### Vérifié : `PnpLean/Article1/`

| Fichier | Contenu | Énoncé de l'article |
|---|---|---|
| `WitnessDet.lean` | `sum_nullity_le` : Σ_t dim ker Φ(t) ≤ nR ; `deficiency_le_nullity_add` : borne de rang par blocs | Cœur du lemme du déterminant témoin |
| `Scheme.lean` | Structure `Scheme` (définition 5.1 de Dias-Lima), `lemma1`, `factR` | Lemme du déterminant témoin (complet), comptage de Hermite |
| `TheoremB.lean` | `avg_lemma1`, `exists_point`, `lemma3`, `theoremB`, `theoremB_bounded` | Lemme du séparateur ; théorème principal avec le lemme d'espérance du déficit en hypothèse (`Lemma2Hyp`) |
| `TheoremA.lean` | `theoremA` | Théorème pour les déformations uniformes en la variable, sans condition |
| `UAvg.lean` | Moyennes uniformes sur des produits finis ; `uavg_pi_mul` (indépendance) | Outil pour les lemmes probabilistes, sans théorie de la mesure |
| `Anticonc.lean` | `bernoulli_anticonc` : c₀·min(V, √V) ≤ E\|X − a\| avec c₀ = (1 − e⁻²)/4 ; symétrisation, 2ᵉ et 4ᵉ moments, Hölder | Lemme d'anti-concentration |
| `Extras.lean` | `theoremB_unconditional'`, `theoremA'`, `exists_point'`, `lemma1_occ`, `WitnessDet.natDegree_det_le_card` | **Les deux théorèmes sans l'hypothèse `n ≤ M`** (borne de degré par le nombre de variables qui apparaissent) |
| `TrivialScheme.lean` | `trivialScheme`, `exists_trivial_scheme_sqrt` | Le schéma trivial de la proposition sur le modèle de type PD : déformations constantes, q = 0, reconstruction de Lagrange, groupes de taille ≤ s avec n ≤ s², ρ_t = n |
| `Lemma2.lean` | `avg_home_ind`, `avg_Icard`, `half_sum_abs_le_Edef`, `anticonc_le_Edef`, **`lemma2`**, **`theoremB_unconditional`**, `theoremB_bounded_unconditional` | Lemme d'espérance du déficit ; **théorème principal sans hypothèse sur le lemme 2** |

Cible : `lake build PnpLean.Article1.Lemma2` (et `PnpLean.Article1.TheoremA`).

**État :** tous les énoncés **déterministes et probabilistes finis** de l'article sont vérifiés : lemme du déterminant témoin, comptage de Hermite, anti-concentration, espérance du déficit, lemme du séparateur, les deux théorèmes. `theoremB_unconditional'` (dans `Extras.lean`) n'a plus que l'hypothèse `hW` (W minore le nombre de séparation pondéré : c'est une définition).

### Non formalisé

Le lemme sur les 3-CNF aléatoires (sep(F) ≥ c·n avec forte probabilité), le corollaire asymptotique, la largeur arborescente et les charges 2^s et 2^n du modèle de type PD (le schéma trivial lui-même est construit), tout ce qui concerne F\*.

---

## Article 2 (en préparation) — erreurs par point dans la méthode d'interpolation d'OpenAI

Pas encore d'article. Sources : `code/gen6/g6d/`, `code/gen7/g7a/`, `code/gen8/g8d/` ; relecture : G10b.

### Vérifié : `PnpLean/Article2/`

| Fichier | Contenu | Énoncé |
|---|---|---|
| `PerPoint.lean` | `sum_nullity_le_of_det_ne_zero`, `sum_rank_ge`, `sum_rank_error_ge`, **`g6d_determinant_lemma`** | Lemme du déterminant par bloc : pour J de degré ≤ d avec J(0) = X^{⊗ℓ}, et J(t) = K_t + E_t avec K_t « faiblement agréable » dans une direction de S : Σ_{t∈T} rang(E_t) ≥ (\|T\| − d)·2^(ℓ−\|S\|). Vaut dans le cas général (non commutatif). |
| `QdPlus2.lean` | `two_splittings`, `four_le_card_bad`, **`qd2_lemma`** | Lemme Q_{d+2} : sur le cube de dimension d+2, avec des polynômes de sommet de degré ≤ d valant 1 en 0 et des points τ distincts non nuls, il y a au moins 4 mauvaises arêtes. |
| `SqrtLog.lean` | **`sum_rank_ge_N`**, `core`, `N_eq` | **Théorème principal (candidat) de l'article 2 :** pour J de degré ≤ d avec J(0) = X^{⊗ℓ}, des points τ_k distincts non nuls (un par direction) et J(τ_k) = K_k + E_k avec K_k[x,y] = 0 si x_k ≠ y_k : si 2m + d + 1 ≤ ℓ, alors N(ℓ, m) ≤ Σ_k rang(E_k), avec N(ℓ, m) = Σ_{j≤m} Σ_{i≤j} C(ℓ, i). Sur ℚ. |
| `Consequences.lean` | `exists_rank_ge`, `avg_le_max_rank`, `sum_choose_le_mul_of_rank_le`, `degree_ge_of_rank_le`, `N_one`, `N_two`, `sum_rank_ge_add_two` | Conséquences finies du théorème principal : forme max/moyenne, bornes sur N, et « si tous les rangs sont ≤ r et ℓ·r < Σ_{i≤m} C(ℓ,i), alors ℓ ≤ 2m + d » |
| `AnyField.lean` | **`sum_rank_ge_N2`**, `sum_rank_ge_N` (tout corps), `sum_rank_ge_N_rat` | **Le théorème principal sur un corps quelconque**, sous la forme à deux paramètres : si ma + mc + d + 1 ≤ ℓ, alors N2(ℓ, ma, mc) = Σ_{|x| ≤ mc} (ma + 1 − |x|) ≤ Σ_k rang(E_k). La preuve n'utilise pas de positivité (un seul vecteur test). Le cas ℚ de `SqrtLog.lean` en est un cas particulier. |
| `MultiPoint.lean` | **`sum_rank_ge`**, `sum_rank_ge_half`, `sum_rank_block01_ge`, `sum_rank_ge_one`, `placement_le`, `sum_rank_ge_balanced` | **La forme la plus forte du théorème, sur tout corps, avec plusieurs points par direction :** pour des points indexés par ι, de directions c : ι → Fin ℓ, Σ_x (\|ι\| − d − 2·μ(x)) ≤ Σ_k rang(E_k), où μ(x) = nombre de points k avec x(c k) = 1. Vaut deux fois la borne de `SqrtLog.lean` quand \|ι\| − d est pair. Seule la moitié « 01 » de la condition sur K est utilisée. `sum_rank_ge_balanced` : forme moyennée (le placement des points n'aide pas). |

Cibles : `lake build PnpLean.Article2.MultiPoint`, `PnpLean.Article2.AnyField`, `PnpLean.Article2.Consequences` et `PnpLean.Article2.QdPlus2`. `PerPoint.lean` importe `PnpLean.Article1.WitnessDet` (même argument de déterminant).

### Non formalisé

Les conséquences asymptotiques avec logarithmes (seuil sur d, temps de calcul dans le modèle de coût), l'énoncé entier « le placement équilibré est le minimum » (énuméré seulement), le théorème au niveau de la sortie (effondrement par L, R′), la réduction de la classe commutative (rang E_k = nombre de mauvaises arêtes), les corollaires (a) et (b) de G6d, la borne ρ_moy ≥ 2^(−d)/(d+2), la restriction du lemme Q_{d+2} à un sous-cube d'un cube plus grand, et le lien avec l'algorithme d'OpenAI lui-même.

---

## Compatibilité avec la version 2 de l'article 1

La version 2 de l'article (envoyée en relecture le 24 novembre 2026) cite `lean/PnpLean/` et la commande `lake build PnpLean.TheoremA`. Les fichiers `PnpLean/WitnessDet.lean`, `Scheme.lean`, `TheoremB.lean` et `TheoremA.lean` sont de simples renvois vers `PnpLean/Article1/`, pour que cette commande marche toujours. L'état exact du dépôt pour la version 2 est le commit `d597db3`.

---

## Ajouter un article

Créer `PnpLean/ArticleN/` (vérifié) et `PnpLean/ArticleN/EnCours/` (travail en cours), puis ajouter une section ici. Ajouter la cible vérifiée dans `PnpLean.lean`.

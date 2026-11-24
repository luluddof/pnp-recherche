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
| `TheoremB.lean` | `avg_lemma1`, `exists_point`, `lemma3`, `theoremB`, `theoremB_bounded` | Lemme du séparateur ; théorème principal, **avec le lemme d'espérance du déficit en hypothèse** (`Lemma2Hyp`) |
| `TheoremA.lean` | `theoremA` | Théorème pour les déformations uniformes en la variable, **sans condition** |

Cible : `lake build PnpLean.Article1.TheoremA`

### En cours : `PnpLean/Article1/EnCours/`

| Fichier | But | État |
|---|---|---|
| `UAvg.lean` | Moyennes uniformes sur des produits finis (pour l'identité « moyenne sur les cartes témoins = choix indépendants ») | Compile ; non relié au théorème principal |
| `Anticonc.lean` | Anti-concentration (lemme d'espérance du déficit), pour enlever l'hypothèse `Lemma2Hyp` | **Ne compile pas encore** |

### Non formalisé

Le lemme d'anti-concentration et le lemme d'espérance du déficit (en cours), le lemme sur les 3-CNF aléatoires, le corollaire asymptotique, la proposition sur le modèle de type PD.

---

## Article 2 (en préparation) — erreurs par point dans la méthode d'interpolation d'OpenAI

Pas encore d'article. Sources : `code/gen6/g6d/`, `code/gen7/g7a/`, `code/gen8/g8d/` ; relecture : G10b.

### Vérifié : `PnpLean/Article2/`

| Fichier | Contenu | Énoncé |
|---|---|---|
| `PerPoint.lean` | `sum_nullity_le_of_det_ne_zero`, `sum_rank_ge`, `sum_rank_error_ge`, **`g6d_determinant_lemma`** | Lemme du déterminant par bloc : pour J de degré ≤ d avec J(0) = X^{⊗ℓ}, et J(t) = K_t + E_t avec K_t « faiblement agréable » dans une direction de S : Σ_{t∈T} rang(E_t) ≥ (\|T\| − d)·2^(ℓ−\|S\|). Vaut dans le cas général (non commutatif). |
| `QdPlus2.lean` | `two_splittings`, `four_le_card_bad`, **`qd2_lemma`** | Lemme Q_{d+2} : sur le cube de dimension d+2, avec des polynômes de sommet de degré ≤ d valant 1 en 0 et des points τ distincts non nuls, il y a au moins 4 mauvaises arêtes. |

Cible : `lake build PnpLean.Article2.QdPlus2`. `PerPoint.lean` importe `PnpLean.Article1.WitnessDet` (même argument de déterminant).

### Non formalisé

La réduction de la classe commutative (rang E_k = nombre de mauvaises arêtes), les corollaires (a) et (b) de G6d, la borne ρ_moy ≥ 2^(−d)/(d+2), la restriction du lemme Q_{d+2} à un sous-cube d'un cube plus grand, et le lien avec l'algorithme d'OpenAI lui-même.

---

## Ajouter un article

Créer `PnpLean/ArticleN/` (vérifié) et `PnpLean/ArticleN/EnCours/` (travail en cours), puis ajouter une section ici. Ajouter la cible vérifiée dans `PnpLean.lean`.

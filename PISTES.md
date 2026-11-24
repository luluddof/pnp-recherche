# Pistes de recherche

Ce fichier liste toutes les pistes : celles en cours, celles pas encore explorées et celles fermées.
Mise à jour : 24 novembre 2026.

**Rappel :** aucune de ces pistes ne mène directement à P = NP. Ce sont des étapes qui peuvent donner des résultats vérifiables.

Légende : 🔄 en cours · 🆕 pas encore explorée · ⏸️ en attente · ✅ terminée · ❌ fermée (mur prouvé ou mesuré)

---

## 0. L'échelle vers P = NP (temps pour 3-SAT)

| Marche | Temps | État |
|---|---|---|
| 0 | 2^(0,39 n) | Record avant octobre 2026 |
| 1 | 2^(n/√log n) | Atteint par OpenAI et Claude (annoncé, non validé) |
| 2 | 2^(n/log n) | Mur des méthodes actuelles (générations 4 à 6) |
| 3 | 2^(n^(1−ε)) | Demande une idée nouvelle (hors structure produit) |
| 4 | 2^(√n) | — |
| 5 | Quasi-polynomial puis polynomial (= P = NP) | La plupart des experts pensent que cette marche n'existe pas |

Priorité recommandée : pistes A et B, puis Q1\*.

## 1. En cours (générations 10 et 11)

| Id | Piste | Question précise | Origine |
|---|---|---|---|
| G11a | 🔄 **Relecture adverse du théorème sep(F\*)** | Nécessaire avant la version 3 de l'article 1 | G10f |
| G11b | 🔄 **Article 2, résultat final : compromis (I′) et (II)** | Prouver l'optimalité de √log n dans la classe commutative, ou trouver une construction qui la contredit | G8d |
| G11c | 🔄 **Article 2 : cas non commutatif** | Étendre le lemme Q_{d+2} et la réduction « rang = mauvaises arêtes » au cas général | G6d, G10b |
| G10d | 🔄 **Relecture de G9c** | Théorèmes A et B sur le comptage des ensembles indépendants | G9c |

**Consigne de l'utilisateur (24/11) :** dès que les générations 8 et 9 sont terminées, lancer une **génération 10** sur les pistes « pas encore explorées » (section 2) et « en attente » (section 3).

## 2. Pas encore explorées

| Piste | Question précise | Difficulté | Origine |
|---|---|---|---|
| 🆕 **Problème de Naslund / Hurwitz** | Prouver P(N) ≥ N^(1+ε) ; équivaut à un plancher ρ > 1/2 dans l'algorithme de Claude | Très difficile (circuits non commutatifs) | N1, G5b |
| ❌ **Tests de parité économiques (Q2)** — fausse sur la famille de référence (G8b, lemme CA) | Tests avec O(n·2^O(k)) incidences au total dans l'annexe D de Claude | Technique | G6c |
| 🆕 **Erreurs par point, cas non commutatif** | Étendre le lemme de G6d au-delà de la classe commutative | Moyenne | G6d, G7a |
| 🆕 **Forme parité et contraction de tenseurs (Limits)** | Analyser les limites de ces deux voies d'OpenAI, non couvertes par Limits ni par notre note | Difficile | G4c, G6a |
| 🆕 **Sparsification avec variables d'extension** | Échappatoire à la proposition L de G6c | Moyenne | G6c, G7d |
| 🆕 **Lemme 4 en Lean** | Graphes aléatoires | Très difficile (peu de choses dans Mathlib) | G9L |

| 🆕 **Formes closes aux niveaux profonds (suite de G8c)** | La valeur des sous-appels de l'algorithme de Claude sur une entrée « épinglée » a une forme close calculable en temps polynomial au niveau 1 (vérifié). Si cela reste vrai aux niveaux profonds, Extended Frege réfuterait toute 3-CNF creuse insatisfaisable en taille 2^(δn). Vérifier aussi l'hypothèse (S) : degré ≤ 1 en chaque table de poids | Moyenne à difficile | G8c |


| 🆕 **Dépasser n^(2/3) dans le modèle produit** | La méthode du déterminant témoin est limitée à n^(2/3) (prop. S de G8a). Il faut une contrainte nouvelle qui couple les points ou qui utilise la structure polynomiale. Ou construire un schéma de coût n^(2/3) : il a besoin de Ω(n^(2/3)) points | Difficile | G8a |
| 🆕 **Un modèle DP valide** | Le modèle corrigé de G7b (largeur W_t) ne minore que l'élimination par groupes entiers (G10a) ; il faut un modèle qui résiste au découpage ligne par ligne, par exemple une largeur sur le graphe des lignes | Moyenne | G10a |
| 🆕 **Borne de Rado pour le modèle DP corrigé** | nullité ≥ max_S(\|S\| − dim Σ U_ω(v)) contient le lemme 1 et le lemme H ; pourrait enlever l'hypothèse de régularité | Moyenne | G8a |

| 🆕 **Compromis (I′) et (II) pour la méthode d'OpenAI** | (I′) ρ(Cδ², δ) ≥ δ^(−a) ; (II) ρ(ℓ₁+ℓ₂, d₁+d₂) ≥ ρ(ℓ₁,d₁)·ρ(ℓ₂,d₂). Ensemble : √log n optimal à (log n)^o(1) près | Difficile | G8d |


| 🆕 **Pièces hors du cadre « un seul espace gradué »** | Seule façon restante d'obtenir ρ < 1/√2 : des pièces qui ne passent pas par les composantes graduées d'un seul espace (ou des profils de croissance non concaves). Lié au problème de Hurwitz : P(N) ≤ N^(1+ε) | Très difficile | G9a |

| 🆕 **Termes non produits (« intriqués »)** | Seule hypothèse non testée du mur 2^(m/log m). Il faut E1 : décomposition jointe d'un tenseur de grille 4×h (termes colonnes ≤ (1/2+o(1))·S bits, termes lignes de rang 2^o(S)) et E2 : évaluation récursive des parties fusionnées. Aucune borne connue n'exclut 2^(S^(1−ε)) dans la classe générale | Très difficile ; c'est la question centrale pour la marche 3 | G9d |

| 🆕 **Pièces sur un corps fini fixe avec ρ < 1/√2** | Seule liberté propre à la décision. Sur GF(2) : meilleur connu 0,811 (classe restreinte) ; recherche exhaustive à λ = 2 non terminée | Difficile | G9b |

| 🆕 **Noyau de Q1\* : graphes bipartis d-réguliers aléatoires** | Compter les ensembles indépendants de taille entre n/(log d·log log d) et n/2 ; toutes les méthodes connues paient n·log log n/√log n ; reformulation : un bit partagé par r groupes doit coûter O(1) | Difficile | G9c |

## 3. En attente

| Piste | Raison de l'attente |
|---|---|
| ⏸️ **Version 3 de l'article** | **Prête à lancer** (G10e est terminé) : (1) corrections de G10f (facteur b superflu ; dire que la borne de Dias-Lima est meilleure pour les petits groupes ; précisions de l'annexe A ; 4 attributions à adoucir) ; (2) le théorème sep(F\*) ≥ c·n\* et son corollaire, **après relecture adverse** ; (3) mise à jour de la section Lean |
| ⏸️ **Relecture adverse du théorème sep(F\*)** (G10f) | Niveau P pour l'instant |
| ⏸️ **Rédaction de l'article 2** (erreurs par point) | **Décision de l'auteur (24/11) : attendre un résultat final** (pas de note de résultats partiels). Résultat final visé : l'optimalité de √log n pour la méthode d'OpenAI (pistes (I′) et (II)), ou au contraire une construction qui bat √log n. Les deux lemmes centraux sont déjà relus et en Lean (`lean/PnpLean/Article2/`) |
| ⏸️ **Note « plus long chemin / plus long cycle en 2^o(n) sur les graphes creux »** (candidat article 4) | Preuve confirmée par G10c ; conditionnel aux prépublications ; doit citer Wu-Wang (arXiv 2611.02548) ; résultat modeste. Décision de l'auteur |
| ⏸️ **Relecture de l'article par un expert humain** | Auteurs naturels à contacter : Dias et Lima. Décision de l'auteur. |

| Piste | Raison de l'attente |
|---|---|

## 4. Terminées (résultats)

| Piste | Résultat |
|---|---|
| ✅ **Lean : théorème B sans hypothèse** (G10e) | Lemme 2 (anti-concentration, espérance du déficit) et identité de moyenne prouvés ; `theoremB_unconditional` |
| ✅ Relecture + nouveauté de G9e (G10c) | Preuve confirmée ; Steiner déjà revendiqué (Wu-Wang) ; plus long chemin/cycle : semble nouveau |
| ✅ Relecture + Lean des lemmes « erreurs par point » (G10b) | Lemme de G6d et lemme Q_{d+2} confirmés et **prouvés en Lean** (`lean/PnpLean/Article2/`) ; la conséquence (c) de G6d n'est pas prouvée en général |
| ✅ Relecture adverse de G8a et G7b (G10a) | B2, lemme H, théorèmes D et F confirmés ; prop. S correcte à 2 constantes près ; une clause du théorème C réfutée ; le modèle DP corrigé n'est que partiellement valide |
| ✅ Relecture des points non relus de l'article (G10f) | Note 7.1(ii) correcte (2 corrections d'équité) ; annexe A vérifiée ; bibliographie sans entrée fausse |
| ✅ **Article en français** (G9F) | `resultats/article/main_fr.tex` + PDF (28 pages) ; traduction fidèle, contrôle automatique des comptes (théorèmes, équations, références) |
| ✅ **Article, version 2** (G9W) | `resultats/article/` : 26 pages, LaTeX + PDF ; les 10 corrections obligatoires du rapport G8R sont traitées ; hypothèse M ≥ n supprimée sur papier (Lean garde n ≤ M) |
| ✅ Limits §9, modèle produit (G6a) | Réponse négative, 2^Ω(n^(2/3)) ; confirmée par relecture adverse ; note rédigée |
| ✅ Auteur humain de la note | **Lucas Requena**, étudiant en informatique (« Computer Science student »), contact : l.requena@outlook.fr (avec une mention claire du travail fait par les agents IA) |
| ✅ Lemme 1 en Lean (G8L) | Prouvé sans `sorry` |
| ✅ Partie déterministe en Lean (G9L) | Schéma, lemme 1 complet, fait R, lemme 3, théorème A (sans condition), théorème B (avec le lemme 2 en hypothèse) ; sans `sorry` |
| ✅ Hurwitz P(N) ≤ N^√2 (G4a) | Correct, mais déjà publié par d'autres |
| ✅ Réfutations Extended Frege (G8c) | **Conditionnel :** sous (S) et des preuves EF courtes de « Pin_F(b) », EF réfute F en taille 2s + poly(n)·2^(δn). Le calcul de l'algorithme est de **profondeur constante** (nouveau) : le certificat IPS peut être une formule de profondeur constante. Pas de route vers TC⁰-Frege (descente séquentielle). |
| ✅ Conséquences de la chute de l'ETH (G9e) | Arbre de Steiner, plus long cycle/chemin sur graphes creux en 2^(O(m log log m/√log m)) (conditionnel, non relu, nouveauté non établie) |
| ✅ #IS sans log d sur des familles (G9c) | Graphes à petit α, graphes (n,d,λ), d-réguliers aléatoires : 2^(O(en)) (niveau P, non relu) |
| ✅ Système de preuve IPS (G4d) | Réfutations IPS de taille 2^O(n/√log n) |
| ✅ Plancher ρ_d ≥ 1/d (G6b) | Confirmé sous conditions |

## 5. Fermées (mur prouvé ou mesuré)

| Piste | Cause |
|---|---|
| ❌ Algorithmes locaux, LP, SOS, spectral, PC, plans de coupe, parités courtes, variables d'extension simples | Bornes inférieures connues ou croissance exponentielle mesurée (générations 0 à 3) |
| ❌ Anneaux Z_m, Z, Q | Réduction modulo p (preuve) |
| ❌ Certificats de Feige-Kim-Ofek trouvables en temps polynomial | Égalent seulement le seuil spectral n^1,5 |
| ❌ A – Meilleures pièces par espaces gradués à croissance lente (G9a) | Comptage de rang : α ≥ 1/2, donc ρ ≥ 1/√2 dans le cadre de l'annexe B ; la localité se casse si on restreint les modes |
| ❌ D – Gadgets à 3 parties, multiplication rapide de matrices, conjecture du rang asymptotique, instances non cliques, réutilisation des égalités (G9d) | Gain constant au mieux ; « loi harmonique » : tout schéma à termes produits reste à 2^Ω(m/log m) (esquisse, dans le modèle) |
| ❌ B – Décider au lieu de compter (G9b) | Les meilleures pièces ne marchent qu'en grande caractéristique ; sur GF(2) le meilleur ρ connu est 0,811 > 1/√2 ; Valiant-Vazirani et le hachage n'aident pas |
| ❌ Géométries de récursion (grilles, boîtes de dimension d) | 2^Ω(m/log m) dans le modèle (G4b, G6b) |
| ❌ Structures non produit pour OpenAI (expanseurs, codes) | Rang ≥ K/2 (G5a) |
| ❌ Combinaisons OpenAI + Claude | Même mur 2^(n/log n) (G4c) |
| ❌ Réfuter Feige (m = O(n)) | Bornes inférieures SOS ; hors de portée |

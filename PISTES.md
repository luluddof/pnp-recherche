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

## 1. En cours (génération 9)

| Id | Piste | Question précise | Origine |
|---|---|---|---|
| G9c | 🔄 **C – Q1\*** | Compter les ensembles indépendants d'un graphe presque creux en 2^O(n/√log n) ; idées : hybride, autre encodage, identité de l'étoile | G8b |

**Consigne de l'utilisateur (24/11) :** dès que les générations 8 et 9 sont terminées, lancer une **génération 10** sur les pistes « pas encore explorées » (section 2) et « en attente » (section 3).

## 2. Pas encore explorées

| Piste | Question précise | Difficulté | Origine |
|---|---|---|---|
| 🆕 **Lemme 4 pondéré pour F\*** | Prouver la borne de séparateur pour la formule F\* de la prépublication et pour les feuilles de la sparsification, pour comparer avec le théorème 8.3 de Dias-Lima | Moyenne | Rapport G8R |
| 🆕 **Problème de Naslund / Hurwitz** | Prouver P(N) ≥ N^(1+ε) ; équivaut à un plancher ρ > 1/2 dans l'algorithme de Claude | Très difficile (circuits non commutatifs) | N1, G5b |
| ❌ **Tests de parité économiques (Q2)** — fausse sur la famille de référence (G8b, lemme CA) | Tests avec O(n·2^O(k)) incidences au total dans l'annexe D de Claude | Technique | G6c |
| 🆕 **Erreurs par point, cas non commutatif** | Étendre le lemme de G6d au-delà de la classe commutative | Moyenne | G6d, G7a |
| 🆕 **Forme parité et contraction de tenseurs (Limits)** | Analyser les limites de ces deux voies d'OpenAI, non couvertes par Limits ni par notre note | Difficile | G4c, G6a |
| 🆕 **Sparsification avec variables d'extension** | Échappatoire à la proposition L de G6c | Moyenne | G6c, G7d |
| 🆕 **Égalité « moyenne sur les cartes témoins = choix indépendants » en Lean** | Comptage sur un produit | Facile | G9L |
| 🆕 **Lemme 2 en Lean** | Anti-concentration avec constante explicite | Difficile (probabilités) | G9L |
| 🆕 **Lemme 4 en Lean** | Graphes aléatoires | Très difficile (peu de choses dans Mathlib) | G9L |

| 🆕 **Formes closes aux niveaux profonds (suite de G8c)** | La valeur des sous-appels de l'algorithme de Claude sur une entrée « épinglée » a une forme close calculable en temps polynomial au niveau 1 (vérifié). Si cela reste vrai aux niveaux profonds, Extended Frege réfuterait toute 3-CNF creuse insatisfaisable en taille 2^(δn). Vérifier aussi l'hypothèse (S) : degré ≤ 1 en chaque table de poids | Moyenne à difficile | G8c |


| 🆕 **Dépasser n^(2/3) dans le modèle produit** | La méthode du déterminant témoin est limitée à n^(2/3) (prop. S de G8a). Il faut une contrainte nouvelle qui couple les points ou qui utilise la structure polynomiale. Ou construire un schéma de coût n^(2/3) : il a besoin de Ω(n^(2/3)) points | Difficile | G8a |
| 🆕 **Borne de Rado pour le modèle DP corrigé** | nullité ≥ max_S(\|S\| − dim Σ U_ω(v)) contient le lemme 1 et le lemme H ; pourrait enlever l'hypothèse de régularité | Moyenne | G8a |
| 🆕 **Relecture adverse de G8a** (théorème B2, prop. S, lemme H) | Niveau P pour l'instant | Facile à lancer | G8a |

| 🆕 **Compromis (I′) et (II) pour la méthode d'OpenAI** | (I′) ρ(Cδ², δ) ≥ δ^(−a) ; (II) ρ(ℓ₁+ℓ₂, d₁+d₂) ≥ ρ(ℓ₁,d₁)·ρ(ℓ₂,d₂). Ensemble : √log n optimal à (log n)^o(1) près | Difficile | G8d |

| 🆕 **Relecture et nouveauté de G9e** | Vérifier la preuve (encodage CNF, identité Σ p^cc) et lire les articles voisins (arXiv 2611.00994, 2611.02548, 2611.02650, dépôt eth-hamiltonicity) pour savoir si Steiner / plus long chemin sont déjà revendiqués | Facile à moyenne | G9e |

| 🆕 **Pièces hors du cadre « un seul espace gradué »** | Seule façon restante d'obtenir ρ < 1/√2 : des pièces qui ne passent pas par les composantes graduées d'un seul espace (ou des profils de croissance non concaves). Lié au problème de Hurwitz : P(N) ≤ N^(1+ε) | Très difficile | G9a |

| 🆕 **Termes non produits (« intriqués »)** | Seule hypothèse non testée du mur 2^(m/log m). Il faut E1 : décomposition jointe d'un tenseur de grille 4×h (termes colonnes ≤ (1/2+o(1))·S bits, termes lignes de rang 2^o(S)) et E2 : évaluation récursive des parties fusionnées. Aucune borne connue n'exclut 2^(S^(1−ε)) dans la classe générale | Très difficile ; c'est la question centrale pour la marche 3 | G9d |

| 🆕 **Pièces sur un corps fini fixe avec ρ < 1/√2** | Seule liberté propre à la décision. Sur GF(2) : meilleur connu 0,811 (classe restreinte) ; recherche exhaustive à λ = 2 non terminée | Difficile | G9b |

## 3. En attente

| Piste | Raison de l'attente |
|---|---|
| ⏸️ **Relecture de l'article par un expert humain** | Auteurs naturels à contacter : Dias et Lima. Décision de l'auteur. |
| ⏸️ **Points nouveaux de la version 2, non relus** | Note de bas de page sur le théorème 7.1(ii) de Dias-Lima ; annexe A (preuves des lemmes A.1/A.2 de Dias-Lima reprises) ; références classiques citées de mémoire |
| ⏸️ **sep(F\*) ≥ c·n\*** | Laissé comme problème ouvert dans l'article |

| Piste | Raison de l'attente |
|---|---|
| ⏸️ **Relecture adverse de G7b** (modèle DP trop généreux, prop. O1) | Le rapport G8R confirme déjà la prop. O1 sur le texte source ; relecture complète à faire avec G8a |

## 4. Terminées (résultats)

| Piste | Résultat |
|---|---|
| ✅ **Article, version 2** (G9W) | `resultats/article/` : 26 pages, LaTeX + PDF ; les 10 corrections obligatoires du rapport G8R sont traitées ; hypothèse M ≥ n supprimée sur papier (Lean garde n ≤ M) |
| ✅ Limits §9, modèle produit (G6a) | Réponse négative, 2^Ω(n^(2/3)) ; confirmée par relecture adverse ; note rédigée |
| ✅ Auteur humain de la note | **Lucas Requena**, étudiant en informatique (« Computer Science student »), contact : l.requena@outlook.fr (avec une mention claire du travail fait par les agents IA) |
| ✅ Lemme 1 en Lean (G8L) | Prouvé sans `sorry` |
| ✅ Partie déterministe en Lean (G9L) | Schéma, lemme 1 complet, fait R, lemme 3, théorème A (sans condition), théorème B (avec le lemme 2 en hypothèse) ; sans `sorry` |
| ✅ Hurwitz P(N) ≤ N^√2 (G4a) | Correct, mais déjà publié par d'autres |
| ✅ Réfutations Extended Frege (G8c) | **Conditionnel :** sous (S) et des preuves EF courtes de « Pin_F(b) », EF réfute F en taille 2s + poly(n)·2^(δn). Le calcul de l'algorithme est de **profondeur constante** (nouveau) : le certificat IPS peut être une formule de profondeur constante. Pas de route vers TC⁰-Frege (descente séquentielle). |
| ✅ Conséquences de la chute de l'ETH (G9e) | Arbre de Steiner, plus long cycle/chemin sur graphes creux en 2^(O(m log log m/√log m)) (conditionnel, non relu, nouveauté non établie) |
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

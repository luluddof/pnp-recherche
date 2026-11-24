# Pistes de recherche

Ce fichier liste toutes les pistes : celles en cours, celles pas encore explorées et celles fermées.
Mise à jour : 24 novembre 2026.

**Rappel :** aucune de ces pistes ne mène directement à P = NP. Ce sont des étapes qui peuvent donner des résultats vérifiables.

Légende : 🔄 en cours · 🆕 pas encore explorée · ⏸️ en attente · ✅ terminée · ❌ fermée (mur prouvé ou mesuré)

---

## 1. En cours (génération 8)

| Id | Piste | Question précise | Origine |
|---|---|---|---|
| G8d | 🔄 **Turán sur l'hypercube** | Prouver ρ(∞, d) ≥ poly(1/d) dans la méthode d'OpenAI (√log n optimal) | G7a |

## 2. Pas encore explorées

| Piste | Question précise | Difficulté | Origine |
|---|---|---|---|
| 🆕 **Lemme 4 pondéré pour F\*** | Prouver la borne de séparateur pour la formule F\* de la prépublication et pour les feuilles de la sparsification, pour comparer avec le théorème 8.3 de Dias-Lima | Moyenne | Rapport G8R |
| 🆕 **Problème de Naslund / Hurwitz** | Prouver P(N) ≥ N^(1+ε) ; équivaut à un plancher ρ > 1/2 dans l'algorithme de Claude | Très difficile (circuits non commutatifs) | N1, G5b |
| 🆕 **Gadgets à 3 parties ou plus** | Seule voie pour passer 2^(m/log m) dans le cadre de Claude | Difficile, aucune idée concrète | G4b |
| 🆕 **Récursion sur des instances non cliques** | Autre façon de sortir du modèle à deux partitions | Difficile | G4b |
| ❌ **Tests de parité économiques (Q2)** — fausse sur la famille de référence (G8b, lemme CA) | Tests avec O(n·2^O(k)) incidences au total dans l'annexe D de Claude | Technique | G6c |
| 🆕 **Erreurs par point, cas non commutatif** | Étendre le lemme de G6d au-delà de la classe commutative | Moyenne | G6d, G7a |
| 🆕 **Forme parité et contraction de tenseurs (Limits)** | Analyser les limites de ces deux voies d'OpenAI, non couvertes par Limits ni par notre note | Difficile | G4c, G6a |
| 🆕 **Sparsification avec variables d'extension** | Échappatoire à la proposition L de G6c | Moyenne | G6c, G7d |
| 🆕 **Égalité « moyenne sur les cartes témoins = choix indépendants » en Lean** | Comptage sur un produit | Facile | G9L |
| 🆕 **Lemme 2 en Lean** | Anti-concentration avec constante explicite | Difficile (probabilités) | G9L |
| 🆕 **Lemme 4 en Lean** | Graphes aléatoires | Très difficile (peu de choses dans Mathlib) | G9L |

| 🆕 **Formes closes aux niveaux profonds (suite de G8c)** | La valeur des sous-appels de l'algorithme de Claude sur une entrée « épinglée » a une forme close calculable en temps polynomial au niveau 1 (vérifié). Si cela reste vrai aux niveaux profonds, Extended Frege réfuterait toute 3-CNF creuse insatisfaisable en taille 2^(δn). Vérifier aussi l'hypothèse (S) : degré ≤ 1 en chaque table de poids | Moyenne à difficile | G8c |

| 🆕 **Q1\* (suite de G8b)** | Compter les ensembles indépendants de tout graphe à n·polylog(n) arêtes en temps 2^O(n/√log n) ; enlèverait le facteur log log n de l'algorithme de Claude | Difficile ; aucune borne inférieure connue ne bloque | G8b |

| 🆕 **Dépasser n^(2/3) dans le modèle produit** | La méthode du déterminant témoin est limitée à n^(2/3) (prop. S de G8a). Il faut une contrainte nouvelle qui couple les points ou qui utilise la structure polynomiale. Ou construire un schéma de coût n^(2/3) : il a besoin de Ω(n^(2/3)) points | Difficile | G8a |
| 🆕 **Borne de Rado pour le modèle DP corrigé** | nullité ≥ max_S(\|S\| − dim Σ U_ω(v)) contient le lemme 1 et le lemme H ; pourrait enlever l'hypothèse de régularité | Moyenne | G8a |
| 🆕 **Relecture adverse de G8a** (théorème B2, prop. S, lemme H) | Niveau P pour l'instant | Facile à lancer | G8a |

## 3. En attente

| Piste | Raison de l'attente |
|---|---|
| ⏸️ **Réécriture de la note** selon le rapport G8R (révision majeure) | Attendre G8a pour intégrer son résultat (G9L est terminé) |
| ⏸️ **Relecture adverse de G7b** (modèle DP trop généreux, prop. O1) | Le rapport G8R confirme déjà la prop. O1 sur le texte source ; relecture complète à faire avec G8a |

## 4. Terminées (résultats)

| Piste | Résultat |
|---|---|
| ✅ Limits §9, modèle produit (G6a) | Réponse négative, 2^Ω(n^(2/3)) ; confirmée par relecture adverse ; note rédigée |
| ✅ Auteur humain de la note | **Lucas Requena**, étudiant en informatique (« Computer Science student »), contact : l.requena@outlook.fr (avec une mention claire du travail fait par les agents IA) |
| ✅ Lemme 1 en Lean (G8L) | Prouvé sans `sorry` |
| ✅ Partie déterministe en Lean (G9L) | Schéma, lemme 1 complet, fait R, lemme 3, théorème A (sans condition), théorème B (avec le lemme 2 en hypothèse) ; sans `sorry` |
| ✅ Hurwitz P(N) ≤ N^√2 (G4a) | Correct, mais déjà publié par d'autres |
| ✅ Réfutations Extended Frege (G8c) | **Conditionnel :** sous (S) et des preuves EF courtes de « Pin_F(b) », EF réfute F en taille 2s + poly(n)·2^(δn). Le calcul de l'algorithme est de **profondeur constante** (nouveau) : le certificat IPS peut être une formule de profondeur constante. Pas de route vers TC⁰-Frege (descente séquentielle). |
| ✅ Système de preuve IPS (G4d) | Réfutations IPS de taille 2^O(n/√log n) |
| ✅ Plancher ρ_d ≥ 1/d (G6b) | Confirmé sous conditions |

## 5. Fermées (mur prouvé ou mesuré)

| Piste | Cause |
|---|---|
| ❌ Algorithmes locaux, LP, SOS, spectral, PC, plans de coupe, parités courtes, variables d'extension simples | Bornes inférieures connues ou croissance exponentielle mesurée (générations 0 à 3) |
| ❌ Anneaux Z_m, Z, Q | Réduction modulo p (preuve) |
| ❌ Certificats de Feige-Kim-Ofek trouvables en temps polynomial | Égalent seulement le seuil spectral n^1,5 |
| ❌ Géométries de récursion (grilles, boîtes de dimension d) | 2^Ω(m/log m) dans le modèle (G4b, G6b) |
| ❌ Structures non produit pour OpenAI (expanseurs, codes) | Rang ≥ K/2 (G5a) |
| ❌ Combinaisons OpenAI + Claude | Même mur 2^(n/log n) (G4c) |
| ❌ Réfuter Feige (m = O(n)) | Bornes inférieures SOS ; hors de portée |

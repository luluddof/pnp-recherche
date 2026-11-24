# Pistes de recherche

Ce fichier liste toutes les pistes : celles en cours, celles pas encore explorées et celles fermées.
Mise à jour : 24 novembre 2026.

**Rappel :** aucune de ces pistes ne mène directement à P = NP. Ce sont des étapes qui peuvent donner des résultats vérifiables.

Légende : 🔄 en cours · 🆕 pas encore explorée · ⏸️ en attente · ✅ terminée · ❌ fermée (mur prouvé ou mesuré)

---

## 1. En cours (génération 8)

| Id | Piste | Question précise | Origine |
|---|---|---|---|
| G8a | 🔄 **Terminer Limits §9** | (a) Modèle DP corrigé : enlever l'hypothèse de régularité. (b) La borne n^(2/3) du modèle produit est-elle optimale ? | G6a, G7b |
| G8b | 🔄 **Compter les ensembles indépendants** (Q1′) | Compter les ensembles indépendants d'un graphe de degré moyen d en temps 2^(O(n/d)) | G7d |
| G8c | 🔄 **Réfutations dans un système faible** | 3-SAT aléatoire a-t-il des réfutations de taille 2^o(n) en Extended Frege ou TC⁰-Frege ? | G4d |
| G8d | 🔄 **Turán sur l'hypercube** | Prouver ρ(∞, d) ≥ poly(1/d) dans la méthode d'OpenAI (√log n optimal) | G7a |

## 2. Pas encore explorées

| Piste | Question précise | Difficulté | Origine |
|---|---|---|---|
| 🆕 **Lemme 4 pondéré pour F\*** | Prouver la borne de séparateur pour la formule F\* de la prépublication et pour les feuilles de la sparsification, pour comparer avec le théorème 8.3 de Dias-Lima | Moyenne | Rapport G8R |
| 🆕 **Problème de Naslund / Hurwitz** | Prouver P(N) ≥ N^(1+ε) ; équivaut à un plancher ρ > 1/2 dans l'algorithme de Claude | Très difficile (circuits non commutatifs) | N1, G5b |
| 🆕 **Gadgets à 3 parties ou plus** | Seule voie pour passer 2^(m/log m) dans le cadre de Claude | Difficile, aucune idée concrète | G4b |
| 🆕 **Récursion sur des instances non cliques** | Autre façon de sortir du modèle à deux partitions | Difficile | G4b |
| 🆕 **Tests de parité économiques (Q2)** | Tests avec O(n·2^O(k)) incidences au total dans l'annexe D de Claude | Technique | G6c |
| 🆕 **Erreurs par point, cas non commutatif** | Étendre le lemme de G6d au-delà de la classe commutative | Moyenne | G6d, G7a |
| 🆕 **Forme parité et contraction de tenseurs (Limits)** | Analyser les limites de ces deux voies d'OpenAI, non couvertes par Limits ni par notre note | Difficile | G4c, G6a |
| 🆕 **Sparsification avec variables d'extension** | Échappatoire à la proposition L de G6c | Moyenne | G6c, G7d |
| 🆕 **Égalité « moyenne sur les cartes témoins = choix indépendants » en Lean** | Comptage sur un produit | Facile | G9L |
| 🆕 **Lemme 2 en Lean** | Anti-concentration avec constante explicite | Difficile (probabilités) | G9L |
| 🆕 **Lemme 4 en Lean** | Graphes aléatoires | Très difficile (peu de choses dans Mathlib) | G9L |

## 3. En attente

| Piste | Raison de l'attente |
|---|---|
| ⏸️ **Réécriture de la note** selon le rapport G8R (révision majeure) | Attendre G8a pour intégrer son résultat (G9L est terminé) |
| ⏸️ **Relecture adverse de G7b** (modèle DP trop généreux, prop. O1) | Le rapport G8R confirme déjà la prop. O1 sur le texte source ; relecture complète à faire avec G8a |

## 4. Terminées (résultats)

| Piste | Résultat |
|---|---|
| ✅ Limits §9, modèle produit (G6a) | Réponse négative, 2^Ω(n^(2/3)) ; confirmée par relecture adverse ; note rédigée |
| ✅ Auteur humain de la note | **Lucas Requena** (avec une mention claire du travail fait par les agents IA) |
| ✅ Lemme 1 en Lean (G8L) | Prouvé sans `sorry` |
| ✅ Partie déterministe en Lean (G9L) | Schéma, lemme 1 complet, fait R, lemme 3, théorème A (sans condition), théorème B (avec le lemme 2 en hypothèse) ; sans `sorry` |
| ✅ Hurwitz P(N) ≤ N^√2 (G4a) | Correct, mais déjà publié par d'autres |
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

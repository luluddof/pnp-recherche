# Journal de recherche : P vs NP

Date de début : 23 novembre 2026.
Objectif choisi par l'utilisateur : chercher **P = NP**, c'est-à-dire un algorithme polynomial pour 3-SAT, avec une preuve.

Ce journal contient toutes les voies explorées, les résultats et la cause de chaque échec.
Les programmes et les données sont dans le dossier `code/`.

> **Rappel honnête.** Aucune voie ci-dessous ne résout P vs NP. Chaque résultat est soit un échec mesuré, soit un échec expliqué par un théorème connu.

---

## ⚠️ Niveau de fiabilité de chaque résultat (à lire avant tout le reste)

Ce journal est un **résumé**. Un résumé ne permet pas de vérifier une preuve. Le tableau ci-dessous dit, pour chaque résultat, ce qui a réellement été vérifié.

**Échelle :**
- **L** = prouvé en Lean (vrai pour **tout** n ; vérification indépendante de nos raisonnements).
- **P+** = preuve complète écrite **et** refaite par un relecteur. Attention : nos relecteurs sont des agents Claude, donc **pas indépendants** (même modèle, souvent mêmes définitions).
- **P** = preuve écrite dans des notes de travail, **non relue**.
- **E** = esquisse ou argument incomplet.
- **N** = seulement des tests numériques. Un test sur n = 12 à 30 ou ℓ = 2 à 10 **détecte des erreurs, mais ne prouve rien** pour les grandes valeurs.
- Aucun résultat n'a encore été relu par un **expert humain**.

| Résultat | Énoncé (portée exacte) | Niveau | Où est la preuve |
|---|---|---|---|
| G6a lemme 1, fait R, lemme 3, théorème A | Énoncés déterministes de la note | **L** + P+ | `lean/PnpLean/*.lean`, note §3 |
| G6a théorème B | **Sans hypothèse** depuis G10e | **L** + P+ | `lean/PnpLean/Article1/Lemma2.lean` |
| G6a lemme 2 (anti-concentration, espérance du déficit) | Constante c₀ = (1−e⁻²)/4 | **L** + P+ | `lean/PnpLean/Article1/Anticonc.lean`, `Lemma2.lean` |
| G6a lemme 4 (3-SAT aléatoire) | sep(F) ≥ cn avec forte probabilité | P+ | note §3.7 |
| G6a corollaire C | **Seulement dans le modèle de coût produit**, pour F lui-même (pas F\*) ; rien sur la programmation dynamique, les tenseurs ou les parités | P+ | note §3.8 |
| G7b prop. O1 | Le modèle DP de Dias-Lima est dégénéré | P+ (confirmé par G8R sur le texte source) | `code/gen7/g7b/notes.md` |
| G8a prop. S | Résultat principal correct ; 2 constantes à corriger | P+ (en partie) | `code/gen10/g10a/` |
| G7b théorèmes D, F, lemme E | Confirmés par G10a. **Théorème C : clause « de plus » réfutée (erreur de 1).** **Le modèle DP corrigé ne minore pas tous les évaluateurs** (seulement l'élimination par groupes entiers) | P+ (en partie) | idem |
| G4a Hurwitz P(N) ≤ N^√2 | Correct, **déjà publié par d'autres** | P+ | `code/gen4/g4a/proof.md` |
| G6b ρ_d ≥ 1/d | Dans un modèle restreint, sous 2 conditions | P+ | rapport G6r |
| G4d réfutations IPS 2^O(n/√log n) | Corollaire de l'article de Claude (non relu par des pairs) | E | `code/gen4/g4d/` |
| G6c prop. R, prop. L | Prop. L seulement pour une classe restreinte de sparsifieurs | P | `code/gen6/g6c/` |
| G6d lemme du déterminant par bloc (général) | Confirmé par G10b ; la conséquence (c) n'est pas prouvée hors de la classe commutative | **L** + P+ | `lean/PnpLean/Article2/PerPoint.lean` |
| G6d conséquences (a), (b) | (b) dépend du modèle de coût | P+ | `code/gen6/g6d/notes.md` |
| G7a famille « par morceaux », barrière A | Classe commutative seulement | P (+ N pour ℓ ≤ 10) | `code/gen7/g7a/notes.md` |
| G7d lemme S, prop. E | Classes restreintes de réductions | E | `code/gen7/g7d/` |
| G8a théorème B2, lemme H | Confirmés par G10a | **P+** | `code/gen8/g8a/notes.md` |
| G8d lemme Q_{d+2} | Classe commutative, un point par direction ; confirmé par G10b | **L** + P+ | `lean/PnpLean/Article2/QdPlus2.lean` |
| G8d corollaire ρ ≥ 2^(−d)/(d+2) | idem | P+ | `code/gen8/g8d/notes.md` |
| G9e conséquences de la chute de l'ETH | **Conditionnel** aux prépublications ; confirmé par G10c ; **Steiner déjà revendiqué par Wu-Wang (arXiv 2611.02548)** ; plus long chemin/cycle : pas de revendication trouvée | P+ | `code/gen9/g9e/consequences.md` |
| G9a obstruction A (α ≥ 1/2 dans le cadre de l'annexe B) | Pièces qui passent par les composantes graduées d'un seul espace | P (+ N) | `code/gen9/g9a/notes.md` |
| G9d loi harmonique, analyse des triangles | Dans un modèle de coût ; termes produits seulement | E | `code/gen9/g9d/g9d_model.py` |
| G9b ρ = 0,811 sur GF(2), lemme de la trace | Petite classe de pièces (puissances tensorielles à un bit) | P (+ N) | `code/gen9/g9b/notes.md` |
| G9c théorèmes A et B (#IS sans log d sur certaines familles) | Confirmés par G10d, **relativement à l'annexe D de Claude (non relue)** ; gain modeste (facteur polylog(d)) | P+ | `code/gen9/g9c/g9c_analysis.md` |
| G10f sep(F\*) ≥ c·n\* et corollaire pour F\* | Δ ≥ 2 ; confirmé par G11a ; **hypothèse : F\* a une copie par occurrence** (inférence, pas une citation) | P+ | `code/gen10/g10f/fstar_separator.md` |
| **G11b théorème Σ rang(E_k) ≥ N_m** | Cadre des erreurs par point, un point par direction, sur ℚ | **L** | `lean/PnpLean/Article2/SqrtLog.lean` |
| G11b conséquences (seuil sur d, optimalité de √log n, temps de calcul) | Asymptotique ; modèle de coût de G6d seulement ; nouveauté non vérifiée | P (G13a en cours) | `code/gen11/g11b/notes.md` |
| G11c lemme Q_{d+2} général, formule r_k = max(rang B01, rang B10) | Un point par direction ; preuves à la main ; les valeurs « trouvées » sont numériques | P (+ N) | `code/gen11/g11c/notes.md` |
| G8b lemme CA, échec de Q2 | Tests construits sur les clauses d'arêtes seulement | P (+ N) | `code/gen8/g8b/` |
| G8c réfutations EF | **Conditionnel** à l'hypothèse (S) non vérifiée | E | `code/gen8/g8c/` |
| G4b, G4c, G5a (murs des méthodes de 2026) | Chaque mur vaut **seulement dans le modèle étudié** | E à P | `code/gen4/`, `code/gen5/` |

**Ce qu'il faut faire pour qu'un résultat devienne fiable :** une preuve complète dans un document autonome (pas un renvoi vers du code), puis une relecture par un expert humain. Aujourd'hui, seule la note sur G6a remplit la première condition.

---

## Bilan après 5 générations (23 novembre 2026)

- **P vs NP n'est pas résolu.** Aucune voie ne mène à un algorithme polynomial pour 3-SAT.
- **20 voies explorées** par 26 sous-agents, plus 2 outils (banc de test, cahier des charges) et 1 relecture adverse.
- **L'obstacle central :** prouver qu'une formule 3-SAT aléatoire est insatisfaisable avec m = O(n) clauses. Toutes les méthodes testées demandent m ≈ n^1,5, ou un temps exponentiel.
- **Le mur des algorithmes de 2026 :** les algorithmes d'OpenAI et de Claude (2^(n/√log n)) utilisent tous les deux une structure produit. Dans ce cadre, on ne descend pas sous environ 2^(n/log n) (G4b, G4c, G5a).
- **Seul résultat mathématique « nouveau » trouvé :** la borne P(N) ≤ N^√2 sur une variante du problème de Hurwitz. Elle est correcte, mais **déjà publiée** par d'autres en novembre 2026.
- **Questions ouvertes les plus concrètes :**
  1. Le schéma « produit de clauses » avec √n groupes sur un expanseur (papier Limits, §9). Il donnerait 2^Õ(√n), toujours pas polynomial.
  2. Une identité à d directions pour les diagonales d'une boîte de dimension d, avec ρ_d < ρ₂^(log₂ d) (G4b).
  3. Des erreurs différentes en chaque point dans la méthode d'OpenAI (G4c, G5a).

---

## 0. Faits de base à retenir

| Fait | Source |
|---|---|
| L'algorithme universel de Levin (1973) trouve les solutions en temps polynomial **si** P = NP. Le code n'est donc pas la difficulté. La difficulté est la **preuve**. | Levin 1973 |
| Pour prouver P = NP, l'algorithme doit être correct sur **toutes** les formules. Les tests trouvent seulement des contre-exemples. | — |
| **Le vrai obstacle mesuré :** prouver qu'une formule 3-SAT aléatoire (densité 4,26) est UNSAT. Trouver des solutions est plus facile en pratique. | Sous-agent 1 |
| Réfuter 3-SAT aléatoire avec m = O(n) clauses réfuterait l'hypothèse de Feige. | Feige 2002 |
| Méthode spectrale : elle réfute seulement si m ≥ environ 7,4 × n^1,5 (mesuré pour n = 30 à 150). | Sous-agent 1 |
| Cahier des charges complet (13 conditions et 14 familles éliminées) : `code/spec/requirements.md` | Étape 2 |

---

## 1. Génération 0 (moi, sans sous-agent)

### 1.1 Algorithme à regard en avant de profondeur k (direction P = NP)
- **Idée :** propagation unitaire + test des littéraux qui échouent jusqu'à la profondeur k + choix glouton, sans retour en arrière. Temps en O(n^(k+2)).
- **Résultat :**

| k | n | SAT trouvées | UNSAT prouvées |
|---|---|---|---|
| 1 | 40 | 15/28 | 0/12 |
| 1 | 120 | 2/25 | 0/15 |
| 2 | 40 | 22/28 | 1/12 |
| 2 | 120 | 5/25 | 0/15 |
| 3 | 30 | 21/21 | 9/9 |
| 3 | 80 | 7/13 | 0/17 |

- **Cause de l'échec (théorème) :** une réponse UNSAT de profondeur k est une preuve par résolution de largeur ≈ k + 3. Ben-Sasson et Wigderson (2001) : la largeur nécessaire est Ω(n).
- **Fichiers :** `code/gen0/look.c`, `code/gen0/results3.txt`

### 1.2 Géométrie de l'ensemble des solutions (physique statistique, topologie)
- **Hypothèse :** beaucoup d'« îles » de solutions ⇒ problème difficile.
- **Résultat :** 3-XORSAT (dans P) a 416 îles en moyenne (α = 0,7), alors que 3-SAT (NP-complet) en a 4 (α = 4,2).
- **Conclusion :** **hypothèse fausse.** La géométrie seule ne mesure pas la difficulté.
- **Fichier :** `code/gen0/geo.c`

### 1.3 Polymorphismes (algèbre universelle)
- **Résultat :** 3-SAT et 1-in-3 SAT n'ont aucun polymorphisme non trivial. 2-SAT a la majorité, XORSAT a x⊕y⊕z, Horn a ET. Cela confirme le théorème de Schaefer.
- **Pourquoi l'idée ne suffit pas :** « aucun polymorphisme ⇒ NP-complet » ne donne pas « ⇒ pas dans P ». L'argument est circulaire.
- **Fichier :** `code/gen0/poly2.py` (`poly.py` contient une première version avec une erreur : une seule contrainte au lieu de toutes)

---

## 2. Génération 1 : 5 sous-agents en parallèle (exploration large)

| # | Direction | Résultat | Cause de l'échec | Dossier |
|---|---|---|---|---|
| 1 | Algorithmes polynomiaux plus forts (P = NP) | 0 % des formules UNSAT aléatoires prouvées dès n = 50, même avec toutes les méthodes ensemble | Polynomial Calculus : degré Ω(n). LP/SOS : degré Ω(n) (Grigoriev, Schoenebeck). Méthode spectrale : seuil n^1,5. | `code/agent1/` |
| 2 | Bornes inférieures sur les circuits | Tailles exactes pour n = 4 identiques à la table de Knuth. 10 mesures testées, meilleure corrélation 0,64. | Mesures limitées à O(n), ou fausses (AND et MAJ), ou bloquées par les preuves naturelles | `code/agent2/` |
| 3 | Complexité des preuves (AC0[p]-Frege) | Résolution exponentielle sur PHP, Tseitin et formules aléatoires | Il manque un objet avec un ordre (ℝ) **et** une approximation des portes MOD (F_p) | `code/agent3/` |
| 4 | Méta-complexité (MCSP, magnification) | Écart à combler : 3,1 N contre N^(1+ε) | Barrière de localité (CHOPRS 2020) | `code/agent4/` |
| 5 | Algèbre, relaxations BLP/AIP, GCT | Les relaxations décident 2-SAT, XORSAT et Horn, mais pas 3-SAT | Bornes seulement pour des modèles limités. Rothvoss 2014 : le couplage est dans P sans petit LP. GCT : obstructions d'occurrence impossibles (BIP 2016). | `code/agent5/` |

---

## 3. Génération 2 : construction par briques simples (direction P = NP)

### Étape 1 – Banc de test commun ✅
- **Outil :** `code/bench/bench.py`
- **Familles de formules :** petites formules (n = 4 à 10), 3-SAT aléatoire, 3-SAT avec solution cachée, Tseitin SAT et UNSAT, tiroirs (PHP), 3-XOR SAT, UNSAT et avec solution cachée.
- **Vraies réponses vérifiées de 3 manières :** recherche exhaustive, CaDiCaL, construction connue.
- **Minimiseur :** il trouve le plus petit contre-exemple d'un candidat.
- **Utilisation :**
  ```bash
  cd code/bench && python3 bench.py --candidate cand.py:solve --suite quick --minimize
  ```
- **Remarque :** le fichier `truth_cache.json` est la mémoire des vraies réponses. Il évite de recalculer.

### Étape 2 – Cahier des charges ✅
- **Fichier :** `code/spec/requirements.md`
- **5 directions qu'aucun théorème connu n'élimine :**
  1. Variables d'extension (Extended Frege, DRAT/PR)
  2. Branchement sur les parités, Res(⊕)
  3. Plans de coupe sur 3-CNF aléatoire
  4. Algèbre sur des anneaux qui ne sont pas des corps
  5. Certificats de Feige-Kim-Ofek (m ≈ n^1,4) — **pas encore explorée**

### Étape 3 – Tests des directions

| Étape | Direction | Résultat | Détail | Dossier |
|---|---|---|---|---|
| 3a | Variables d'extension (XOR, AND, OR, BVA, aléatoire) | ❌ Exponentiel, environ 2^(n/21) pour n = 30 à 200 | Toutes les variantes restent à ±15 % de la base. Le contrôle aléatoire est aussi bon que les variantes « intelligentes ». | `code/step3a/` |
| 3b | Branchement sur les parités (DPLL + élimination GF(2)) | ❌ Même base exponentielle que DPLL normal (≈ 1,04^n). Les heuristiques « parité pure » sont pires (1,12^n à 1,30^n). | Contrôle Tseitin : le branchement sur la parité des clauses passe de 2^(n/2) à ≈ 2n feuilles, donc l'outil fonctionne. Les parités aident seulement s'il existe une structure XOR, et 3-SAT aléatoire n'en a pas. Limite : n ≤ 60. | `code/step3b/` |
| 3c | Plans de coupe (CBC, HiGHS, Gomory) | ❌ Exponentiel, environ e^(0,1 n), plus lent que CDCL | La relaxation LP est toujours satisfaite avec toutes les variables à 1/2. Contrôle : PHP résolu immédiatement, donc l'outil fonctionne. | `code/step3c/` |
| 3d | Anneaux Z₄, Z₆, Z_m | ❌ Éliminé par une preuve | Réduction modulo p : deg(Z_m) ≥ deg(GF(p)) pour tout p qui divise m. Expérience : même degré pour tous les anneaux sur chaque formule. | `code/step3d/` |

**Limites importantes à retenir :**
- 3a et 3c mesurent seulement les preuves que les **solveurs trouvent**. Des preuves plus courtes peuvent exister.
- 3d n'a pas testé l'anneau des entiers Z.

---

## 4. Génération 3 (en cours) : nouvelles voies inspirées par l'actualité

| Voie | Question | État |
|---|---|---|
| N1 | Facteur ρ de l'algorithme de Claude (tâche corrigée après lecture de l'article) | ❌ pour P = NP. Détails ci-dessous. Dossier : `code/gen3/n1/` |

**Détails de N1 :**
- Identité des huit carrés vérifiée sur les 4 096 cas : 0 erreur.
- Recherche exhaustive pour λ = 1 : meilleur ρ = 0,896 (identité à un bit).
- **Plancher prouvé : ρ ≥ 1/2** pour toute décomposition (car P(N) ≤ #pièces × N^(2ρ) et P(N) ≥ N). Atteindre ρ → 1/2 demande des compositions de type Hurwitz presque optimales. Prouver ρ > 1/2 est au moins aussi difficile que le problème ouvert de Naslund.
- **Limite de la méthode :** temps 2^(m/(log m)^(log₂(1/ρ))), donc au mieux 2^(O(m/log m)). Pour atteindre 2^(n^(1−ε)), il faudrait environ log m niveaux, mais la grille élève K au carré à chaque niveau. Il faut changer la géométrie de la récursion, pas seulement les pièces. Le temps quasi-polynomial est hors de portée de cette méthode.
- Inclusion-exclusion simple tronquée : la troncature doit être t ≈ 1,26 n, donc environ 2^(3,6 n) termes. Inutile seule.
- 🔎 **Observation à vérifier (peut-être nouvelle) :** l'inégalité P(N) ≤ #pièces × N^(2ρ), appliquée aux pièces de l'annexe B (ρ = 1/√2), donne **P(N) ≤ N^(√2 + o(1)) ≈ N^1,414**. Cela améliorerait la borne N^1,62 de Hrubeš (2024) pour une variante du problème de Hurwitz. L'article donne seulement P(N) = O(N^(11/6)) pour les huit carrés. **Non vérifié** : il faut contrôler que les pièces de l'annexe B satisfont exactement la définition de P(N) (fonctions rationnelles, somme exacte, nombre de pièces polynomial), et que ce résultat n'est pas déjà connu.
| N2 | Certificats courts de Feige-Kim-Ofek (m ≈ n^1,4) : peut-on les trouver en temps polynomial et battre le seuil spectral ? | ❌ **Non.** Les « couvertures paires » trouvables en temps polynomial sont des cycles de longueur O(log n). Elles existent en nombre suffisant seulement si m = Ω(n^1,5) : au mieux, on égale la méthode spectrale. Sous n^1,5, il faut des couvertures de taille ≈ n^0,2 (dépendances creuses sur GF(2)) : aucune méthode polynomiale connue. Kikuchi : temps n^O(r), donc sous-exponentiel pour r = n^0,2. À n = 3000, aucun LP de couvertures ne peut réussir. Obstacle théorique : SOS de degré r échoue si m ≲ n^1,5/√r (Kothari-Mori-O'Donnell-Witmer). Dossier : `code/gen3/n2/` |
| N3 | Nullstellensatz sur les entiers Z et les rationnels Q | ❌ **Éliminé.** Preuve : deg(Z) ≥ deg(Q), et deg(Q) = deg(GF(p)) pour presque tout p (un mineur non nul reste non nul modulo presque tout p). Mesure : 130 formules (n = 8 à 20), 0 différence entre GF(2), GF(3), GF(5), GF(1000003), Q et Z. Contrôle Tseitin : l'outil voit bien les différences entre anneaux. Dossier : `code/gen3/n3/` |

## 4 ter. Génération 4 (autonomie, en cours)

Articles disponibles en texte (pour les sous-agents) : article de Claude, exposé d'O'Donnell sur l'article d'OpenAI, article « Limits and Constants ».

| Voie | Question | État |
|---|---|---|
| G4a | Vérifier l'observation P(N) ≤ N^(√2+o(1)) (variante du problème de Hurwitz) | ✅ **Vérifiée** (sous réserve des lemmes 9 et 10 de l'article de Claude). Voir ci-dessous. Dossier : `code/gen4/g4a/` |
| G4b | Autres géométries de récursion pour l'algorithme de Claude | ❌ Aucune géométrie n'atteint 2^(n^(1−ε)). **Lemme 1 (dans le modèle) :** pour réécrire une fraction constante des bits, il faut élever K au carré. **Théorème 1 (esquisse, vérifiée par programmation dynamique) :** dans le modèle à deux partitions, le temps est au moins 2^(Ω(m/log m)). Boîtes de dimension d (spéculatif) : au mieux 2^(m/2^(O(√log m))). Pour aller plus loin, il faudrait des gadgets qui touchent 3 parties ou plus, ou une récursion sur des instances qui ne sont pas des cliques. Dossier : `code/gen4/g4b/` |

**Génération 5 :**
- G5a = structure non produit sur des expanseurs (la piste des « Limits » pour 2^Õ(√n)) : ❌ pour le mécanisme d'OpenAI ; question ouverte pour l'autre.
  - **Question ouverte citée** (Limits, §9) : un schéma sur un expanseur avec √n groupes, largeur O(√n), ordre de Taylor O(√n) et peu de mesures par groupe donnerait 2^Õ(√n). « From below, one must show that measurements cannot absorb the cut. »
  - **Lemme A :** l'erreur de rang minimal est une complétion de rang minimal de l'échange C sur les paires éloignées.
  - **Lemme B :** rang ≥ plus grand anticode de diamètre < D − d ≥ K/(1 + N_d). Cela généralise le lemme de G4c à toute structure d'indices.
  - **Corollaire (prouvé) :** sur un code bon, un expanseur ou un 2-design, le rang est ≥ K/2, sauf si le gain dépasse une constante. L'interpolation demande de la localité, le contraire de l'expansion. Expériences MILP : le cube (structure produit) est toujours le meilleur.
  - **Reste ouvert :** la voie « produit de clauses » avec √n groupes (le problème ouvert du papier Limits), des erreurs différentes en chaque point, et l'optimalité générale du cube.
  - Dossier : `code/gen5/g5a/`.
- G5b = relecture adverse : ✅ terminée.
  - G4a (Hurwitz) : **correct mais déjà publié** (Siemaszko, Zenodo, 7 novembre 2026 ; Debnath, novembre 2026).
  - Lemme de G4c : **correct**, mais c'est une observation courte et élémentaire. Le cas symétrique découle de résultats connus (Aspnes-Beigel-Furst-Rudich, 1994). « √log n est optimal pour la méthode d'OpenAI » est vrai seulement sous des hypothèses précises (erreur E commune, un point par coordonnée, ℓ ≤ log n/log log n).
  - Dossier : `code/gen5/g5b/`.
| G4c | Méthode d'OpenAI : comparaison et combinaison avec celle de Claude | ❌ pour aller plus loin. **Même idée** dans les deux articles : développer un « échange croisé » sur une structure produit (Claude : 2 directions, lignes ou colonnes, ρ = 1/√2 ; OpenAI : ℓ directions, interpolation d'Hermite + erreur de petit rang). Les deux donnent α = 1/2, donc 2^(n/√log n). 3 combinaisons testées : aucune ne change α. **Lemme nouveau (non relu) :** si l'erreur E est la même en tous les points, rang(E) ≥ taille d'une boule de Hamming de rayon (ℓ−d−1)/2. Donc α = 1/2 est optimal pour la méthode d'OpenAI (vérifié numériquement pour ℓ ≤ 100). **Mur commun :** une structure produit d'au plus log n coordonnées ne descend pas sous 2^(n/log n). Piste : une structure **non produit**, par exemple sur des expanseurs avec √n groupes (donnerait 2^Õ(√n)). Dossier : `code/gen4/g4c/` |
| G4d | Quel système de preuve correspond à l'algorithme de comptage ? | ✅ **IPS** (Ideal Proof System, Grochow-Pitassi) sur Q. Toute 3-CNF insatisfaisable a une réfutation IPS de taille 2^(O(m/√log m)) ; pour 3-SAT aléatoire (m = O(n)) : **2^(O(n/√log n))**. Cela sépare quantitativement IPS de PC, Nullstellensatz et SOS (2^Ω(n)) sur ces formules. Aucune voie trouvée vers Frege ou Extended Frege. Corollaire direct, mais aucune borne 2^o(n) antérieure trouvée pour un système standard à densité constante. Pour des réfutations **polynomiales**, il faudrait une hypothèse de type VP = VNP. Expérience : certificat construit et vérifié pour n = 8 à 32. Dossier : `code/gen4/g4d/` |

### Résultat G4a : un corollaire sur le problème de Hurwitz (correct, mais déjà publié par d'autres)

**Ce n'est pas un progrès sur P vs NP.** C'est un résultat annexe, sur une question d'algèbre classique.

- **Énoncé :** P(N) ≤ N^√2 · (log N)^O(1). Conséquence : pour la variante du problème de Hurwitz (Σa_i²)(Σb_j²) = Σ_{k≤s} f_k f'_k (formes bilinéaires sur Q), on a s ≤ N^(√2+o(1)). Sur C, cela donne une formule de composition avec 2·N^(√2+o(1)) carrés.
- **Meilleure borne publiée connue :** O(N^1,62) (Hrubeš, CCC 2024).
- **Preuve (complète dans `code/gen4/g4a/proof.md`) :**
  1. Si Σ A_s = d·[x=v][y=u] avec d ≠ 0, alors P ≤ Σ_s min(rang_lignes(A_s), rang_colonnes(A_s)).
  2. Les 4h+1 pièces de l'annexe B (algèbre de Heisenberg, indexées par les partitions de h) satisfont Σ A_s = h!·[X=V][Y=U] exactement (lemmes 9 et 10).
  3. Chaque pièce passe par des espaces de poids : min(rangs) ≤ Σ_{w≤2h} p(w) ≤ (2h+1)·p(2h).
  4. Hardy-Ramanujan : p(2h) = p(h)^√2 · Θ(h^(√2−1)). Avec M = p(h), on obtient P(M) ≤ M^√2 · polylog(M).
- **Vérifications :** calcul exact en entiers pour h = 2 à 6 (la somme des pièces égale h!·croisement en tous les points ; rangs exacts proches des bornes).
- **Limite :** asymptotique seulement. La borne bat N^1,62 seulement pour N ≥ environ 2^71.
- **Nouveauté : NON.** La relecture adverse (G5b) a trouvé que ce résultat est **déjà publié** : J. P. Siemaszko (7 novembre 2026, Zenodo, DOI 10.5281/zenodo.24024689), repris par L. Debnath (dépôt GitHub puterfaction/hurwitz-sqrt2, novembre 2026). Les deux le tirent aussi du théorème B de l'article de Claude.
- **Relecture adverse (G5b) :** mathématiquement **confirmé** avec un code indépendant (identité exacte pour h = 1 à 6, rangs exacts identiques, identité de composition vérifiée symboliquement pour h = 3 et 4).
- **Portée exacte :** améliore la variante « produit de formes bilinéaires » sur Q (caractéristique 0), et σ_C(N) ≤ 2·N^(√2+o(1)) sur C. N'améliore **pas** : la caractéristique positive fixée (division par h!), les sommes de carrés sur R, les formules entières.

## 4 quater. Génération 6 (en cours) : problèmes ouverts précis, comme étapes

Stratégie : résoudre d'abord des problèmes ouverts plus petits, cités dans les articles de 2026. Chaque résultat positif passe ensuite par une relecture adverse.

| Voie | Problème ouvert | État |
|---|---|---|
| G6a | Limits §9 : schéma « produit de clauses » avec √n groupes sur un expanseur (2^Õ(√n)) — construire ou prouver impossible | ✅ **Confirmé par la relecture adverse G6r** (chaque lemme vérifié ligne par ligne, plus des tests indépendants). Aucune publication antérieure trouvée. **Portée exacte :** contredit le schéma du §9 si les groupes ont une taille s_t = Õ(√n) ; le gain porte sur le **modèle de coût produit** (2^Ω(√n) → 2^Ω(n^(2/3))) ; le modèle de programmation dynamique reste ouvert ; il faut sep(F) ≥ cn (vrai pour 3-SAT aléatoire par une adaptation du lemme 8.2 ; supposé pour les incidences 3-régulières aléatoires). Détails : Lemme 1 (déterminant témoin) : Σ_t (déficience_t − ρ_t)₊ ≤ nR. **Théorème B :** en un point t, ρ_t + q_t + 1 > min(0,034·sep(F)/√s_t, (2n/3 − s_t)/2). **Théorème A** (déformations qui dépendent seulement de la variable) : aucune absorption possible. **Conséquences :** sur un expanseur, le schéma du §9 ne peut pas exister ; dans le modèle de coût « produit », tout schéma coûte au moins 2^Ω(n^(2/3)) (le papier avait 2^Ω(√n)). Reste ouvert : le modèle de programmation dynamique, et la forme parité / contraction de tenseurs (non couvertes). Vérifications exactes pour n = 12 à 30. Dossier : `code/gen6/g6a/` |
| G6b | Identité à d directions (d = 3) pour les diagonales d'une boîte, avec un gain meilleur que les grilles | ❌ Aucune décomposition ne bat la grille. **Nouveau plancher (prouvé dans le modèle) : ρ_d ≥ 1/d** (au lieu de 2^(1−d) dans G4b), par sous-additivité du rang. **Relecture G6r : partiellement confirmé** — correction : le rang de l'aplatissement vaut 2^(2^(d−1)·λ), pas 2^(dλ) (c'est l'origine de l'ancien plancher 2^(1−d)) ; le plancher 1/d tient si une seule pondération q est commune à toutes les pièces et si log|I| = o(λ). Conséquence : ρ_d^(1/log₂ d) ≥ 1/2, donc les boîtes de dimension d ne battent jamais la grille au plancher. La borne spéculative 2^(m/2^O(√log m)) de G4b est **éliminée** : ce chemin aussi est limité à 2^Ω(m/log m). Recherche à λ = 1 : meilleur ρ₃ = 0,914, loin de la cible 0,577. Dossier : `code/gen6/g6b/` |
| G6c | Article de Claude : enlever le facteur log log n de la borne en n | ⚠️ Non résolu, mais la question est réduite. **Prop. R (prouvée) :** le coût dépend du log de la densité **moyenne** m/n, pas du degré maximal. La question devient : **Q1** — existe-t-il une sparsification qui préserve le comptage, en 2^(O(n/√log n)) morceaux de densité moyenne 2^O(k) ? Rééquilibrer ε ou brancher sur les variables de grand degré ne suffit pas (prouvé). **Prop. L :** les sparsifieurs « par parties, sans nouvelles variables » ont besoin d'une densité ≥ 1/(19ε) (sur des cliques disjointes). Échappatoire possible : des variables d'extension. **Q2** (ouverte) : des tests de parité avec O(n·2^O(k)) incidences au total. Dossier : `code/gen6/g6c/` |
| G6d | Méthode d'OpenAI avec une erreur différente en chaque point | ⚠️ **Ouvert, avec une piste positive.** **Lemme (prouvé) :** pour tout ensemble S de s directions, Σ rang(E_k) ≥ (N_S − d)·2^(ℓ−s) (argument du déterminant sur un bloc). Conséquence : le cadre ne descend jamais sous ≈ 2^(n (log log n)³/log n). **Les erreurs par point aident vraiment :** famille explicite de degré 1 avec ρ = 1/4 (contre → 1/2 avec une erreur commune), vérifiée exactement pour ℓ = 2 à 10. **Question ouverte précise :** existe-t-il des familles avec ρ = 2^(−Θ(d)) uniformément en ℓ ? Cela donnerait #SAT en 2^(O(n (log log n)³/log n)), **mieux que 2^(n/√log n)**. MILP : valeurs (ℓ−d)·2^(ℓ−d−1) pour les petits cas, mais le motif se casse à (5,2). Dossier : `code/gen6/g6d/` |

## 4 quinquies. Génération 7 (en cours)

| Voie | Objectif | État |
|---|---|---|
| G7a | Construire des familles à erreurs par point avec ρ = 2^(−Θ(d)) : **accélérer l'algorithme d'OpenAI** au-delà de 2^(n/√log n) | ❌ pour battre √log n, mais avec des résultats. **Nouvelle famille « OpenAI par morceaux »** (prouvée, vérifiée exactement pour ℓ ≤ 10) : ρ → 1/(2(d+1)) uniformément en ℓ ; avec un bloc central, ρ ≈ √(2/(πℓ))·e^(−d²/2ℓ), environ ℓ/d fois mieux que l'erreur commune d'OpenAI ; mais le temps reste 2^Θ(n (log log n)²/√log n) (seule la constante change). **Barrière A (prouvée) :** une seule fonction globale F ne gagne rien (ρ ≥ 1/4 − O(d/√ℓ), argument de Smolensky). **Barrière B (condition nécessaire) :** lien avec le **problème de Turán sur l'hypercube** : battre √log n demanderait des ensembles qui touchent tous les sous-cubes de dimension d+1 avec une densité ≈ ℓ^(−2), ce que personne ne connaît ; ρ = 2^(−Θ(d)) trancherait l'écart d'Alon-Krech-Szabó. **Correction de G6d :** sa conséquence 2(c) n'est justifiée que dans la classe commutative. Dossier : `code/gen7/g7a/` |
| G7b | Étendre le résultat G6a (Limits §9) au modèle de coût « programmation dynamique » | ⚠️ **Le modèle DP de l'article est trop généreux.** **Prop. O1 (prouvée) :** sans déformation, avec q_t = 0 et des groupes de taille √n entièrement mesurés, le coût DP « littéral » vaut ≤ 2^⌈√n⌉ ; avec le théorème 8.1 de l'article, l'optimum sur expanseur est donc 2^Θ(√n), de façon triviale. Cause : le modèle compte un facteur 2 par lecture de groupe, quel que soit son nombre de lignes, alors qu'une vraie table DP a au moins 2^(τ+1−s) entrées. **Modèle corrigé proposé** (largeur pondérée par le nombre de lignes W_t) : théorèmes C, D, F prouvés ; **corollaire G** : sous une hypothèse de régularité, coût ≥ 2^Ω(n^(2/3)), l'analogue DP du résultat G6a. **Ouvert :** une borne sans hypothèse de régularité (elle échoue pour 32 à 44 % des groupes dans les schémas de G6a). Vérifications exactes pour n = 12 à 30. **Pas encore relu.** Dossier : `code/gen7/g7b/` |
| G7c | Rédiger une note mathématique propre sur le résultat G6a, pour une relecture par des experts humains | ✅ **Note rédigée** : `resultats/note_limits_s9.md` (environ 6 pages, en anglais, Markdown + LaTeX). Elle contient aussi le lemme 4 (sep(F) ≥ c·n pour 3-SAT aléatoire de densité Δ > 3) et un séparateur pondéré pour les degrés non bornés. **Petites corrections trouvées :** le second terme du théorème B est 2n/3 − s_t (pas (2n/3 − s_t)/2) ; la constante exacte est c₃ = 2(1−e^(−2))/(36√2) ≈ 0,03397 ; c₀ = (1−e^(−2))/4 est prouvée, pas seulement vérifiée numériquement ; « équivalemment, √s lignes par groupe » vaut seulement pour des groupes équilibrés. |
| G7d | Question Q1 : sparsification qui garde le comptage avec une densité moyenne bornée (enlever log log n) | ❌ Non résolu. **Q1 sans restriction est équivalente à l'objectif lui-même** (prouvé) : ce n'est donc pas une étape plus facile. **Lemme S (prouvé) :** ajouter des variables ou des clauses n'aide jamais ; découper les variables en copies redonne la borne en m, en pire. **Prop. E (comptage) :** pour presque tout graphe de degré moyen d, toute CNF (même avec variables auxiliaires) dont la projection donne ses ensembles indépendants a au moins nd/(8k') clauses ; mais elle ne couvre qu'un seul morceau. Obstruction explicite hors de portée (elle donnerait une borne inférieure de circuit ω(n)). **Nouvelle question concrète Q1′ :** compter les ensembles indépendants d'un graphe de degré moyen d en temps 2^(O(n/d)). Dossier : `code/gen7/g7d/` |

## 4 sexies. Vérification Lean du lemme 1 (G8L) ✅

Fichier : `lean/PnpLean/WitnessDet.lean` (234 lignes, Lean 4.37.0-rc1 + Mathlib). Compile sans erreur ; aucun `sorry`, `admit`, `axiom` ni `native_decide` ; `#print axioms` = [propext, Classical.choice, Quot.sound] (les axiomes standards). Vérifié par le coordinateur après le sous-agent.

| Théorème Lean | Contenu |
|---|---|
| `det_ne_zero_of_eval_zero` | Φ(0) = I ⇒ det Φ ≠ 0 |
| `natDegree_det_le` | entrées de degré ≤ R ⇒ deg det Φ ≤ nR |
| `pow_dvd_det_of_nullity` | (X − t)^(dim ker Φ(t)) divise det Φ |
| `sum_nullity_le` | **Σ_{t∈T} dim ker Φ(t) ≤ nR** pour tout ensemble fini T |
| `deficiency_le_nullity_add` | borne de rang par blocs : Σ_γ (a_γ − \|I_γ\|)₊ ≤ dim ker A + Σ_γ dim R_γ (déficience − ρ ≤ dim ker) |

**Extension G9L ✅ (partie déterministe, vérifiée par le coordinateur) :** fichiers `Scheme.lean`, `TheoremB.lean`, `TheoremA.lean` (1 001 lignes au total avec WitnessDet). `lake build PnpLean.TheoremA` réussit ; aucun `sorry`, `admit`, `axiom` ni `native_decide` ; tous les théorèmes dépendent seulement de [propext, Classical.choice, Quot.sound].
- Structure `Scheme` fidèle à la définition 5.1 (déformations P, points, ordres q, fonctionnelle Λ de reconstruction d'Hermite, groupes, α, sous-espaces R, condition (3)).
- `lemma1` (complet, avec la traduction schéma → matrices), `factR`, `avg_lemma1` (moyenne sur toutes les cartes témoins), `exists_point`, **`lemma3` (prouvé)**, `theoremB` et `theoremB_bounded` (avec le **lemme 2 comme hypothèse explicite**, énoncé exactement comme dans la note), **`theoremA` (sans condition)**.
- Écarts mineurs : groupes indexés par un type fini fixe ; sep(F) et W(F) donnés comme bornes inférieures valables pour tout séparateur (équivalent) ; l'égalité « moyenne uniforme sur les cartes témoins = choix indépendants » est argumentée, pas prouvée en Lean.
- **Reste hors de Lean :** le lemme 2 (anti-concentration), le lemme 4 (3-SAT aléatoire), le corollaire C (asymptotique).

**Non formalisé dans la première version :** le passage du schéma de la note (polynômes P_o, carte témoin ω) à ces matrices (fait ensuite par G9L, voir ci-dessus). C'est une traduction directe, mais elle reste vérifiée seulement à la main. Les lemmes 2 à 4 (probabilités, séparateurs, 3-SAT aléatoire) ne sont pas formalisés.

**Réorganisation du 24 novembre :** les fichiers Lean sont classés par article. `lean/PnpLean/Article1/` (WitnessDet, Scheme, TheoremB, TheoremA : vérifiés) ; `Article1/EnCours/` (UAvg, Anticonc : lemme 2, en cours) ; `Article2/EnCours/PerPoint.lean` (lemme du déterminant par bloc pour les erreurs par point). Guide : `lean/README.md`. Cible : `lake build PnpLean.Article1.TheoremA`. Les deux versions de l'article citent les nouveaux chemins.

## 4 septies. Relecture indépendante de la note (G8R) et génération 8

**Rapport G8R** (`resultats/rapport_relecture_G8R.md`) : **révision majeure**. Les mathématiques tiennent (toutes les preuves refaites ligne par ligne ; c₀ = 0,21617, c₃ = 0,033967). Corrections obligatoires principales :
1. Le modèle DP de Dias-Lima est **dégénéré** (pas « ouvert ») : un schéma trivial coûte 2^√n (confirme la prop. O1 de G7b).
2. La borne 2^Ω(n^(2/3)) porte sur F lui-même, pas sur F\* ni sur les feuilles de la sparsification ; comparer avec le théorème 8.1, le lemme 8.2 et le théorème 7.1(ii) de Dias-Lima.
3. Énoncés trop forts : le titre (« cannot absorb ») et « ils résolvent le cas d'une ligne ».
4. Fidélité à la source : arêtes du graphe d'interaction, réserves du modèle de coût, définition de p_{e,v,v'}.
5. Vérifications numériques : 18 tests réels et non 36 ; l'étape de moyenne du théorème B n'est jamais vraiment testée ; pas de commandes reproductibles.
6. Il faut un auteur humain responsable.
7. Sections manquantes : introduction avec contributions, travaux liés, conclusion, bibliographie.
8. Retirer le matériel interne (« Issues found », références au journal).
9. Conflits de notation (X, Pr, Φ, w).

Simplification proposée : l'hypothèse M ≥ n peut être supprimée.

**Auteur humain de la note :** Lucas Requena, étudiant en informatique, contact dans la note : l.requena@outlook.fr (décision du 24 novembre 2026), avec une mention claire du travail fait par les agents IA.

**Génération 8 (en cours) :** G8a (terminer Limits §9), G8b (compter les ensembles indépendants), G8c (réfutations Extended Frege / TC⁰-Frege), G8d (Turán sur l'hypercube). Liste complète des pistes : `PISTES.md`.

- **G8c ✅ (conditionnel) :** aucune borne 2^o(n) prouvée en EF, TC⁰-Frege ou Frege. Observation nouvelle : pour δ fixé, la récursion de Claude a une profondeur constante ; Count(F) est une formule arithmétique de profondeur constante c(δ,Δ) et de taille 2^(δn)·poly(n). Énoncé conditionnel : sous l'hypothèse (S) (degré ≤ 1 en chaque table de poids, à vérifier) et des preuves EF de taille s de Pin_F(b), EF réfute F en taille 2s + poly(n)·2^(δn). Point ouvert : les formes closes des sous-appels aux niveaux profonds (vérifié au niveau 1 : 2 000 cas + les 4 096 cas de l'identité des huit carrés). Dossier : `code/gen8/g8c/`.
- **G8b ❌ (avec résultats) :** Q1′ (compter les ensembles indépendants en 2^O(n/d)) est nécessaire mais pas suffisante ; la bonne cible est **Q1\*** : compter les ensembles indépendants de tout graphe à n·polylog(n) arêtes en 2^O(n/√log n). Prouvé : le réseau direct est pire (U = Θ(nd)) ; **la voie des tests de parité (Q2) est fausse** sur la famille de référence (lemme CA). Aucune borne inférieure ne bloque Q1\* depuis la réfutation de l'ETH. Expériences (n ≤ 110) : travail ≈ 0,22·n·log₂²d/d. Dossier : `code/gen8/g8b/`.
- **G8a ⚠️ (niveau P, non relu) :** ni (a) ni (b) n'est fermé. **Théorème B2** (témoin adverse) : K = max_t(ρ_t+q_t+1) > min(n/4, c_b·sep(F)/√|X|) ; donc un schéma qui atteint n^(2/3) a besoin de |X| = Ω(n^(2/3)) points et de plans qui changent d'un point à l'autre. **Proposition S (barrière) :** la méthode du déterminant témoin ne peut pas dépasser n^(2/3) ; la perte √s est intrinsèque à la méthode. Aucun schéma réel de coût n^(2/3) trouvé (les familles testées retombent à Θ(n)). Modèle DP corrigé : **lemme H** (les régions en surplus portent leur surplus sur leur frontière) et corollaire H1 ; l'hypothèse de régularité n'est pas enlevée ; proposition P (barrière par point) seulement esquissée. Piste suivante : la borne de Rado. Tests : 5 160 vérifications du lemme H (n = 24 et 30), qui détectent des erreurs mais ne prouvent rien. Dossier : `code/gen8/g8a/`.
- **G8d ⚠️ (niveau P, non relu ; classe commutative seulement) :** l'optimalité de √log n n'est **pas** prouvée. **Lemme nouveau (Q_{d+2}) :** tout sous-cube de dimension d+2 contient au moins 4 mauvaises arêtes (l'argument combinatoire seul en donne 3) ; preuve : le même polynôme Ψ_Q calculé par deux découpages. Corollaire : ρ ≥ 2^(−d)/(d+2), un facteur ≈ 2 de mieux que G6d, sans conséquence nouvelle sur le temps. Critère exact de faisabilité (alternative de Fredholm) ; cas d = 1 : relèvement à un revêtement Z^(ℓ−1). **Correction de la cible :** une borne ρ ≥ c/poly(d) pour tout ℓ est fausse (la construction d'OpenAI la contredit). Ce qu'il faudrait : (I′) ρ(Cδ², δ) ≥ δ^(−a) et (II) super-multiplicativité ; ensemble, elles donneraient √log n optimal à (log n)^o(1) près. Littérature : Ellis-Ivan-Leader (2024) limite tout argument purement combinatoire à 2^(−Θ(d)) (non vérifié dans le texte). Dossier : `code/gen8/g8d/`.

**Génération 9 :**
- **G9e ⚠️ (niveau P, non relu ; tout est CONDITIONNEL aux deux prépublications non validées) :** conséquences nouvelles de la chute de l'ETH. Règle : une réduction vers M clauses garde un gain seulement si M = n·o(√log n). **Meilleur candidat :** problèmes de connexité sur les graphes creux par comptage — nombre d'arbres de Steiner à ℓ arêtes, de cycles et de chemins de longueur ℓ modulo p en 2^(O(m log p/√log m)) (déterministe) ; arbre de Steiner, plus long cycle, plus long chemin (décision) en 2^(O(m log log m/√log m)) (aléatoire, erreur d'un seul côté). Le cycle hamiltonien est déjà revendiqué par d'autres (dépôt xennygrimmato/eth-hamiltonicity) ; Steiner et plus long chemin semblent nouveaux (nouveauté **non établie** : articles voisins non lus). Autres : taux explicites pour #CNF-SAT de largeur non bornée, #CSP, Vertex Cover paramétré, Max-Cut. Gain perdu pour : k-Clique, nombre chromatique, subset sum, etc. Identité de comptage testée sur 6 petits graphes seulement. Dossier : `code/gen9/g9e/`.
- **G9a ❌ (piste A ; niveau P, non relu) :** pas d'amélioration sous ρ = 1/√2. Dans le cadre de l'annexe B, α = 1/2 (les partitions) est **optimal** parmi les croissances en loi de puissance. **Obstruction A (comptage de rang, prouvée) :** le croisement a un rang N², ce qui force L(4h) ≥ 2·L(h), donc α ≥ 1/2 ; en général ρ ≥ max(2^(α−1), 2^(−α)) ≥ 1/√2. La phrase de l'article (« ρ = 2^(α−1) ») est juste comme calcul, mais un α < 1/2 est impossible : carrés, puissances k-ièmes, puissances de 2, parts bornées sont tous exclus. **Obstruction B :** la localité (lemme 9) se casse si on restreint les modes (prouvé pour un oscillateur). Tests : identité fausse dans les 23 cas restreints. Reste non exclu (spéculatif) : des profils de croissance non concaves « en escalier ». Dossier : `code/gen9/g9a/`.
- **G9d ❌ (piste D ; niveau E, esquisses dans un modèle de coût) :** aucun schéma ne passe le mur. Classe 𝒢 = sommes signées de réseaux de tenseurs, chacun contracté selon son arbre. **Planchers connus :** termes de rang 1 : 2^(S/2) ; un réseau, un arbre : 2^(4S/9) ; **toute la classe 𝒢 : seulement 2^Θ(√m)**, donc 2^(S^(1−ε)) n'est exclu par aucune borne connue. **Triangles et multiplication de matrices :** gain d'un facteur constant seulement (e₀ = 2ω/9), mur inchangé. La conjecture du rang asymptotique ne donne pas 2^o(n) pour 3-SAT. **« Loi harmonique » (nouvelle, esquisse) :** si chaque terme global est un produit de pièces choisies indépendamment, alors u(K) ≤ Σ u(g_j) avec u = 1/e, donc temps ≥ 2^Ω(m/log m), pour toute arité de gadget et tout nombre de directions. **Ce qui reste :** des termes **non produits (« intriqués »)**. Il faudrait à la fois E1 (une décomposition jointe d'un tenseur de grille 4×h) et E2 (évaluation récursive des parties fusionnées). Lié au problème de Naslund (N ≤ P(N) ≤ N^√2). Dossier : `code/gen9/g9d/`.
- **G9b ❌ (piste B ; niveau P, non relu) :** décider n'apporte aucun gain dans ces cadres. La récursion de Claude marche sur tout anneau commutatif, mais les pièces de l'annexe B demandent une caractéristique p > h (division par h!), donc elles ne sont pas disponibles en caractéristique fixe ; l'identité des huit carrés s'annule sur GF(2). Plancher ρ ≥ 1/2 sur tout corps (prouvé). **Lemme de la trace (prouvé) :** explique précisément le phénomène de la caractéristique 2. **Sur GF(2), le meilleur ρ trouvé est 0,811** (≈ H(1/4), puissances tensorielles de l'identité à un bit, pièces nilpotentes ; formule de Wilson) : mieux que 0,896 sur Q dans cette petite classe, mais **pire que 1/√2 ≈ 0,707**. Valiant-Vazirani : inutile (un premier aléatoire p ~ n² ou des poids aléatoires suffisent) et coûteux. Hachage des égalités : non fiable. Pièces positives ou monotones : ρ = 1 (prouvé), la méthode a besoin d'annulations. **Question ouverte :** des pièces avec ρ < 1/√2 sur un corps fini fixe. Dossier : `code/gen9/g9b/`.
- **G9c ⚠️ (piste C ; niveau P, non relu ; relatif à l'annexe D de l'article de Claude, non revérifiée) :** le facteur log log n n'est pas enlevé en général. **Résultats positifs nouveaux :** théorème A (inclusion-exclusion tronquée exacte sans tests de parité si n̄ ≥ max e(S)) et théorème B (i_t(G) en 2^(O(e·h_t·log(4m/h_t)))). Corollaires : #IS en 2^(O(en)), **sans log d**, pour tout graphe avec α ≤ n/(log d·log log d), pour les graphes (n,d,λ) avec λ ≤ d/log²d, et pour les **graphes d-réguliers aléatoires** : la famille de référence de G7d/G8b n'est donc **pas** difficile. **Obstructions :** hybride (prop. H1 : trop de feuilles), zêta/Möbius (le rang ne dépend pas de la base), identité de l'étoile (Θ(nd) bits). **Noyau ouvert, plus étroit :** ensembles indépendants de taille entre n/(log d·log log d) et n/2 ; nouvelle famille de référence : graphes **bipartis** d-réguliers aléatoires. Dossier : `code/gen9/g9c/`.

**Décision de l'auteur (24 novembre) :** l'article 2 (erreurs par point) ne sera écrit que s'il y a un **résultat final** ; pas de note de résultats partiels.

**Relecture humaine (24 novembre) :** l'auteur a envoyé la version 2 de l'article 1 à F. Dias pour vérification. **La version 3 ne sera faite qu'après son retour.** Les corrections et ajouts prévus sont listés dans `PISTES.md` (section « En attente »).

**Génération 10 :**
- **G10f ✅ Relecture des points non relus de l'article + borne pour F\* :**
  - (a) Note sur le théorème 7.1(ii) de Dias-Lima : **correcte**, mais deux corrections d'équité à faire : le facteur b est superflu (q_t + μ_t·s_t ≥ sep(F)) ; la borne de Dias-Lima est **meilleure que la nôtre pour les petits groupes** (s_t ≲ 867 pour b = 3), et l'article doit le dire.
  - (b) Annexe A : **chaque étape et chaque constante vérifiée** ; deux précisions mineures à ajouter (C_Δ ≥ 3 ; η indépendant de Δ).
  - (c) Bibliographie : **aucune entrée fausse** ; 15 entrées confirmées par recherche, 7 crues correctes de mémoire ; 4 attributions à adoucir (méthode du 4ᵉ moment = simple Hölder ; « introduits par Lipton et Tarjan » → « remontent au moins à » ; Bod93 ; GLR82).
  - **Nouveau (niveau P, non relu) : sep(F\*) ≥ c·n\* est prouvé** pour Δ ≥ 2, avec η\* = 2,418·10⁻⁴, par un lemme de transfert et un nouveau lemme de premier moment (la réduction naïve par contraction échoue). **Corollaire :** la borne 2^Ω(n\*^(2/3)) vaut aussi pour le schéma appliqué à F\*, ce qui rend la comparaison avec le théorème 8.3(a) de Dias-Lima **directe**. Réserve : F\* est pris selon la description de Dias-Lima, non vérifiée dans la prépublication d'OpenAI. Dossier : `code/gen10/g10f/`.
- **G10a ✅ Relecture adverse de G8a et G7b** (preuves refaites, code indépendant) :
  - G8a : **théorème B2 confirmé** (c₃ = 0,0772 ; demande b ≥ 2). **Proposition S : en partie** — le résultat principal (s = n^(2/3)) est correct ; 2 erreurs de constantes ((iii) fausse pour s proche de n ; « r ≤ 6s » demande s ≥ 3b ln n) ; c'est une configuration et non un schéma ; « barrière pour la méthode » est une description juste ; « ne dépasse pas n^(2/3) » vaut pour l'exposant seulement (facteur √log). **Lemme H confirmé** (6 360 vérifications). Corollaire H1 : en partie (vaut pour une décomposition W-optimale, pas pour « toute » décomposition). Proposition P : en partie (esquisse ; l'affirmation n^(3/4) non vérifiée).
  - G7b : **théorème C : borne principale confirmée, mais la clause « de plus » est RÉFUTÉE** (erreur de 1 : il faut τ + 1 − s − |Z| et non τ + 2 ; contre-exemple avec n = 6). **Théorèmes D, F et lemme E confirmés.** Corollaire G : en partie (hérite de l'erreur de 1, sans effet asymptotique ; doit préciser la décomposition W-optimale).
  - **Le modèle DP corrigé est seulement en partie valide :** 2^(W_t) minore seulement l'élimination « par groupes entiers ». Un éliminateur ligne par ligne fait mieux (exemples exacts : 12 contre 6, 15 contre 7, 16 contre 10). Le 2^(w_t) de l'article de Dias-Lima résiste à ce découpage, pas W_t.
  - Dossier : `code/gen10/g10a/`.
- **G10b ✅ Relecture et Lean des lemmes « erreurs par point » (base de l'article 2) :**
  - **Lemme du déterminant de G6d : confirmé**, dans le cas général (non commutatif) ; τ ≠ 0 n'est pas nécessaire ; plusieurs points par direction sont permis. Conséquence (a) confirmée (ℓ ≥ d+1) ; (b) confirmée à l'arrondi près (ρ ≤ 1/4) ; **(c) non prouvée en général** (« les rangs s'additionnent sur des blocs disjoints » est faux ; vrai seulement dans la classe commutative).
  - **Lemme Q_{d+2} de G8d : confirmé** pour tout d et tous τ distincts non nuls (classe commutative, un point par direction), avec un argument plus court. Corollaire ρ_moy ≥ 2^(−d)/(d+2) confirmé.
  - **Lean (vérifié par le coordinateur, sans sorry) :** `lean/PnpLean/Article2/PerPoint.lean` (`g6d_determinant_lemma`) et `QdPlus2.lean` (`qd2_lemma`), 522 lignes. Non formalisé : la réduction de la classe commutative, les corollaires, le lien avec l'algorithme d'OpenAI.
  - Dossier : `code/gen10/g10b/`.
- **G10e ✅ Lean : le théorème B est maintenant prouvé SANS hypothèse sur le lemme 2** (vérifié par le coordinateur ; fichiers promus dans `lean/PnpLean/Article1/`) : `UAvg.lean` (moyennes uniformes sur des produits finis, indépendance), `Anticonc.lean` (`bernoulli_anticonc` : c₀·min(V, √V) ≤ E|X − a| avec c₀ = (1 − e⁻²)/4 ; symétrisation, 2ᵉ et 4ᵉ moments, Hölder ; forme finie, sans théorie de la mesure), `Lemma2.lean` (`avg_home_ind`, `avg_Icard`, `lemma2 : Lemma2Hyp c0`, **`theoremB_unconditional`**, `theoremB_bounded_unconditional`). Total Lean des deux articles : 2 497 lignes, 0 `sorry`. Reste hors de Lean pour l'article 1 : le lemme sur les 3-CNF aléatoires, le corollaire asymptotique, la proposition sur le modèle PD.
- **G10c ✅ Relecture et nouveauté de G9e** (code indépendant : ≈ 570 000 vérifications de l'identité de comptage, 1 400 formules comptées) : la conséquence principale est **mathématiquement confirmée** (identité, taille de l'encodage sans perte de log m, analyse d'erreur, exposants). **Nouveauté : l'arbre de Steiner est DÉJÀ REVENDIQUÉ** (arXiv 2611.02548, Wu-Wang, 2 novembre 2026 : comptage exact et déterministe, plus fort, à partir du théorème d'OpenAI). **Plus long chemin / cycle de longueur ℓ (non induits) : aucune revendication trouvée** ; le dépôt eth-hamiltonicity les liste comme « conjecturés ». Gain honnête : premier 2^o(n) pour le plus long chemin/cycle sur les graphes à m = O(n), conditionnel, aléatoire, avec une perte log log. Dossier : `code/gen10/g10c/`.
- **G10d ✅ Relecture adverse de G9c** (simulation indépendante du réseau, 90 cas ; 42 exécutions du pipeline du théorème B) : **aucune affirmation réfutée.** Théorèmes A et B et corollaires B1 à B4 **confirmés**, relativement au lemme 8 et au théorème D(a) de l'article de Claude (non relus). Proposition H1 confirmée (269 cas). Deux résumés trop forts : prop. S (« Θ(nd) pour tout regroupement » est faux tel quel : il faut des groupes équilibrés ; contre-exemple donné) et l'idée 2 (borne heuristique). **Contexte honnête :** l'énumération simple donne déjà 2^(O(n log²d/d)) sur les graphes réguliers aléatoires ; le gain à d ≈ √log n est donc un facteur polylog(d). Dossier : `code/gen10/g10d/`.

**Génération 11 :**
- **G11a ✅ Relecture adverse du théorème sep(F\*) (G10f) : aucune erreur mathématique.** Obstruction confirmée (exemple n = 17, m = 18). Lemme de transfert confirmé (405 752 vérifications exhaustives, borne atteinte). Lemme de premier moment confirmé (κ = 0,392437). Théorème et corollaire confirmés (η\* = 2,41827·10⁻⁴). **Point faible : la définition de F\*.** Dias-Lima disent seulement « copies de variables sur des cycles de Tovey, au plus trois occurrences » ; « une copie par occurrence » est notre inférence : à énoncer comme hypothèse. Corrections : une note de G10f est fausse pour Δ ≤ 3 (variante à copie partielle) ; la remarque sur les feuilles couvre seulement la première feuille et demande une ligne de plus ; conflit de notation n₀ ; `check.py` teste une inégalité plus faible et ne doit pas être cité. Dossier : `code/gen11/g11a/`.
- **G11b 🔎 THÉORÈME ANNONCÉ pour l'article 2 (niveau P ; relu par le coordinateur, sans erreur trouvée ; relecture adverse G13a et Lean G13b en cours) :** dans le cadre général (non commutatif) des erreurs par point, avec un point par direction : **Σ_k rang(E_k) ≥ N_m := Σ_{j≤m} C(ℓ, ≤ j)**, avec m = ⌊(ℓ−d−1)/2⌋. Preuve en 5 étapes (argument de type Smolensky sur les certificats duaux) : un espace A_m de vecteurs lignes polynomiaux a_x(t) = ∏_{k : x_k=1}(1 − t/τ_k)·r_x(t) de dimension N_m ; si la somme des rangs est plus petite, un a ≠ 0 annule toutes les erreurs ; alors Φ(t) = a(t)·J(t)·c(t)ᵀ (c_y = a_ȳ) a un degré ≤ ℓ−1 et ℓ racines, mais Φ(0) = Σ a_x(0)² > 0. **Conséquences annoncées :** ρ_max ≤ ℓ^(−C) force d ≥ (1−o(1))·√((2C−1)·ℓ·ln ℓ), atteint par la construction en couches de G7a ; **la perte √log n est donc optimale** dans ce cadre (un point par direction ou multiplicités égales) ; les erreurs par point ne gagnent qu'un facteur au plus ℓ par rapport à une erreur commune. **Réserves :** l'étape de positivité demande un corps ordonné (ℚ, ℝ) : elle échoue sur les corps finis, et sur ℂ avec des τ complexes (à préciser) ; l'énoncé sur le temps de calcul vaut dans le modèle de coût de G6d seulement ; nouveauté non vérifiée. Dossier : `code/gen11/g11b/`.
- **G13b ✅ Lean : le théorème de G11b est prouvé** (vérifié par le coordinateur ; `lean/PnpLean/Article2/SqrtLog.lean`, 348 lignes, sans `sorry`, axiomes standards) : `SqrtLog.sum_rank_ge_N`, exactement l'énoncé fixé à l'avance, sans hypothèse ajoutée ; la division par t est remplacée par la plus grande puissance X^v qui divise tous les a_x ; `N_eq` donne la forme binomiale de N. Le théorème est donc **certain sur ℚ**. Restent hors de Lean : les conséquences asymptotiques et l'énoncé sur le temps de calcul.
- **G11c ⚠️ Cas général (non commutatif) (niveau P, non relu ; bornes supérieures numériques seulement) :** formule exacte : le minimum de rang(J(τ_k) − K) sur les K faiblement agréables vaut r_k = max(rang B01, rang B10). **Lemme Q_{d+2} étendu au cas général** (théorème A, caractéristique ≠ 2) : pour ℓ = d+2, Σ_k r_k ≥ 4, atteint ; lemme de la trace ; théorème de rigidité B′. **La classe commutative n'est PAS optimale :** à (ℓ, d) = (4, 1), le cas général atteint un total de 11 (numériquement), contre 12 dans la classe commutative ; minorant prouvé : 10. Le corollaire normalisé du lemme commutatif est faux par ensemble de d+2 directions ; ρ_moy ≥ 2^(−d)/(d+2) n'est pas prouvé en général. **Cohérent avec le théorème de G11b** (N_m = 6 ≤ 10). Dossier : `code/gen11/g11c/`.
- **G9W ✅ Article, version 2 :** `resultats/article/main.tex` et le PDF (26 pages). Titre : « Measurement Rows Absorb at Most a √s Factor of the Cut: A Determinant Bound for Clause-Product Interpolation Schemes ». Auteur : Lucas Requena. Les 10 corrections obligatoires du rapport G8R sont traitées. Changements mathématiques : hypothèse M ≥ n supprimée (deg det Φ ≤ n_occ·R ≤ MR ; Lean garde n ≤ M, écart signalé) ; corollaire partie 3 : 2^(min(n/6, …)) ; remarque « équilibrée » avec γ ≤ 1/6. Tests refaits : les 4 scripts passent ; Σ(E def − ρ)₊ = 0 dans les 12 exécutions (l'étape de moyenne n'est pas réellement testée) ; 18 vérifications du relecteur. Vérifié par le coordinateur : compilation sans erreur, en-tête, résumé, section IA, bibliographie (2 identifiants ajoutés). **Non relu :** une note de bas de page sur le théorème 7.1(ii), l'annexe A, les références classiques (citées de mémoire).

## 4 bis. Pistes restantes pour chercher plus loin

| Piste | Pourquoi elle reste ouverte | Première expérience simple |
|---|---|---|
| Certificats de Feige-Kim-Ofek | Des réfutations courtes existent à m ≈ n^1,4, sous le seuil spectral n^1,5 | Construire ces certificats pour n ≤ 200 et voir quelles parties une recherche peut trouver |
| Algèbre sur l'anneau Z (entiers) | L'argument « réduction modulo p » ne s'applique pas directement | Degré de Nullstellensatz sur Z pour n = 8 à 16 |
| Extensions **structurées** (pas gloutonnes) | 3a a testé seulement des définitions simples | Chercher des définitions qui raccourcissent une preuve connue sur de petites formules |
| Res(⊕) avec des parités **longues** | 3b a testé seulement des parités de 2 ou 3 variables, et n ≤ 60 | Brancher sur des parités de log n variables, avec des masques plus grands que 64 bits |

**Leçon commune aux étapes 3a à 3d :** chaque outil (extension, parité, coupe, anneau) aide seulement quand la formule a la structure que cet outil sait utiliser (XOR pour Gauss, cliques pour les coupes). 3-SAT aléatoire n'a aucune de ces structures. Un algorithme qui prouve P = NP doit donc trouver une structure **cachée** dans les formules aléatoires. Selon l'hypothèse de Feige, cette structure n'existe pas quand m = O(n).

---

## 5. Contexte externe (novembre 2026)

### 5.1 Actualité directement liée à P vs NP (source : recherche web ; les articles sont bloqués par le réseau de la session)
- **P vs NP reste ouvert.** Le Clay Institute le liste toujours comme non résolu.
- **L'hypothèse du temps exponentiel (ETH) est réfutée (annonce, octobre 2026) :**
  - **OpenAI**, « The Exponential Time Hypothesis Is False » (17 octobre 2026, révisé le 25) : compte les solutions d'une k-CNF en temps 2^(O(n log log n / √log n)), et 2^(O(n/√log n)) à densité bornée. Méthode : contraction de tenseurs et interpolation.
  - **Anthropic (Claude)**, « A 2^{o(n)} Time Algorithm for 3-SAT » (26 octobre 2026, 49 pages, lu en entier pour l'aperçu) : algorithme déterministe, compte aussi les solutions. Trouvé le 16 octobre par un modèle Claude de recherche interne (14 M de tokens, sans aide humaine). Temps : 2^(O(min(k² n log log n/√log n, k m/√log m))).
    - **Méthode réelle (correction : les extraits web étaient faux) :** #3-SAT devient une somme de produits sur le graphe des clauses (un réseau de tenseurs). La récursion place K nœuds sur une grille g × g. Pour chaque rectangle, elle réécrit le produit des deux tests d'égalité diagonaux avec l'**identité des huit carrés** (Degen 1818, Graves 1843, octonions) : 2[x=v][y=u] = [x⊕y⊕u⊕v=0](σ(x,y)σ(v,u) + σ(x,u)σ(v,y)). Chaque terme est fusionné par lignes ou par colonnes, puis équilibré, puis la récursion continue.
    - **Goulot :** à chaque niveau, l'exposant est multiplié par ρ = 11/12 (huit carrés) ou ρ = 1/√2 (pièces construites avec l'algèbre de Heisenberg, annexe B). D'après l'article, ce type de construction donne au mieux ρ → 1/2, donc environ 2^(m/log m) : **loin du temps polynomial**.
    - L'inclusion-exclusion sur des clauses violées sert seulement à la borne en n (annexe D), après la sparsification.
    - Gain réel seulement pour log m > 8000 (m = n) : aucun effet pratique.
  - Vérification : formalisations en Lean annoncées. Évaluation par les pairs en cours.
  - **Limites :** le résultat est seulement asymptotique (une analyse estime un gain seulement au-delà de n ≈ 10^126). Il ne donne pas P = NP, car 2^(n/√log n) reste plus que polynomial.
- **Emanuele Viola** (blog CACM, 17 novembre 2026) : les résultats de l'IA ne prouvent pas P = NP, mais ils renforcent son opinion que P = NP.
- **Feige (version quantique) :** Chan et Xu (octobre 2026) réfutent fortement 3-SAT **quantique** aléatoire à densité constante, en temps polynomial classique. La version classique de l'hypothèse de Feige reste ouverte.

### 5.2 Autres résultats récents en complexité
- **Williams (STOC 2025) :** un calcul en temps t peut être simulé en espace O(√(t log t)). C'est un petit progrès sur P vs PSPACE.
- **Res(⊕) :** bornes exponentielles seulement pour des versions limitées (régulière : Efremenko-Garlik-Itsykson, STOC 2024 ; fragment : Alekseev-Itsykson, STOC 2025). Le cas général reste ouvert.
- **Circuits algébriques :** bornes super-polynomiales pour la profondeur constante (Limaye-Srinivasan-Tavenas, JACM 2025).

### 5.3 Histoire (rappel)
Cook (1971) et Levin (1973) : NP-complétude. Karp (1972) : 21 problèmes NP-complets. Baker-Gill-Solovay (1975) : relativisation. Razborov-Rudich (1997) : preuves naturelles. Aaronson-Wigderson (2008) : algébrisation. Mulmuley-Sohoni (2001) : GCT. Williams (2011) : NEXP ⊄ ACC⁰.

### 5.4 Autres problèmes du millénaire
- Hypothèse de Riemann généralisée : preuve **annoncée** par Levent Alpöge avec l'aide de Claude (novembre 2026), en Lean. Pas encore validée par la communauté.
- Navier-Stokes : preuves **annoncées** par OpenAI (septembre 2026, critiquée) et par Buckmaster et Alpöge (20 novembre 2026).
- Le Clay Institute indique : « solution claimed; prize pending ».
- Sources : blog de Terence Tao (20/11/2026), King's College London, arXiv 2610.08144.

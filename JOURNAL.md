# Journal de recherche : P vs NP

Date de début : 23 novembre 2026.
Objectif choisi par l'utilisateur : chercher **P = NP**, c'est-à-dire un algorithme polynomial pour 3-SAT, avec une preuve.

Ce journal contient toutes les voies explorées, les résultats et la cause de chaque échec.
Les programmes et les données sont dans le dossier `code/`.

> **Rappel honnête.** Aucune voie ci-dessous ne résout P vs NP. Chaque résultat est soit un échec mesuré, soit un échec expliqué par un théorème connu.

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
| G7a | Construire des familles à erreurs par point avec ρ = 2^(−Θ(d)) : **accélérer l'algorithme d'OpenAI** au-delà de 2^(n/√log n) | en cours |
| G7b | Étendre le résultat G6a (Limits §9) au modèle de coût « programmation dynamique » | en cours |
| G7c | Rédiger une note mathématique propre sur le résultat G6a, pour une relecture par des experts humains | ✅ **Note rédigée** : `resultats/note_limits_s9.md` (environ 6 pages, en anglais, Markdown + LaTeX). Elle contient aussi le lemme 4 (sep(F) ≥ c·n pour 3-SAT aléatoire de densité Δ > 3) et un séparateur pondéré pour les degrés non bornés. **Petites corrections trouvées :** le second terme du théorème B est 2n/3 − s_t (pas (2n/3 − s_t)/2) ; la constante exacte est c₃ = 2(1−e^(−2))/(36√2) ≈ 0,03397 ; c₀ = (1−e^(−2))/4 est prouvée, pas seulement vérifiée numériquement ; « équivalemment, √s lignes par groupe » vaut seulement pour des groupes équilibrés. |
| G7d | Question Q1 : sparsification qui garde le comptage avec une densité moyenne bornée (enlever log log n) | ❌ Non résolu. **Q1 sans restriction est équivalente à l'objectif lui-même** (prouvé) : ce n'est donc pas une étape plus facile. **Lemme S (prouvé) :** ajouter des variables ou des clauses n'aide jamais ; découper les variables en copies redonne la borne en m, en pire. **Prop. E (comptage) :** pour presque tout graphe de degré moyen d, toute CNF (même avec variables auxiliaires) dont la projection donne ses ensembles indépendants a au moins nd/(8k') clauses ; mais elle ne couvre qu'un seul morceau. Obstruction explicite hors de portée (elle donnerait une borne inférieure de circuit ω(n)). **Nouvelle question concrète Q1′ :** compter les ensembles indépendants d'un graphe de degré moyen d en temps 2^(O(n/d)). Dossier : `code/gen7/g7d/` |

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

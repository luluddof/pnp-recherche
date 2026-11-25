# Pistes de recherche

Ce fichier liste toutes les pistes : celles en cours, celles pas encore explorées et celles fermées.
Mise à jour : 24 novembre 2026.

**Rappel :** aucune de ces pistes ne mène directement à P = NP. Ce sont des étapes qui peuvent donner des résultats vérifiables.

Légende : 🔄 en cours · 🆕 pas encore explorée · ⏸️ en attente · ✅ terminée · ❌ fermée (mur prouvé ou mesuré)

---

## Plan des articles (décision de l'auteur, 24 novembre 2026)

**Le détail de la matière de chaque article, et les autres candidats, sont dans `ARTICLES.md`.**

**Règle : on attend la validation de l'article 1 (relecture de F. Dias) avant d'écrire en entier les articles suivants.** La recherche, les relectures et Lean peuvent continuer.

| Article | Sujet | État | Portée |
|---|---|---|---|
| **1** | Schémas « produit de clauses » (question du §9 de Dias-Lima) | Version 2 envoyée à F. Dias ; pré-version 3 en brouillon | Borne 2^Ω(n^(2/3)) dans le **modèle de coût produit** ; rien sur la programmation dynamique (sauf : modèle dégénéré), ni sur les autres voies |
| **2** | Erreurs par point dans la voie « interpolation » d'OpenAI | **En attente** (accord de l'auteur requis ; attendre G14a et G14d, puis la validation de l'article 1). Théorème prouvé en Lean sur tout corps | Σ_k rang(E_k) ≥ N(ℓ, m) pour un point par direction ; √log n optimal **dans le modèle de coût de l'exposé d'O'Donnell** ; **État des échappatoires :** corps finis : fermé (Lean) ; multiplicités inégales : fermé (**Lean**, G16L) ; effondrement de rang par L/R′ : fermé sur papier (G14a, non relu) ; degré réduit D < nd : fermé en partie (reste : degré certifié par annulation). Ne couvre pas : la contraction de tenseurs (futur article 3), l'algorithme de Claude **Ajout possible (G18b/G18R) :** borne supérieure en loglog n par un évaluateur modulo 3, contre borne inférieure en √loglog n pour une charge constante. |
| **3** (futur) | Limites de la voie « contraction de tenseurs » d'OpenAI (annexe D) | **Résultat central obtenu et vérifié par Lean (G18e + G18S + G19d, 24/11) : théorème A″ ; corollaire N ≥ 2^(0,564·M/√u) sur papier (P+).** Avant (G17a) : idée seulement : théorème A confirmé mais très faible ; théorème B en partie (axiome du modèle VT réfuté) ; un lemme propre (précision optimale en boîte noire) ; il manque une borne de rang 2^(M/√u) dans le modèle DF | Visé : le √log est optimal pour la famille des « effilages » dépendant de la valuation ; borne d'aplatissement dans le modèle DF. Ne couvrirait pas : d'autres règles de reconstruction, plusieurs niveaux, l'algorithme de Claude |

Avec les trois articles, les trois voies de la prépublication d'OpenAI seraient étudiées. Aucun ne dit rien sur P vs NP ni sur la validité de la prépublication.

## 0. L'échelle vers P = NP (temps pour 3-SAT)

**Précision d'échelle (G21c, 24/11) :** les deux prépublications donnent 2^(O(m/√log m)) en le nombre de clauses m. En le nombre de variables n, la borne générale est 2^(O(n·loglog n/√log n)) ; elle est 2^(O(n/√log n)) seulement à densité bornée. Dans ce fichier, « n/√log n » et « n/log n » se lisent à densité bornée (ou en m).

**Échelons intermédiaires (G21a, 24/11) :** avec un gain ρ par niveau, le temps de l'algorithme de Claude est 2^(O(m/(log m)^α′)), α′ = log₂(1/ρ). ρ = 1/√2 donne la puissance 0,5 (échelon 1) ; ρ = 2/3 donnerait 0,585 ; ρ = 1/2 donnerait 1 (échelon 2). Tout ρ < 1/√2 est donc un progrès sur la puissance du logarithme, pas seulement sur une constante.

| Marche | Temps | État |
|---|---|---|
| 0 | 2^(0,39 n) | Record avant octobre 2026 |
| 1 | 2^(n/√log n) | Atteint par OpenAI et Claude (annoncé, non validé) |
| 2 | 2^(n/log n) | Mur des méthodes actuelles (générations 4 à 6) |
| 3 | 2^(n^(1−ε)) | Demande une idée nouvelle (hors structure produit) |
| 4 | 2^(√n) | — |
| 5 | Quasi-polynomial puis polynomial (= P = NP) | La plupart des experts pensent que cette marche n'existe pas |

Priorité recommandée : pistes A et B, puis Q1\*.

## 1. En cours

| Id | Piste | Question précise | Origine |
|---|---|---|---|
| G27-B1 | ❌ **Brainstorming : constructions** (fini) | Rien dans le cercle ; aucun coin. Ansatz d'orbite Z2^k (découpage en N blocs) ; nouveau point numérique (62,12)@16, non exact. Dans les modèles par blocs, un coin demande Q ≥ N^((1+2δ)/3). `code/gen27/b1/notes.md` | Étape 2 |
| G27-B2 | ✅ **Routes hors du point équilibré** (fini, relu par G27-R1 : P+ avec corrections) | **ρ → 1/2 n'est pas nécessaire : étape 2 ⇔ pente finie au coin θ = 2** (niveaux déséquilibrés, formule α*). Nouveau problème : rigidité de rang de Kronecker du swap, δ(η) ≥ c·η ? `code/gen27/b2/notes.md` | Étape 2 |
| G27-T2 | ⚠️ **Frontière exacte à N = 10 à 12** (fini tôt) | Cinq nouveaux points exacts dont (25,5)@10 et (19,6)@8 ; famille A_a ⊗ A_b. Tendance indéterminée ; le coin n'est pas échantillonné aux petits N. `code/gen27/t2/notes.md` | G26-P11 |
| G27-R1 | ✅ **Relecture de G27-B2** (fini) | Résultat principal confirmé (niveau conception) : α* = 1 ⇔ γ = 1 ; pente finie c au coin ⇒ temps 2^(O(c·m·log log m/log m)). Aucune famille connue près du coin. `code/gen27/r1/notes.md` | G27-B2 |
| G27-T1 | ❌ **Annexe B avec états de grand poids** (fini tôt) | k_max(h, ∞) égal à une prédiction de comptage triviale, borné ; loi des spans de modes prouvée ; θ ne change pas. Non fini : h ≥ 8. `code/gen27/t1/notes.md` | G26-P1 |
| G27-T3 | ⚠️ **Problème du coin** (fini tôt) | Aucun coin linéaire ; toutes les familles (Bell, nombre de singulets, Pauli) ont un coin quadratique, pas meilleur que l'annexe B. Pas de borne δ ≤ C·η² (sauf théorème S pour les constructions de Bell symétriques, esquisse). Route ouverte, aussi difficile que la droite. `code/gen27/t3/notes.md` | G27-B2 |
| G25-B1 | ✅ **Brainstorming : constructions** (fini) | 12 idées, 7 tuées par calcul ; rien près du cercle. Restes : profil à saut convexe (2 %), familles au bord en ε (1,5 %), catalyseurs (1 %). La barrière de Cardy s'étendrait aux familles c → ∞ (argument non relu). `code/gen25/b1/notes.md` | consigne de l'auteur |
| G25-B2 | ✅ **Brainstorming : techniques de preuve** (fini) | 13 idées. Prouvé : N² ≤ p‖R‖_op + q‖Fl2 C‖_op. Top 3 : conjecture d'étalement (donnerait 2θ + T ≥ 4, donc ρ ≥ 2/3 comme mur) ; modèle à blocs ; angle de Friedrichs. `code/gen25/b2/notes.md` | consigne de l'auteur |
| G25-B3 | ✅ **Brainstorming : autres routes** (fini) | 10 idées ; top 3 : sparsification pour la décision seule (enlèverait le facteur loglog en n) ; quotient par grappes ; hypercliques récents comme gadgets à 3 directions. `code/gen25/b3/notes.md` | consigne de l'auteur |
| G25-B4 | ✅ **Brainstorming : rang inférieur** (fini) | « Meilleures pièces à chaque taille » : ne gagne rien (λ ≥ 10⁴ dans l'algorithme). Gains réels dans la comptabilité de la récursion : base de balayage + constantes exactes (seuil log₂ m* : 1,37·10⁵ → 2,4·10⁴) ; annexe D : 4·10¹⁰ → 1,6·10⁶. Reste astronomique. `code/gen25/b4/notes.md` | consigne de l'auteur |

Dossiers de travail de la génération 18 : `scratchpad/gen18/g18a` à `g18i`.

**État (24/11) :** la génération 18 est terminée : onze rapports (G18a à G18i, relecteurs G18R et G18S).

Note : les dossiers de travail des agents G17b à G17e s'appellent `g18a` à `g18d` (erreur de numéro au lancement ; la génération est bien la 17).

**Consigne de l'utilisateur (25/11, nuit) : nouveau brainstorming, traiter les sous-problèmes. Fenêtres d'utilisation : travailler jusqu'à la limite, reprendre à 2h40 (heure de France, 01:40 UTC) et continuer (brainstorming, voies prometteuses pour améliorer les algorithmes de Claude ou d'OpenAI, ou pour l'étape 2) ; à la réinitialisation suivante (7h40 France, 06:40 UTC), NE RIEN lancer : pause jusqu'au retour de l'auteur.** Rappels programmés à 01:42 UTC et 05:50 UTC.

**Consigne de l'utilisateur (génération 25) : faire un brainstorming d'agents ; si l'attaque principale bloque, chercher un meilleur algorithme ou une meilleure résolution au « rang inférieur », pour construire les étapes de calcul.** Fiche commune des agents : `code/gen25/BRIEF.md`.

**Consigne de l'utilisateur (24/11, soir) : continuer les tâches sans le modèle Fable** (limite Fable atteinte jusqu'au 26/11). Tous les sous-agents sont sur Opus 5.5.

**Consigne de l'utilisateur (24/11, génération 22) : attaquer de front le problème de Hurwitz et Naslund, puisque toutes les routes vers l'échelon 2 s'y réduisent.** Le coordinateur a recommandé d'arrêter et de consolider ; l'auteur a choisi de continuer.

**Consigne de l'utilisateur (24/11, génération 20) : les sous-agents utilisent le modèle Opus 5.5.** Le modèle Fable sert seulement quand il est absolument nécessaire. L'auteur laisse le coordinateur choisir les pistes vers l'échelon 2.

**Consigne de l'utilisateur (24/11) :** dès que les générations 8 et 9 sont terminées, lancer une **génération 10** sur les pistes « pas encore explorées » (section 2) et « en attente » (section 3).

## 2. Pas encore explorées

| Piste | Question précise | Difficulté | Origine |
|---|---|---|---|
| **Problème de Naslund / Hurwitz (reformulé par G18g)** | Courbe T(θ) des paires réalisables : entre la droite 2 − θ (prouvé) et le cercle √(4 − θ²) (annexe B de Claude). Cible : une borne déséquilibrée p·q ≥ N^(2+ε) en un seul θ (donne ρ > 1/2), ou une construction sous le cercle (donne ρ < 1/√2). Petit cas : P(3) = 4 ? | Très difficile (circuits non commutatifs) | N1, G5b |
| ❌ **Tests de parité économiques (Q2)** — fausse sur la famille de référence (G8b, lemme CA) | Tests avec O(n·2^O(k)) incidences au total dans l'annexe D de Claude | Technique | G6c |
| 🆕 **Valeur exacte dans le cas général à petit (ℓ, d)** | T_gen(4,1) ∈ {10, 11} (certificat exact manquant) ; (5,2) non déterminé ; relecture des théorèmes A et B′ de G11c | Moyenne | G11c |
| 🆕 **Contraction de tenseurs : N_DF ≥ 2^(cM/u)** | Étendre le plancher ρ_d ≥ 1/d aux termes intriqués ; prouvé 2^(2M/4^u), atteint 2^(O(M/√u)) ; relecture adverse des théorèmes A et B de G14b ; plusieurs niveaux de récursion | Difficile (lié à Naslund pour u = 2) | G14b |
| 🆕 **Sparsification avec variables d'extension** | Échappatoire à la proposition L de G6c | Moyenne | G6c, G7d |
| 🆕 **Lemme 4 en Lean** | Graphes aléatoires | Très difficile (peu de choses dans Mathlib) | G9L |

| 🆕 **Formes closes aux niveaux profonds (suite de G8c)** | La valeur des sous-appels de l'algorithme de Claude sur une entrée « épinglée » a une forme close calculable en temps polynomial au niveau 1 (vérifié). Si cela reste vrai aux niveaux profonds, Extended Frege réfuterait toute 3-CNF creuse insatisfaisable en taille 2^(δn). Vérifier aussi l'hypothèse (S) : degré ≤ 1 en chaque table de poids | Moyenne à difficile | G8c |


| 🆕 **Dépasser n^(2/3) dans le modèle produit** | La méthode du déterminant témoin est limitée à n^(2/3) (prop. S de G8a). Il faut une contrainte nouvelle qui couple les points ou qui utilise la structure polynomiale. Ou construire un schéma de coût n^(2/3) : il a besoin de Ω(n^(2/3)) points | Difficile | G8a |
| 🆕 **Dépasser √n pour les vrais programmes dynamiques** | Le modèle « rang d'interface » de G14c est valide et non dégénéré, mais la borne est √n (comme Dias-Lima). Il faut une distribution de témoins non produit, ou une inégalité de rang sur des coupes emboîtées. Relecture adverse de G14c à faire | Difficile | G14c |
| 🆕 **Borne de Rado pour le modèle DP corrigé** | nullité ≥ max_S(\|S\| − dim Σ U_ω(v)) contient le lemme 1 et le lemme H ; pourrait enlever l'hypothèse de régularité | Moyenne | G8a |



| 🆕 **Pièces hors du cadre « un seul espace gradué »** | Seule façon restante d'obtenir ρ < 1/√2 : des pièces qui ne passent pas par les composantes graduées d'un seul espace (ou des profils de croissance non concaves). Lié au problème de Hurwitz : P(N) ≤ N^(1+ε) | Très difficile | G9a |

| **Termes intriqués, E1 : ce qui reste après G19c** | Le plus petit cas (trois étoiles + rang r, N = 3) ne montre aucun gain : seuil = valeur produit (N(N−1)/2)³. Reste : une preuve de ce seuil ; N = 4 (≈ 10 h de calcul) ; d'autres formes (plus de trois étoiles, liaisons plus grandes, corps ℂ) ; les grilles g, h ≥ 3 (seulement un plancher ≈ 1/2). Sans idée structurelle nouvelle, piste à mettre en veille. | Très difficile | G18f, G19c |

| 🆕 **Pièces sur un corps fini fixe avec ρ < 1/√2** | Seule liberté propre à la décision. Sur GF(2) : meilleur connu 0,811 (classe restreinte) ; recherche exhaustive à λ = 2 non terminée | Difficile | G9b |

| 🆕 **Noyau de Q1\* (forme la plus précise)** | Calculer S_t = Σ_{\|T\| = τN/d} 2^(−\|N(T)\|), τ = Θ(1), en 2^(O(N/d)) sur les graphes bipartis d-réguliers aléatoires ; seule échappatoire connue : des messages à annulation | Difficile ; 3 générations sans progrès (G8b, G9c, G14f) | G14f |

| 🆕 **Fermer l'écart d'un facteur log ℓ (article 2)** | Construction non symétrique en poids qui remplit les bandes, ou preuve qu'il n'en existe pas | Moyenne à difficile | G14d |

| Contraction de tenseurs : borne de rang dans le modèle DF | Prouvé : 2^(M/4^u). Visé : 2^(M/√u). L'écart est exponentiel en u. Aussi : réalisabilité du témoin F_m par des noyaux du réseau ; incompressibilité de la somme sur β. C'est la condition pour un article 3. | G17a |
| 🆕 **Degré certifié par annulation (article 2)** | Dernière échappatoire ouverte : 1_F comme produit de o(n) facteurs creux (polynômes probabilistes) | Difficile | G14a |

**Nouvelles pistes issues de G17c (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|
| Théorème 4 avec le coût de la prépublication | Refaire le théorème du temps de calcul avec la charge réelle (2D+1)^lignes et la base de la prépublication. La borne reste-t-elle Ω(n(loglog n)²/√log n), ou tombe-t-elle à Ω(n√(loglog n/log n)) ? | G17c, Q2 et Q5 |
| Borne pour l'interpolation à plusieurs niveaux (annexe D) | Notre borne est vide pour l'annexe D. Existe-t-il une borne pour les modes 2-adiques et la récursion ? C'est la vraie question pour l'échelon 2. | G17c, Q3 |

**Nouvelles pistes issues de G18c et G18d (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|
| Noyau à profil de rang 2-adique sous-linéaire | Seul levier non constant dans la route d'OpenAI. Le noyau de Cauchy a un profil linéaire en le poids. Existe-t-il un autre noyau (non Cauchy, plusieurs variables) avec un profil sous-linéaire ? | G18d |
| Conjecture C1 et formule de Smith des blocs de Cauchy | Prouver les invariants 2^((a−i+1)+2k+2v₂(k!)) et l'égalité rang d'aplatissement = nombre de monômes de modes. Rend l'énumération optimale parmi les développements aveugles. | G18d |
| Gain constant sur l'annexe D | Modes 4 fois redondants (facteur 2), effilage optimal (1/√3) : écrire l'algorithme amélioré et vérifier le gain total (≈ 1/5 de l'exposant, non vérifié). | G18d |
| Pièces à D directions avec ρ(D) < D^(−1/2) | Loi commune des deux algorithmes : gain = D^(−1/2). Le plancher 1/D donnerait n/log n. | G18c |

**Nouvelles pistes issues de G18h (24/11) — les trois objets qui passeraient sous 1/√2 :**

| Piste | Question précise | Origine |
|---|---|---|
| Profil en escalier (après G19b) | Tel que défini par G18h : exclu par le théorème 2. Réparé : équivaut à « L linéaire sans (V) », c'est-à-dire à toute la droite. Objet à construire ou exclure : pour tout θ, un sous-espace de M_N de codimension N^(2−θ+o(1)) où la transposition est une somme de N^(θ+o(1)) multiplications des deux côtés. Calcul suivant : b(N) pour N = 9 à 12. | G18h, G19b |
| Opérateurs presque commutants à croissance exponentielle | Familles N×N sur un espace gradué avec log dim ~ E (pas √E) et une fenêtre H = o(E₀). Sans l'hypothèse (V), le compte donne ρ ≥ 2^(−α). | G18h |
| Identité de localité à D directions | Une identité sur un espace avec log dim ~ E^((D−1)/D). Donnerait ρ(D) ≈ D^(−(D−1)/D) (spéculatif). Petits cas : D = 3 : 1/2 ; D = 4 : 3/8. | G18h, G18c |

**Nouvelles pistes issues de G18S (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|
| **Le « paquet en escalier » (reste de la piste des charges non uniformes, G19a)** | Étiquetage circulaire : u faisceaux entre classes antipodales sur le 2u-cycle. Encadrement prouvé : φ′M/(2u) − u ≤ log₂ N ≤ min_i M_i. Quel est N ? u = 2 est le problème de Hurwitz ; u = 3 à 3 fils : N = 4 dans le sous-modèle GF(2). Une identité à 2^(o(M)) termes sur un expanseur serait l'échelon 2 ; une borne 2^Ω(M) fermerait la piste. Il faut une méthode autre que l'aplatissement. | G19a |

**Nouvelles pistes issues de G20d (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|
| Borne inférieure au-delà de la droite (après G21d) | Pas par des équations de bas degré ni par des aplatissements (exclus jusqu'au degré 7 à N = 3). Routes restantes : une quantité analytique (distance exacte 2 à N = 3 ; structure des points les plus proches) ; une version « au bord » des arguments spectral et de Segre. Aucune n'a de premier pas concret. | G20d, G21d |
| Barrière pour les aplatissements de Young | G21d : à N = 3, aucun des 168 aplatissements de Koszul/Young ne donne une équation de J(2,2). Reste : un théorème général (étendre le théorème B aux Φ linéaires générales). | G20d, G21d |

**Nouvelle cible issue de G20c (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|
| **Profil à plateau (après G21e et G21a)** | Les contraintes de dimension permettent ρ → 1/2 (G21e), mais l'objet qu'il faut est un seul module cyclique avec dim ≈ N^(1+o(1)) sur [E0, E0 + H] : aussi fort que P(N) = N^(1+o(1)) (problème de Hurwitz-Naslund). Ni construit, ni exclu. Restes non testés : familles non unitaires avec c < 0 ; étiquettes de poids mixtes ; algèbres non libres (W-algèbres, affines à niveau fractionnaire). | G21e, G21a |

**Pistes issues de G22e (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|

**Piste issue de G22d (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|
| Masse captée c(p, q) (après G23d) | Prouvé : c(N;1,1) ≤ 8,081 et c(N;p,1) ≤ 11,7p + O(√p) pour tout N. Reste : la valeur exacte 4 (reformulée en une inégalité matricielle explicite) ; et surtout q ≥ 2, où la preuve ne s'étend pas (les vecteurs singuliers ne sont plus des produits). Sans une borne sous-linéaire en p·q au point équilibré, pas de gain d'exposant. | G22d, G23d |

**Pistes issues de G22b-bis (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|
| Modèle monomial après G23a | Les ansatz de groupe sont exclus (obstruction d'anisotropie). Reste : tableaux ST génériques sous N²/8 ((28,28) à N = 16 non décidé) ; ZT avec τ non linéaire à N = 64. Exposant du modèle : entre 4/3 (prouvé, Lean) et 2 (meilleur connu). | G22b-bis, G23a |
| Relecture adverse de la borne du modèle monomial | Lemme clé (une couleur T de taille h force |S| ≥ h(h−1)/2) ; comptage ; usage du lemme 7.4 de HWY. | G22b-bis |
| Décider N = 12 (18,18) et N = 16 (28,28) dans le modèle ST | SAT arrêté à 600 s. | G22b-bis |

**Piste issue de G23c (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|

**Question décisive issue de G24c (24/11) :**

| Piste | Question précise | Origine |
|---|---|---|

## 3. En attente

| Piste | Raison de l'attente |
|---|---|
| ⏸️ **Version 3 de l'article 1** | **EN ATTENTE — décision de l'auteur (24/11) : ne pas faire la version 3 avant le retour de F. Dias**, à qui la version 2 a été envoyée pour vérification. Contenu préparé pour la version 3 (à compléter avec ses remarques) : (1) section Lean mise à jour : théorème principal prouvé sans hypothèse (G10e) ; (2) corrections de G10f : facteur b superflu, dire que la borne de Dias-Lima est meilleure pour les petits groupes, précisions de l'annexe A, 4 attributions à adoucir ; (3) théorème sep(F\*) ≥ c·n\* et son corollaire, si la relecture G11a le confirme ; (4) nouveaux chemins Lean (`lean/PnpLean/Article1/`) **Ajout G17b :** 7 corrections tirées de la prépublication elle-même (`code/gen17/g17b/rapport.md`). |
| ⏸️ **Rédaction de l'article 2** (erreurs par point) | Condition de l'auteur : un résultat final. **État (24/11) : théorème prouvé en Lean sur tout corps, confirmé par relecture adverse, nouveauté : semble nouveau (technique classique).** Réserves à écrire : vaut dans un modèle de coût ; ne couvre pas la contraction de tenseurs ni D < nd ni L/R′ ni les multiplicités inégales (G14a, G14d en cours). **En attente de l'accord de l'auteur pour rédiger.** |
| ⏸️ **Note « plus long chemin / plus long cycle en 2^o(n) sur les graphes creux »** (candidat article 4) | Preuve confirmée par G10c ; conditionnel aux prépublications ; doit citer Wu-Wang (arXiv 2611.02548) ; résultat modeste. Décision de l'auteur |
| ⏸️ **Licence du dépôt** | Sans licence, personne n'a le droit de réutiliser la bibliothèque `PnpLean.Outils`. Décision de l'auteur (Mathlib utilise Apache 2.0). |
| ⏸️ **Relecture de l'article 1 par un expert humain** | **En cours : l'auteur a envoyé la version 2 à F. Dias (24/11).** Attendre sa réponse. |

| Piste | Raison de l'attente |
|---|---|

## 4. Terminées (résultats)

| Piste | Résultat |
|---|---|
| ✅ **Constante de Claude : Freedman et plancher** (G26-P5) | c1 4,05 → 3,82 (R2F + R3″) ; structure S2 : ≈ 3,31 après correction par P10 (esquisse). Seuil contre 2^n : 1 943, puis 1 514. Plancher des lemmes uniformes ≈ 3,36 ; α ≥ 1 pour toute conception. Relu (P+, sauf S2). `code/gen26/p5/notes.md` |
| ✅ **Relecture courte des lemmes de G26-P4** (G26-P16) | Quatre lemmes et relaxation GF(2) confirmés (P+) ; aucune violation sur 3 600 instances. |
| ⚠️ **Recherche SAT exacte dans le modèle MC** (G26-P4) | (9,3)@5 et (6,4)@5 impossibles dans MC (en partie kissat seul) ; N = 4 complet avec DRAT ; (12,3)@6 retrouvé ; (5,5)@5 et (11,3)@6 non décidés. Lemmes implicites relus par P16 (P+). `code/gen26/p4/notes.md` |
| ✅ **Relecture de G26-P13** (G26-P14) | E, C, D confirmés ; B corrigé (terme de concentration aux bases) : c1 ≈ 3,18 à 3,21 (P+), seuil ≈ 1 370 à 1 400. `code/gen26/p14/notes.md` |
| ✅ **Vérification des 9 points de G26-P11** (G26-P15) | Tous confirmés (rang exact sur ℚ, deux méthodes) ; tous nouveaux et non dominés. `code/gen26/p15/notes.md` |
| ❌ **Conception spectrale de constructions** (G26-P11) | Aucun point au-delà de la droite. Optimums : spectre {1, −1, −2}, blocs-cliques (p = nombre de blocs). 9 nouveaux points exacts (frontière N = 7 complète). Modèle à blocs-cliques : T ≥ (5 − 2θ)/3. `code/gen26/p11/notes.md` |
| ✅ **S2 rigoureux, découpage, covariance** (G26-P13) | c1 = 3,312 entièrement compté ; découper + covariance (iii) : **c1 = 3,166**, seuil 1 352 contre 2^n à densité 4,26. Relu par P14 : c1 ≈ 3,18 à 3,21 après correction. `code/gen26/p13/notes.md` |
| ⚠️ **Relecture de G26-P8 ; attaque de (A)** (G26-P12) | Théorèmes U et A2 confirmés (P+). (A) réfuté (a = 6, k = 2). C1′ faux pour des sous-blocs ; ouvert pour les blocs complets (a ≥ 5). Forme faible de (A) pour k = 2 prouvée. `code/gen26/p12/notes.md` |
| ❌ **Enveloppe de l'annexe B à grand h** (G26-P1) | Dans l'ansatz structuré (degré 1, D = h), k maximal = 3 pour h = 4 à 11 ; loi des fonctionnelles prédite puis confirmée à h = 10, 11. Gain constant, pas linéaire : θ ne change pas. Reste ouvert : états de poids ≫ h. `code/gen26/p1/notes.md` |
| ✅ **Seuils d'OpenAI avec les modes de Smith ; C1′** (G26-P8) | Seuils ligne H : 2,25e5 (contre 2^n) et 1,52e6 (contre PPSZ) à densité 4,26 ; constante → 37,5. Théorème U : C1′ ⟸ (A) ; théorème A2 : C1′ au niveau a − 2. U et A2 relus (P+) ; (A) réfuté par P12. `code/gen26/p8/notes.md` |
| ✅ **Relecture de G26-P5** (G26-P10) | R2F, R3″, c1 = 3,82 et seuils confirmés (P+). Erreur d'un facteur 2 dans S2 : corrigé c1 ≈ 3,31, seuil 1 514 (esquisse). Planchers : avec réserves. `code/gen26/p10/notes.md` |
| ✅ **Relecture de G26-P3** (G26-P9) | Aucune erreur qui casse les résultats : P3 passe au niveau P+. Réserves : matrices signées ±1, caractéristique 0, bornes vides pour N ≲ 20. `code/gen26/p9/notes.md` |
| ⚠️ **Au-delà de la droite dans le modèle monomial-C** (G26-P3) | Pas de gain d'exposant. Nouveau : p ≥ (1/2 − o(1))·N² si q = o(log N) ; p ≥ c·N²·ln N/q ; théorème de groupe ; gain d'exposant conditionnel (énergie des mots de longueur 4) ; la méthode des polynômes est limitée à N²·log N/q. Relu (P+). `code/gen26/p3/notes.md` |
| ✅ **Relecture de la formule de Smith** (G26-P7) | Théorèmes 1 et 2 et réfutation de C1 confirmés (P+). Constante d'OpenAI 37,5 rigoureuse relativement à l'annexe D et aux lignes A-G de P3. `code/gen26/p7/notes.md` |
| ✅ **Formule de Smith des blocs de Cauchy** (G26-P6) | Prouvée (non relue). Constante d'OpenAI 75 → 37,5 rigoureuse (relativement à l'annexe D). La conjecture C1 de G18d est fausse telle qu'écrite ; C1′ proposée. `code/gen26/p6/notes.md` |
| ✅ **Relecture de B7** (G26-P2) | R1 à R4 confirmés (P+). Constante de Claude 4,05 (constantes du relecteur). Seuil contre la force brute à densité 4,26 : log₂ m* = 2 237 (papier : 276 000). Mémoire polynomiale : 6,16, seuil 5 123. `code/gen26/p2/notes.md` |
| ⚠️ **Arbre des sous-problèmes** (G26-B5) | Nouveau modèle monomial-C qui contient tous les optimums exacts ; lemme MC2 : p·q ≳ ½·N²·ln q (au-delà de la droite d'un log). Pas de gain linéaire à q = 3 (numérique, k ≤ 5). `code/gen26/b5/notes.md` |
| ✅ **Rang inférieur, deuxième tour** (G26-B7) | Constante de Claude c1 = 4,12 (déterministe et Las Vegas ; papier : 42,4) ; mémoire polynomiale : 6,16 ; compromis temps-mémoire (1 − M)/3. Non relu. `code/gen26/b7/notes.md` |
| ❌ **Perturbations de l'annexe B** (G26-B6) | Rien ne courbe la frontière. Lemme 1-D : le seuillage est optimal pour une paire ; un gain demande des relations entre paires. Enveloppe exacte jusqu'à N = 15 : sur la droite triviale. `code/gen26/b6/notes.md` |
| ✅ **Synthèse des résultats** (G25-P10) | Encadrement pour l'article 3 (écart constant 22,7 avec P3) ; conséquences explicites (Max-Cut, Max-3-SAT, #CSP, plus long chemin : constante 5,3) ; **la fenêtre du module de classe se réduit à Hurwitz-Naslund** ; « pinceau + petit rang ». `code/gen25/p10/notes.md` |
| ✅ **Relecture du théorème A explicite** (G25-P8) | Tout confirmé (P+). Le balayage est un ordre d'élimination optimal. Réserve : mémoire exponentielle (≈ 2^(0,95 × temps)). `code/gen25/p8/notes.md` |
| ⚠️ **Balayage à chaque niveau** (G25-P9) | Aucun levier ne change la puissance du logarithme : loi de Pythagore 1/e′² ≤ 1/e_a² + 1/e_b² pour tout arbre de niveaux (Prop. P). Changer la puissance demande ρ < 1/√2 (Hurwitz-Naslund encore) ou une feuille avec e0 = o(1/√log K0). Constante : 7,47 → 4,82 (Las Vegas), 10,6 → 6,81 (déterministe). `code/gen25/p9/notes.md` |
| ⚠️ **Conjecture d'étalement** (G25-P7) | Ni réfutée ni soutenue. L'annexe B est extrêmement non étalée (annulations énormes). L'étalement serait au mieux une propriété des points. Une version limitée à q ≤ N^β reste à tester. `code/gen25/p7/notes.md` |
| ✅ **Annexe D d'OpenAI explicite** (G25-P3) | ε légal maximal 1/960 (papier : 2⁻¹⁷, 136 fois trop petit). Exposant E ≤ 75·k·m/√log₂ s (rigoureux, densité bornée) ; 37,5 avec une conjecture. Seuil contre 2^n : 7·10¹¹ (papier) → 9·10⁵. Corrections au dossier : b4 sous-estimait le papier d'un facteur 20 ; le gain « effilage LP 1/√3 » venait seulement de la constante 3. `code/gen25/p3/notes.md` |
| ❌ **Triangularisation pour la fenêtre du module de classe** (G25-P4) | Non : l'hypothèse du lemme T échoue pour de vraies solutions (certificats exacts). La route de preuve tombe, pas la borne. Fenêtre E0 ≥ 2 toujours ouverte. `code/gen25/p4/notes.md` |
| ✅ **Théorème A de Claude explicite** (G25-P2) | Constante de l'exposant 42,4 → 10,6 (déterministe) / 7,47 (Las Vegas). Lemme de balayage prouvé et vérifié. Seuil contre la force brute ÷ 29 à 58 (densité 4,26 : log₂ m* ≈ 9 375 ; densité 1 : 151). Toujours sans effet aux tailles réalistes. `code/gen25/p2/notes.md` |
| ⚠️ **Sparsification pour la décision** (G25-P1) | Pas de construction. Sans restriction, la question équivaut au problème lui-même (Prop. 0). Pour les sparsificateurs modulaires « restriction + sous-formule », le facteur loglog n reste (corollaire de L′). Ouvert : arbres avec autarcies. `code/gen25/p1/notes.md` |
| ❌ **Hypercliques de novembre 2026** (G25-P5) | Pas d'identité à 3 directions : Fischer = croisement D = 2 + multiplication rapide de matrices ; Chan-Itakura-Ma-Yu = construction de Hrubeš + symétrie du 4-cycle. Max-3-SAT en 2^(0,99n) (leur résultat). `code/gen25/p5/notes.md` |
| ✅ **Lean : lemme T et « 2 produits = 1 »** (G25-P6) | `lean/PnpLean/Outils/Triang.lean` : lemme T, lemme P, f(N,1) = f(N,2) = N(N−1)/2 sur tout corps. Niveau L. |
| ✅ **Lemme d'échange fort** (G24g) | **Théorème 1 : f(N,2) = f(N,1) = N(N−1)/2 pour tout N, tout corps** (deux produits = un seul). Lemme de triangularisation T. r ≥ 3 ouvert. Critère pour fermer la fenêtre du module de classe : la famille de G24f se triangularise-t-elle ? `code/gen24/g24g/notes.md` |
| ✅ **Petits cas par SAT** (G24d) | Nouveaux UNSAT avec preuves DRAT : ST (7,7) à N = 7, ST (10,10) à N = 9, GF(2) (3,3) et (5,2) à N = 4 (pas de gain en caractéristique 2). Encadrements ST pour N = 9 à 16. Grandes tailles non décidées. `code/gen24/g24d/notes.md` |
| ❌ **(N, N) réalisable ?** (G24e) | Non au-delà de 8 (numérique) : (7,7) à N = 7 non réalisable ; saut net après N = 8 (motif de Hurwitz). (5,5) à N = 5 et (6,6) à N = 6 : seulement au bord. Bornes exactes inchangées. `code/gen24/g24e/notes.md` |
| ⚠️ **Fenêtre du module de classe** (G24f) | Aucune construction sous 2√2. Bande E0 = 1 fermée par une preuve (4ρ ≥ 4 ; lemme d'échange = théorème 3 de G19b). Reste ouvert : E0 ≥ 2 avec H < 0,83·E0 ; il faudrait le lemme d'échange avec un terme de rang de Schmidt r (c'est notre théorème C dans une autre forme). Heuristique : toute la fenêtre à 4ρ ≥ 4. `code/gen24/g24f/notes.md` |
| ⚠️ **Plusieurs canaux, algèbres de sommets** (G24b) | Fermé pour A1^n et (A1+A1)^n avec étiquettes dans V′ (nouvelle borne N² ≤ rang F). Nouveau cas ouvert : étiquettes dans le module de classe V_{α/2+L}, borne 2,674 < 2√2, aucune construction. `code/gen24/g24b/notes.md` |
| ⚠️ **Distance pour q ≥ 2** (G24c) | Borne indépendante de N seulement pour une sous-classe (c ≤ 13p + 4q) ; lemme L2 manquant. Données : (6,6) réalisable numériquement à N = 6 (b(6) = 6 ?). Question décisive : (N, N) reste-t-il réalisable ? Si oui, la droite est la vérité et ρ ≤ 0,683. `code/gen24/g24c/notes.md` |
| ❌ **Pièces de Hrubeš, lecture non symétrique** (G24a) | Courbe exacte θ = 1 + 2H(b), T = 1 + 2H(0,110 − b) : strictement hors du cercle partout. Limites prouvées du cadre de Clifford : θ, T ≥ 1 et max(θ, T) ≥ 1,5. L'hybride avec l'annexe B ne gagne rien. `code/gen24/g24a/notes.md` |
| ❌ **Ansatz de groupe dans le modèle monomial** (G23a) | Aucune famille invariante par un groupe ne passe sous la somme directe N²/8 (N = 16, 27, 32, 64). Obstruction prouvée (non relue) : l'ansatz ZT force une forme anisotrope, donc des blocs de taille ≤ 4 (Chevalley-Warning). Le point Q8 est un accident de la taille 8. `code/gen23/g23a/notes.md` |
| ⚠️ **Loi de distance pour q = 1** (G23d) | c(N;1,1) = 4 non prouvé ; prouvé : c(N;1,1) ≤ 8,081 pour tout N, et c(N;p,1) ≤ 11,7·p + O(√p) (indépendant de N). Rien pour q ≥ 2. `code/gen23/g23d/notes.md` |
| ✅ **Lean : borne du modèle monomial** (G23e) | `lean/PnpLean/Outils/Intercalate.lean` : lemme clé, comptage, compromis réel. Niveau L (modèle ST). |
| ❌ **Creux des algèbres de sommets, vrai rang** (G23c) | Aucune valeur sous 2√2 avec le vrai rang. Règle de sélection (ξ intérieur) et dichotomie (ξ sur le bord) prouvées sur papier : elles tuent les creux pour toute algèbre fixe. Reste ouvert : constructions à plusieurs canaux dans des familles de puissances tensorielles (écart 2,799 – 3,030). `code/gen23/g23c/notes.md` |
| ✅ **Relecture du modèle monomial** (G23b) | Tout confirmé (P+) ; hypothèse inutile retirée ; énumérations exhaustives sans échec ; les optimums exacts à N = 5 sont hors du modèle. `code/gen23/g23b/notes.md` |
| ✅ **Modèle monomial** (G22b, G22b-bis) | Correspondance exacte avec les matrices d'intercalates de HWY (formules entières de sommes de carrés). **Borne dans le modèle au-delà de la droite : p²·q ≳ N⁴/8 ; point équilibré ≥ N^(4/3)/2** (papier, non relu). Meilleur connu dans le modèle : exposant 2. Le modèle n'est pas exclu pour battre √2 (4/3 < √2). `code/gen22/g22b/notes.md` |
| ⚠️ **Distance de l'échange à J(p, q)** (G22d) | La distance (masse captée c = N² − d², toujours entière) dépasse les méthodes de rang à tous les points calculés ; loi exacte pour q = 1 : d² = N² − m² − 2j. Théorème C quantitatif prouvé (faible). Pour gagner un exposant il faudrait c(p, q) sous-linéaire en p·q : non prouvé ; q = 1 est linéaire. `code/gen22/g22d/notes.md` |
| ⚠️ **Barrière de Cardy** (G22e) | Théorème A (papier, non relu) : √2 est optimal pour toute algèbre de sommets fixe avec croissance √E. Pour les familles : ouvert. Échappement possible seulement avec un « creux » de L(E)/√E et un nombre de paires aussi grand que la borne lâche. Deux routes (regraduation du réseau A1 ; modules de sl₂) donnent 2,81 < 2√2 sous la borne lâche, mais 2,90 à 3,03 avec les vraies bornes. `code/gen22/g22e/notes.md` |
| ✅ **Littérature Hurwitz-Naslund lue** (G22a) | Notre problème = problème 12 de Naslund = question 1.2 de Karam = largeur de HWY = PT-rank de Rossman-Zhu. Record : √2 (annexe B) ; publié : 1,62 (Hrubeš). **Le théorème C n'est dans aucun article lu ; la question du compromis (p, q) est nouvelle.** Techniques non essayées : comptage d'intercalates de HWY, pièces de Hrubeš, caractéristique 2. `code/gen22/g22a/notes.md` |
| ❌ **Familles structurées jusqu'à N = 256** (G22c) | Aucune ne bat l'exposant 2 ni la somme directe (N²/8). Sédénions : défaut de rang Θ(N²). Groupes, Paley, tournois : N² − O(N). `code/gen22/g22c/notes.md` |
| ❌ **Équation au-delà des méthodes de rang** (G21d) | Aucun témoin. À N = 3 : aucune équation de J(2,2) en degré ≤ 6, aucun aplatissement de Koszul/Young (168 testés), aucune équation séparante en degré 7. Séparer l'échange est un phénomène de degré ≥ 8. Distance² de l'échange à J(2,2) = 2 (numérique). Reformulation symétrique à trois aplatissements. `code/gen21/g21d/notes.md` |
| ❌ **Caractéristique p** (G21b) | Rien de nouveau. La récursion de Claude tourne bien sur F_{p^k} pour la décision (deux étapes réparées), mais l'annexe B sur F_p donne le même cercle (le facteur 1 − 1/p se simplifie). Frobenius : aucun gain. Dvir-Edelman ne s'applique pas au croisement. Aucun petit cas décidé n'est meilleur sur F_p. Seul reste : (3,3) à N = 4 sur GF(2), non décidé. `code/gen21/g21b/notes.md` |
| ❌ **Modules à espace fondamental dégénéré** (G21a, modèle Fable) | Aucune construction. Un espace fondamental est une perte pure (ρ ≥ 1/2 + L_M(0)/(4λ)). Neuf candidats (Fock, tordus, NS, Ramond, symplectiques), tous au-dessus du cercle ; Ramond : 1,4279 contre 1,4142. Le lemme de localité s'étend aux modules ; c'est le profil qui échoue. Un plateau demande un module cyclique aussi fort que P(N) = N^(1+o(1)). Précision : le temps est 2^(O(m/(log m)^α′)) avec α′ = log₂(1/ρ). `code/gen21/g21a/notes.md` |
| ✅ **Relecture du théorème B et de la borne 2/3** (G21e) | Théorème B confirmé (P+), deux corrections de portée (compositions non locales hors classe ; « exactement 1/4 » seulement comme barrière). Borne ρ ≥ 2/3 : vraie seulement avec la concavité, qui est une hypothèse. **Sans concavité, les contraintes de dimension permettent ρ → 1/2** (profil rampe – plateau – rampe). Question décisive : les lemmes 9 et 10 du papier de Claude valent-ils pour des modules à grand espace de plus bas poids ? `code/gen21/g21e/notes.md` |
| ❌ **L'échelon 2 par d'autres problèmes ouverts** (G21c) | Aucune route. Fermé : conjecture du rang asymptotique et ω = 2 (théorème R : constante seulement) ; le log de Marx est une quantité de routage sans lien avec l'échelon 2 (T2, T3) ; k-Clique est le pire cas des algorithmes de 2026. Positif conditionnel : CSP en s^(O(tw/√log tw)) sur les expanseurs de degré borné (T1, esquisse). Correction : les bornes de 2026 sont en m ; en n elles ont un facteur loglog n. `code/gen21/g21c/notes.md` |
| ⚠️ **Pièces à D = 3 directions** (G20a) | Aucun signe que D = 3 bat D^(−1/2) au-delà du plus petit cas : exposants 0,40 → 0,50 → 0,65 pour N = 2, 3, 4 (sphère : 0,577 ; plancher : 0,333). La frontière des triplets est nécessaire mais pas suffisante. Théorème G (généralise le théorème C, non relu) : exposant 1/4 seulement au point équilibré. Solveur non fiable à l'intérieur. `code/gen20/g20a/notes.md` |
| ⚠️ **Annexe B à N fini et frontière b(N) jusqu'à N = 12** (G20b) | L'annexe B de Claude est triviale (b = N²) pour N ≤ 256 ; son exposant passe sous 1,6 vers N ≈ 6·10¹¹. Les petits cas ne disent donc rien sur le cercle. b(N) = N + O(1) jusqu'à N = 8, saut à N = 9 (motif de Hurwitz). Aucune famille d'exposant < √2. Non tranché. `code/gen20/g20b/notes.md` |
| ❌ **Constructions algébriques sous le cercle** (G20c) | Aucune. Sept familles testées, toutes au-dessus du cercle ou dominées. La non-rigidité d'Alman-Williams ne donne rien (la partie C coûte la parcimonie de Fourier, pas le rang). Dans le cadre gradué avec profil concave : ρ ≥ 2/3, atteint seulement par un profil « à coude » (N, N^(4/3), N²) qu'aucune algèbre connue ne réalise. `code/gen20/g20c/notes.md` |
| ❌ **Pièces adaptées à l'instance** (G20e) | Pas de gain au-delà d'une constante sur les expanseurs. Le croisement restreint aux états atteignables est encore un échange complet (proposition R), donc le théorème C s'applique. Conductance mesurée ≥ 0,115. Brancher d'abord ne gagne rien. `code/gen20/g20e/notes.md` |
| ⚠️ **Barrière pour les méthodes de rang** (G20d) | Théorème B (papier, non relu) : aucune méthode de rang de classe B (aplatissements suivis de multiplications) n'exclut une paire avec p·q ≥ N²/4. Le théorème C est donc optimal pour sa classe. L'échange est spécial, pas générique (frontière générique : p + q ≈ 0,59·N²). Non couvert : aplatissements de Young. `code/gen20/g20d/notes.md` |
| ✅ **Lean : théorème C** (G20g) | `lean/PnpLean/Outils/SwapRank.lean` : parties (a), (b), (c), formes miroirs et forme produit 2·Q·M² ≤ 9·r·b, tout corps. Niveau L. |
| ✅ **Relecture adverse de G19a et G19b** (G20f) | Dix affirmations confirmées (niveau P+). Théorème E : preuve complète. Théorème 5 : faux tel qu'écrit (q = 1, N = 2, k = 3), corrigé avec 1 ≤ k ≤ N. Les quatre décompositions nouvelles sont des identités sur ℤ. Cinq candidats pour Lean. `code/gen20/g20f/notes.md` |
| ⚠️ **Paires réalisables exactes, profil en escalier** (G19b) | Rien sous le cercle, rien au-dessus de la droite dans l'exposant. Le profil en escalier de G18h tel que défini est exclu ; réparé, il couvre toute la droite. Frontière exacte pour N ≤ 3, partielle pour N = 4 à 6 ; P(3) = 4, P(4) ≥ 6. Théorèmes 3, 4, 5 : gains de constante seulement. Droite contre cercle : indécidable à ces tailles. Les petits gains ne se multiplient pas. `code/gen19/g19b/notes.md` |
| ❌ **Termes intriqués, le coin ouvert à N = 3** (G19c) | Aucun gain (numérique, étalonné). Seuil r* = 27 = (N(N−1)/2)³ : la valeur du point produit. r = 26 : 0 succès sur 6. Le « gain » à N = 2 est un accident. Même seuil produit sur trois cas à dimensions mixtes. `code/gen19/g19c/notes.md` |
| ✅ **Lean : théorème A″** (G19d) | `lean/PnpLean/Article3/DFRank.lean` : A″ sous trois formes, plus la forme modulo un idéal propre. Niveau L. Reste sur papier : lemme de restriction, constante 1/√π, lien avec l'annexe D. |
| ✅ **Relecture et antériorité du théorème F et du théorème 2** (G19e) | Les deux sont justes et sont un seul théorème (théorème C, échange avec spectateur) : r ≥ kQ(n − kb), constante exacte 2/9, plus (c) b < mn ⇒ r ≥ Q (nouveau). Flanders inutile. Quatre errata dans G18f. P(3) = 4 sur F₂. Antériorité : problème 12 de Naslund ; cas k = 1 connu (Karam) ; probablement du niveau « folklore ». Ce théorème ne peut jamais passer la droite. `code/gen19/g19e/notes.md` |
| ⚠️ **Charges non uniformes sur une instance difficile** (G19a) | L'ouverture « Hamming 1 » est fermée sur les expanseurs (inégalité de séparateur : coût ≥ 2^(φ′M)). Théorème E : log₂ N ≥ min_i M_i/(2u) − u pour tout étiquetage. Mais l'aplatissement ne peut pas donner M/√u en général (étiquetage circulaire : A″ = M/u exactement). Ouvert : le « paquet en escalier ». `code/gen19/g19a/notes.md` |
| ✅ **Borne de rang dans le modèle DF** (G18e, relu par G18S) | Théorème A′ confirmé (P+), forme plus forte A″ : N ≥ 2^((0,564−o(1))·M/√u) pour des charges presque uniformes. Avec la borne supérieure de l'annexe D : 2^Θ(M/√u). « √log est optimal » est un théorème (côté inférieur) dans une famille précise. **Faux pour des charges non uniformes : N = 2^(M/u) est alors atteint.** Énoncé prêt pour Lean. `code/gen18/g18e/`, `code/gen18/g18s/` |
| ⚠️ **Termes intriqués, question E1** (G18f) | Non tranchée. Théorème F (papier, non relu) : r ≥ s·Q·(M − b·s) ; négatif pour un croisement (ferme le point ouvert de G14e) et pour 2×3 hors d'un coin. E1 n'est pas une question de rang asymptotique. Le calcul λ = 2 n'est pas fait (jours de calcul). Plus petit cas ouvert : P_N = (rang ≤ N^(2−δ)) + trois étoiles. `code/gen18/g18f/notes.md` |
| ❌ **Réfutation des 3-CNF aléatoires** (G18i) | L'obstacle central ne bouge pas. Carte temps contre densité mise à jour : les algorithmes de 2026 sont les meilleurs seulement pour Δ ≤ c·log n (2^(O(n/√(Δ log n)))). Énoncé candidat : réfutations IPS de taille 2^(O(km/√log m)) (une étape non contrôlée). Rien pour Extended Frege. `code/gen18/g18i/notes.md` |
| ✅ **Évaluateur modulo 3 pour la route d'interpolation** (G18b, relu par G18R) | Théorème conditionnel confirmé (P+) : décision aléatoire et #F mod 3 pour les 3-CNF à occurrences bornées en 2^((12,6+o(1))·n·lnln n/√ln n), un facteur loglog sous la route d'interpolation de la prépublication. Dominé par le résultat principal de la prépublication : résultat sur la méthode seulement. Complète l'article 2 (borne supérieure loglog contre borne inférieure √loglog pour une charge constante). `code/gen18/g18b/`, `code/gen18/g18r/` |
| ⚠️ **Problème de Hurwitz et Naslund** (G18g) | Aucune borne améliorée. Reformulation : la variante utile est « paires réalisables (p, q) », pas les sommes de carrés ; 4ρ* = sup_θ (θ + T(θ)). Borne élémentaire p·q ≥ 3N²/16. Notre ancienne « équivalence » avec P(N) ≥ N^(1+ε) est corrigée : un plancher ρ > 1/2 suit d'une borne déséquilibrée p·q ≥ N^(2+ε) en un seul θ. `code/gen18/g18g/notes.md` |
| ⚠️ **Pièces hors du cadre** (G18h) | Aucune famille avec ρ < 1/√2. Mais le plancher 1/√2 vient de l'hypothèse (V) et de la forme du profil de croissance, pas de « un seul espace gradué ». Seul plancher prouvé en général : 1/2. Trois objets donneraient un gain (opérateurs presque commutants à croissance exponentielle ; profil en escalier ; identité à D directions). `code/gen18/g18h/notes.md` |
| ❌ **Précision 2-adique des vrais noyaux** (G18a) | Pas de gain au-delà d'une constante. Rang plein ⇒ même précision que la boîte noire (17 cas). Le « 4 contre 14 » de G17a est un effet de nœud de coupe sur peu de paquets. Constante de l'exposant : facteur 2 à 3,5 possible (κ = 24 → 2 à 5,9). Numérique seulement. `code/gen18/g18a/notes.md` |
| ⚠️ **Limite de l'interpolation récursive** (G18d) | n/√log n dans toutes les familles analysées. La récursion ne compose pas les gains : la racine domine. La profondeur achète seulement √loglog. Gain constant possible (modes 4 fois redondants). Borne inférieure seulement « au niveau de la conception » ; famille générale : ouvert. Conjecture C1 (formule de Smith des blocs de Cauchy). `code/gen18/g18d/notes.md` |
| ❌ **Hybride des deux algorithmes de 2026** (G18c) | Aucun gain sur six substitutions. Les deux algorithmes suivent la même loi heuristique : gain = (nombre de directions)^(−1/2). Un vrai gain demande des pièces à D directions avec ρ(D) < D^(−1/2) (plancher 1/D → échelon 2). `code/gen18/g18c/notes.md` |
| ❌ **Les deux points ouverts de l'article 2** (G17e) | Aucun n'est fermé. Écart log ℓ : ouvert dans la fenêtre ; optimum à 1,1–2,3 fois la borne sur de petits cas ; l'argument à deux moitiés n'aide pas ; nouveaux lemmes (dualité de quadrature, élagage). Degré par annulation : générique mais pas forcé ; un exemple d'annulation existe. `code/gen17/g17e/notes.md` |
| ⚠️ **Relecture adverse de G14b** (G17a) | Reconstruction de l'annexe D fidèle. Théorème A confirmé mais très faible. Théorème B en partie : axiome du modèle VT réfuté. Nouveau lemme : précision Θ(M·D_t/u) optimale pour la reconstruction boîte noire. Potentiel « réel mais modeste » : pas une barrière. `code/gen17/g17a/rapport.md` |
| ✅ **Relecture adverse de G14a + Lean** (G17d) | Théorème au niveau de la sortie confirmé et formalisé (`Article2/OutputLevel.lean`) : rang de A(0)[budget] ≤ Σ rang E_k, tout corps, sans factorisation. L'hypothèse « mineur non singulier » est nécessaire en général (contre-exemple), mais se retire à coût 1/(2Δ−1) pour les matrices creuses (papier). Degré réduit confirmé. |
| ⚠️ **Article 2 contre la prépublication d'OpenAI** (G17c) | Le cadre existe dans la prépublication (après traduction ; vérifié exactement pour Q ≤ 6). Nouveauté confirmée. **Mais** le coût de l'erreur diffère (tuples de mesures, pas n^rang), et l'annexe D échappe au cadre (modes 2-adiques, récursion). L'énoncé « perte optimale » vaut seulement dans le « cadre X » à un niveau sur un corps : `code/gen17/g17c/rapport.md` |
| ✅ **Article 1 contre la prépublication d'OpenAI** (G17b) | Aucune hypothèse contredite. Règle de copie = §2.8 ; elle s'applique aux feuilles de la sparsification, pas à F. Définition 5.1 = généralisation fidèle. Coût réel = type produit. Le plan réel n'est pas uniforme en la variable. Sur le plan réel, le théorème B est vrai par q seul : il contraint les nouveaux plans à petits groupes. 7 corrections pour la v3 : `code/gen17/g17b/rapport.md` |
| ✅ **Bibliothèque d'outils réutilisables** | `lean/PnpLean/Outils.lean` : 13 énoncés génériques vérifiés par Lean et absents de Mathlib. Reste : généraliser de ℚ à tout corps, remplacer `uavg` par `Finset.expect`, licence (décision de l'auteur). Voir `lean/OUTILS.md` |
| ✅ **Lean : théorème à plusieurs points** (G16L) | `MultiPoint.sum_rank_ge` : forme la plus forte du théorème de l'article 2 (borne doublée, tout corps, placement quelconque) |
| ✅ Modèle DP valide (G14c) | Modèle « rang d'interface » : valide pour les DP à tables d'interface, non dégénéré ; borne 2^Ω(√n) seulement (non relu) |
| ✅ Échappatoires de l'article 2 (G14a) | Corps finis : fermé ; L/R′ : fermé sur papier ; D < nd : en partie ; la caractéristique 2 gagne seulement une constante sur de petits cas |
| ✅ Théorème à plusieurs points pour l'article 2 (G14d) | Borne doublée, tout corps ; pas d'échappatoire par placement ; écart log ℓ non fermé (non relu) |
| ✅ Limites de la contraction de tenseurs (G14b) | Annexe D reconstruite ; √log optimal pour les effilages dépendant de la valuation (théorème B, non relu) ; base possible d'un article 3 |
| ✅ Relecture + nouveauté du théorème de l'article 2 (G13a) | Confirmé ; vaut sur tout corps ; semble nouveau |
| ✅ **Lean : théorème de l'article 2 sur tout corps** (G15L) | `AnyField.sum_rank_ge_N2`, forme à deux paramètres |
| ✅ Lean : compléments (G14g) | Article 1 sans l'hypothèse n ≤ M ; schéma trivial ; conséquences finies de l'article 2 |
| ✅ **Pré-version 3 (brouillon)** (G12W) | `resultats/brouillons/pre-v3/` : 41 pages, marqué DRAFT ; section F\*, section Lean à jour, section « Scope of the claims » avec 8 questions ; à faire vérifier par une autre IA (voir `LISEZ-MOI.md`) |
| ✅ Cas général des erreurs par point (G11c) | Lemme Q_{d+2} étendu ; la classe commutative n'est pas optimale (11 contre 12 à (4,1), numérique) |
| ✅ **Lean : théorème principal candidat de l'article 2** (G13b) | `sum_rank_ge_N` : Σ rang(E_k) ≥ N(ℓ, m), prouvé sur ℚ sans hypothèse ajoutée |
| ✅ Relecture du théorème sep(F\*) (G11a) | Confirmé ; la définition de F\* doit être énoncée comme hypothèse ; 6 corrections listées |
| ✅ Relecture de G9c (G10d) | Théorèmes A et B confirmés (relativement à l'annexe D de Claude) ; 2 résumés d'obstructions à corriger |
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
| ❌ Encodages directs pour Q1\* (CSP à un côté, Fourier, gadgets) (G14f) | Coût ≥ 2N(d−1) ; rang ≥ C(s+s′, s) pour toute fonction ayant le support du test de conflit |
| ❌ Intrication sur les plus petits cas (grille 2×3, λ = 1) (G14e) | Faisabilité exactement sur la frontière produit (numérique) ; la question reste ouverte pour λ ≥ 2 |
| ❌ Erreurs par point sur un corps fini (piste 1 de la génération 14) | Le théorème de l'article 2 vaut sur tout corps (G13a) : pas de faille |
| ❌ Géométries de récursion (grilles, boîtes de dimension d) | 2^Ω(m/log m) dans le modèle (G4b, G6b) |
| ❌ Structures non produit pour OpenAI (expanseurs, codes) | Rang ≥ K/2 (G5a) |
| ❌ Combinaisons OpenAI + Claude | Même mur 2^(n/log n) (G4c) |
| ❌ Réfuter Feige (m = O(n)) | Bornes inférieures SOS ; hors de portée |

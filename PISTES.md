# Pistes de recherche

Ce fichier liste toutes les pistes : celles en cours, celles pas encore explorées et celles fermées.
Mise à jour : 24 novembre 2026.

**Rappel :** aucune de ces pistes ne mène directement à P = NP. Ce sont des étapes qui peuvent donner des résultats vérifiables.

Légende : 🔄 en cours · 🆕 pas encore explorée · ⏸️ en attente · ✅ terminée · ❌ fermée (mur prouvé ou mesuré)

---

## Plan des articles (décision de l'auteur, 24 novembre 2026)

**Règle : on attend la validation de l'article 1 (relecture de F. Dias) avant d'écrire en entier les articles suivants.** La recherche, les relectures et Lean peuvent continuer.

| Article | Sujet | État | Portée |
|---|---|---|---|
| **1** | Schémas « produit de clauses » (question du §9 de Dias-Lima) | Version 2 envoyée à F. Dias ; pré-version 3 en brouillon | Borne 2^Ω(n^(2/3)) dans le **modèle de coût produit** ; rien sur la programmation dynamique (sauf : modèle dégénéré), ni sur les autres voies |
| **2** | Erreurs par point dans la voie « interpolation » d'OpenAI | **En attente** (accord de l'auteur requis ; attendre G14a et G14d, puis la validation de l'article 1). Théorème prouvé en Lean sur tout corps | Σ_k rang(E_k) ≥ N(ℓ, m) pour un point par direction ; √log n optimal **dans le modèle de coût de l'exposé d'O'Donnell** ; **État des échappatoires :** corps finis : fermé (Lean) ; multiplicités inégales : fermé (**Lean**, G16L) ; effondrement de rang par L/R′ : fermé sur papier (G14a, non relu) ; degré réduit D < nd : fermé en partie (reste : degré certifié par annulation). Ne couvre pas : la contraction de tenseurs (futur article 3), l'algorithme de Claude |
| **3** (futur) | Limites de la voie « contraction de tenseurs » d'OpenAI (annexe D) | **Idée seulement, pas assez pour un article (G17a, 24/11)** : théorème A confirmé mais très faible ; théorème B en partie (axiome du modèle VT réfuté) ; un lemme propre (précision optimale en boîte noire) ; il manque une borne de rang 2^(M/√u) dans le modèle DF | Visé : le √log est optimal pour la famille des « effilages » dépendant de la valuation ; borne d'aplatissement dans le modèle DF. Ne couvrirait pas : d'autres règles de reconstruction, plusieurs niveaux, l'algorithme de Claude |

Avec les trois articles, les trois voies de la prépublication d'OpenAI seraient étudiées. Aucun ne dit rien sur P vs NP ni sur la validité de la prépublication.

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

## 1. En cours

| Id | Piste | Question précise | Origine |
|---|---|---|---|
| G18a | 🔄 **Précision réelle des noyaux 2-adiques** | L'analyse d'OpenAI paie 40 bits là où 10 suffisent (14 contre 4 sur un paquet creux). Le gain dépasse-t-il une constante ? Effet sur l'exposant. | G17a |
| G18R | 🔄 **Relecture adverse de G18b** | G18b annonce un gain d'un facteur loglog n sur la route d'interpolation (décision / #F mod 3, occurrences bornées), par un évaluateur modulo 3. Esquisse. La caractéristique 3 casse-t-elle un lemme de la prépublication ? Test complet indépendant. | G18b (`code/gen18/g18b/`) |
| G18S | 🔄 **Relecture adverse de G18e** | G18e annonce le théorème A′ : N ≥ 2^(0,399·M/√u) dans le modèle DF, donc 2^Θ(M/√u) avec la borne supérieure de l'annexe D. C'était la condition manquante pour l'article 3. Juste ? Charges non uniformes ? Énoncé prêt pour Lean. | G18e (`code/gen18/g18e/`) |
| G18f | 🔄 **Termes intriqués, λ = 2** (question E1) | Seule piste vers 2^(n^(1−ε)). Calcul à λ = 2 ou question réduite ; E1 est-elle une question de rang asymptotique d'un tenseur explicite ? | G9d, G14e |
| G18g | 🔄 **Problème de Hurwitz et Naslund** | P(N) ≤ N^(1+o(1)) ou P(N) ≥ N^(1+ε) ? Variante exacte utile à l'algorithme ; état de l'art 2026 ; une tentative originale. | N1, G4a, G5b |
| G18h | 🔄 **Pièces hors du cadre « un seul espace gradué »** | Le plancher sur ρ est-il une vraie limite ? Retirer chaque hypothèse du cadre, une par une. | G9a, G9b |
| G18i | 🔄 **Réfutation des 3-CNF aléatoires à densité constante** | Carte temps contre densité avec les algorithmes de 2026 ; lecture en complexité des preuves ; une idée nouvelle testée. | obstacle central |

Dossiers de travail de la génération 18 : `scratchpad/gen18/g18a` à `g18i`.

**État (24/11, 10 h 45 UTC) :** les neuf agents se sont arrêtés à cause de la limite d'utilisation, presque au début. Relance en deux vagues. **En marche :** G18a, G18f, G18g, G18h, G18i, G18R, G18S. **Terminés :** G18b (à relire), G18c, G18d, G18e (à relire).

Note : les dossiers de travail des agents G17b à G17e s'appellent `g18a` à `g18d` (erreur de numéro au lancement ; la génération est bien la 17).

**Consigne de l'utilisateur (24/11) :** dès que les générations 8 et 9 sont terminées, lancer une **génération 10** sur les pistes « pas encore explorées » (section 2) et « en attente » (section 3).

## 2. Pas encore explorées

| Piste | Question précise | Difficulté | Origine |
|---|---|---|---|
| 🆕 **Problème de Naslund / Hurwitz** | Prouver P(N) ≥ N^(1+ε) ; équivaut à un plancher ρ > 1/2 dans l'algorithme de Claude | Très difficile (circuits non commutatifs) | N1, G5b |
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

| 🆕 **Termes « intriqués » : question E1 précise** | Un développement avec 2^o(S) termes et κ_r + κ_c < 1, qui couvre plusieurs rectangles (g ou h ≥ 3). Plus petits cas (2×3, λ = 1) : aucun gain (numérique) ; λ = 1 ne peut pas trancher. Il faut λ ≥ 2 ou une idée structurelle | Très difficile ; calculs lourds | G9d, G14e |

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

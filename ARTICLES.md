# Registre des articles possibles

Ce fichier liste **tout ce qui peut entrer dans un article**, présent ou futur. Il est mis à jour à chaque génération.

**Règle de l'auteur (24 novembre 2026).** On attend la validation de l'article 1 (relecture de F. Dias) avant d'écrire en entier les articles suivants. La recherche, les relectures et Lean continuent.

**Règle de ce registre.** Rien n'est inventé. Un résultat entre dans la partie A ou B seulement s'il existe, avec un fichier dans le dépôt. La partie C liste ce qui n'est **pas** assez pour un article, pour que personne ne le présente comme tel.

**Niveaux de fiabilité** (comme dans `JOURNAL.md`) : **L** = vérifié par Lean ; **P+** = preuve refaite par un relecteur IA de la même famille (pas indépendant) ; **P** = preuve écrite, non relue ; **E** = esquisse ; **N** = calculs seulement. **Aucun résultat n'est relu par un expert humain.**

**Limite commune.** Tous ces résultats parlent des méthodes de deux prépublications d'octobre 2026 (OpenAI, Claude/Anthropic), qui ne sont pas vérifiées. Aucun ne dit rien sur P contre NP.

---

## A. Articles en cours ou prévus

### Article 1 — Schémas « produit de clauses » (réponse au §9 de Dias-Lima)

| | |
|---|---|
| **État** | Version 2 envoyée à F. Dias. Version 3 **en attente** de son retour. |
| **Fichiers** | `resultats/article/` (v2, anglais et français) ; `resultats/brouillons/pre-v3/` (brouillon) ; `lean/PnpLean/Article1/` |
| **Résultat central** | Théorème B : ρ_t + q_t + 1 > min(c_b·sep(F)/√s_t, 2n/3 − s_t). Corollaire : coût 2^Ω(n^(2/3)) dans le modèle de type produit. **Niveau L** pour les deux théorèmes et tous leurs lemmes. |
| **Sur papier seulement** | Lemme sur les 3-CNF aléatoires, corollaire asymptotique, section F\* (P+). |

**Matière prête pour la version 3 :**

| Élément | Niveau | Fichier |
|---|---|---|
| Les deux théorèmes sans l'hypothèse n ≤ M | L | `lean/PnpLean/Article1/Extras.lean` |
| Lemme d'espérance du déficit, anti-concentration | L | `Lemma2.lean`, `Anticonc.lean` |
| Schéma trivial du modèle de type PD | L | `TrivialScheme.lean` |
| Section F\* (η\* = 2,418·10⁻⁴) | P+ | brouillon pré-v3, §6 |
| 7 corrections tirées de la prépublication elle-même | lecture de la source | `code/gen17/g17b/rapport.md` |
| Point de portée : sur le plan réel d'OpenAI le théorème B est vrai par q seul ; il contraint les nouveaux plans à petits groupes | lecture de la source | `code/gen17/g17b/rapport.md` |
| Modèle « rang d'interface » pour les vrais programmes dynamiques (borne √n) | P, non relu | `code/gen14/g14c/notes.md` |

### Article 2 — Erreurs par point dans la route d'interpolation d'OpenAI

| | |
|---|---|
| **État** | Pas de texte. Rédaction en attente (validation de l'article 1). |
| **Fichiers** | `lean/PnpLean/Article2/` ; `code/gen6/g6d/`, `code/gen11/g11b/`, `code/gen14/g14a/`, `code/gen14/g14d/`, `code/gen17/g17c/`, `code/gen17/g17d/`, `code/gen17/g17e/`, `code/gen18/g18b/`, `code/gen18/g18r/` |
| **Résultat central** | Pour J(t) de degré ≤ d avec J(0) = X^{⊗ℓ} et J(τ_k) = K_k + E_k : Σ_k rang E_k ≥ Σ_x max(0, N − d − 2μ(x)). Tout corps, plusieurs points par direction. **Niveau L.** |

**Matière disponible :**

| Élément | Niveau | Fichier |
|---|---|---|
| Théorème à plusieurs points, tout corps | L | `MultiPoint.lean`, `AnyField.lean` |
| Théorème au niveau de la sortie (effondrement par L, R′), degrés non uniformes | L (partie déterministe) ; P (moyenne sur S, lemme glouton 1/(2Δ−1)) | `OutputLevel.lean`, `code/gen17/g17d/` |
| Lemme Q_{d+2} (4 mauvaises arêtes) | L | `QdPlus2.lean` |
| Nouveauté : ni la prépublication, ni O'Donnell, ni Dias-Lima n'ont cette borne | lecture des sources | `code/gen17/g17c/rapport.md` |
| Définition du « cadre X » (6 conditions) où « perte √log optimale » est un énoncé juste | lecture de la source | `code/gen17/g17c/rapport.md` |
| Théorème du temps de calcul (hypothèses H1–H4), exposant Ω(n(loglog n)²/√log n) | P | `code/gen14/g14d/notes.md` |
| Borne inférieure selon la charge B^(Mρ) : Ω(n√(loglog n/log n))·min(ln B, loglog n) | P+ (conditionnel à G14d) | `code/gen18/g18b/`, `code/gen18/g18r/` |
| **Borne supérieure nouvelle :** évaluateur modulo 3, exposant (12,6 + o(1))·n·lnln n/√ln n, pour la décision aléatoire et #F mod 3 sur les 3-CNF à occurrences bornées | P+ (conditionnel aux §3–4 de la prépublication) | `code/gen18/g18b/`, `code/gen18/g18r/` |
| Dualité de quadrature, lemme d'élagage | P | `code/gen17/g17e/notes.md` |

**Ce qui manque.** L'écart O(log ℓ) entre la borne et la construction (ouvert). Le degré certifié par annulation (ouvert). Le théorème du temps de calcul refait avec le coût réel de la prépublication. L'écart √loglog n entre la borne supérieure et la borne inférieure pour une charge constante.

**Portée à écrire.** La conclusion vaut dans un cadre à un niveau sur un corps. La route principale de la prépublication (annexe D) est hors de ce cadre.

### Article 3 — Limites de la route « contraction de tenseurs » d'OpenAI (annexe D)

| | |
|---|---|
| **État** | Pas de texte. Le résultat central existe depuis la génération 18 et il est vérifié par Lean depuis la génération 19. |
| **Fichiers** | `code/gen14/g14b/`, `code/gen17/g17a/`, `code/gen18/g18a/`, `code/gen18/g18d/`, `code/gen18/g18e/`, `code/gen18/g18s/`, `code/gen19/g19a/` ; `lean/PnpLean/Article3/DFRank.lean` |
| **Résultat central** | Théorème A″ : dans le modèle DF, log₂ N ≥ min_i ½·Σ_g |e1_i(g) − e2_i(g)|. **Niveau L** (y compris la forme modulo un idéal propre). Corollaire pour des charges presque uniformes : N ≥ 2^((0,564 − o(1))·M/√u), **niveau P+** (identité de comptage sur papier). |

**Matière disponible :**

| Élément | Niveau | Fichier |
|---|---|---|
| Reconstruction fidèle de l'annexe D ; modèle DF | P+ | `code/gen14/g14b/`, `code/gen17/g17a/` |
| Théorème A (ancien, 2^(2M/4^u)) | P+ | `code/gen17/g17a/rapport.md` |
| Théorèmes A′ et A″ (trois formes, et la forme modulo un idéal propre) | **L** | `lean/PnpLean/Article3/DFRank.lean` |
| Corollaire : constante 1/√π sur des charges uniformes ; lemme de restriction | P+ | `code/gen18/g18s/notes.md` |
| Famille précise où « √log est optimal » est un théorème (côté inférieur) | P+ | `code/gen18/g18s/notes.md` |
| Contre-exemple : avec des charges non uniformes, N = 2^(M/u) est atteint | P+ (petits cas exacts) | `code/gen18/g18s/notes.md` |
| Précision Θ(M·D_t/u) optimale pour la reconstruction boîte noire (constante ≈ 12) | P+ | `code/gen17/g17a/rapport.md` |
| La structure des noyaux n'ajoute rien à la précision (rang plein ⇒ même besoin) | N (17 cas) | `code/gen18/g18a/notes.md` |
| Récurrence de l'annexe D : la racine domine ; la profondeur achète √loglog | P | `code/gen18/g18d/notes.md` |
| Proposition B′ (borne au niveau de la conception), lemme T | P | `code/gen18/g18d/notes.md` |
| Invariants de Smith des blocs de Cauchy : 2^((a−i+1)+v₂(ζ)+2k+2v₂(k!)) | **P+ (prouvé par G26-P6, relu par G26-P7)** ; la conjecture C1 est fausse telle qu'écrite, remplacée par C1′ | `code/gen18/g18d/notes.md`, `code/gen26/p6/notes.md` |
| Théorème U (G26-P8, relu par G26-P12 : P+, hypothèse v₂ ζ ≥ −2 explicite) : Smith(M(E)) = union disjointe sur k de Smith(Per_k(E)) (Γ unitriangulaire pour l'inclusion, det Γ = 1). Théorème A2 (C1′ côté droit au niveau a − 2) : P+. **(A) est RÉFUTÉ** (G26-P12 : a = 6, i = 1, k = 2, invariant 115 > 114) ; C1′ (côté droit) est faux pour des sous-blocs de Cauchy (6 sur 156), donc une preuve doit utiliser tout l'ensemble d'indices du bloc de croisement ; pour les blocs complets à a ≥ 5 : ouvert. Prouvé : (A) au niveau a − 2 pour tout k ; forme faible de (A) pour k = 2 à tous les niveaux. L'optimalité de l'énumération des modes à √2 près reste une conjecture | P+ (U, A2) ; P (le reste) | `code/gen26/p8/notes.md`, `code/gen26/p12/notes.md` |
| Gain constant possible sur l'annexe D (modes 4 fois redondants, effilage 1/√3) | E, non vérifié | `code/gen18/g18d/notes.md` |
| Les contraintes de charge viennent de la conception (lemme D.3, (D.6)), pas de l'instance | lecture de la source | `code/gen19/g19a/notes.md` |
| Inégalité de séparateur : chaque direction coupe ≥ φ·M·min(1, 2(1−β)) fils (le coût 2^(φ′M) vaut pour l'identité de conditionnement) | P+ (G20f) | `code/gen19/g19a/notes.md` |
| Théorème E : log₂ N ≥ min_i M_i/(2u) − u pour tout étiquetage (preuve complète ; forme plus forte par direction) | P+ (G20f) | `code/gen19/g19a/notes.md` |
| Limite de la méthode : sur l'étiquetage circulaire, A″ = M/u exactement ; l'aplatissement ne donne pas M/√u en général | P+ (G20f) | `code/gen19/g19a/notes.md` |

**Encadrement (G25-P10) :** à un appel, l'annexe D coûte entre 0,564·M/√u (théorème A″) et 12,8·M/√u avec les constantes rigoureuses de G25-P3 : écart d'un facteur constant 22,7 seulement (79 pour le papier tel qu'écrit).

**Ce qui manque.** Une borne 2^Ω(M/√u) pour les étiquetages non uniformes : G19a prouve seulement M/(2u) − u, et montre que l'aplatissement ne peut pas faire mieux en général (problème ouvert du « paquet en escalier »). Les schémas à plusieurs niveaux et la re-sommation non énumérative. La borne supérieure de l'annexe D n'est pas vérifiée.

---

## B. Candidats à un autre article (matière réelle, mais incomplète)

### B1 — Note « plus long chemin et plus long cycle en 2^(o(n)) sur les graphes creux »

| | |
|---|---|
| **Résultat** | Conséquence des prépublications : algorithmes sous-exponentiels pour le plus long chemin et le plus long cycle sur les graphes creux. |
| **Niveau** | P+ (preuve confirmée par G10c, environ 570 000 vérifications de l'identité de comptage). |
| **Fichiers** | `code/gen9/g9e/consequences.md`, `code/gen10/g10c/` |
| **Limites** | Conditionnel aux deux prépublications. Résultat modeste. La conséquence pour l'arbre de Steiner est déjà revendiquée par Wu et Wang (arXiv:2611.02548), à citer. |
| **État** | En attente : décision de l'auteur. |

### B2 — Compromis de rang pour l'opérateur d'échange (théorème C)

| | |
|---|---|
| **Résultat** | Théorème C : P est l'échange E⊗G⊗W → G⊗E⊗W avec un spectateur (dimensions m, n, Q). Si P = A + B avec rang A ≤ r et rang de Schmidt d'opérateur de B ≤ b, alors pour tout entier k ≤ m : (a) r ≥ k·Q·(n − k·b) ; (b) b·Q ≥ k·(n·Q − k·r) ; (c) b < m·n ⇒ r ≥ Q. Corollaire : r·b ≥ (2/9)·Q·M². Tout corps. |
| **Niveau** | **L** depuis G20g : `lean/PnpLean/Outils/SwapRank.lean` (parties (a), (b), (c), formes miroirs, forme produit). Avant : deux preuves indépendantes (G18g, G18f), unifiées par G19e ; environ 15 000 décompositions exactes sans violation. |
| **Fichiers** | `code/gen19/g19e/notes.md` (énoncé commun, preuve courte, errata) ; `code/gen18/g18g/notes.md`, `code/gen18/g18f/notes.md` |
| **Compléments (G19b, relus par G20f : niveau P+)** | Théorème 3 : q = 1 ⇒ p ≥ N(N−1)/2 (optimal). Théorème 4 : q = 2 ⇒ p ≥ N, donc P(3) = 4. Théorème 5 (corrigé : il faut 1 ≤ k ≤ N) : p·q ≥ (1 − o(1))·(N²/4)·q/(q − ½). Les quatre décompositions nouvelles sont des identités sur ℤ. Table exacte des paires réalisables pour N ≤ 3, partielle pour N = 4 à 6, avec des décompositions à coefficients dans {0, ±1} vérifiées sur ℚ. Valables pour les décompositions exactes en caractéristique ≠ 2, pas au sens « bord ». `code/gen19/g19b/notes.md` |
| **Barrière associée (G20d, relue par G21e : niveau P+)** | Théorème B : pour Φ = somme par blocs de A_t·Fl_{S_t}(X)·B_t sur un corps de caractéristique 0, rang Φ(X)² ≤ N²·c1·c2 ; donc aucune telle application n'exclut une paire avec p·q ≥ N²/4, pour aucun opérateur, exact ou au bord. Ne couvre pas : les compositions avec des applications non locales, les Φ linéaires générales, la caractéristique p. Le théorème C (a), (b) est dans cette classe ; en exposants, ces méthodes s'arrêtent à la droite. Frontière générique calculée pour N ≤ 4 : l'opérateur d'échange est spécial, pas générique. `code/gen20/g20d/notes.md` |
| **Limite des méthodes algébriques de bas degré (G21d, calculs exacts)** | À N = 3 : dim J(2,2) = 64 ; aucun des 168 aplatissements de Koszul/Young, et aucun polynôme de degré ≤ 6, ne s'annule sur J(2,2) ; aucune équation séparante en degré 7. L'échange est à distance² = 2 de J(2,2) (numérique). Complète la barrière : même hors de la classe B, il n'y a pas de témoin de bas degré au plus petit cas. `code/gen21/g21d/notes.md` |
| **Invariant de distance (G22d, G23d)** | Masse captée c = N² − d² entière à tous les points calculés ; loi q = 1 : d² = N² − m² − 2j (numérique) ; prouvé (papier) : c ne décroît pas avec N, est sur-additive, symétrique en (p, q), et **c(N;1,1) ≤ 8,081, c(N;p,1) ≤ 11,7·p + O(√p) pour tout N** — premières bornes analytiques indépendantes de N, alors que les méthodes de rang perdent Θ(N(p+q)). Pour q ≥ 2 : borne c ≤ 13p + 4q prouvée seulement pour une sous-classe (G24c). `code/gen22/g22d/notes.md`, `code/gen24/g24c/notes.md` |
| **Modèle monomial (G22b-bis, relu par G23b ; modèle ST vérifié par Lean, G23e : niveau L)** | Théorème de correspondance exact : décompositions monomiales symétriques ⟺ matrices d'intercalates de Hrubeš-Wigderson-Yehudayoff à graphe d'intercalation biparti (formules entières de sommes de carrés), avec p et q lus sur les couleurs. **Dans ce modèle, borne au-delà de la droite : p²·q ≳ N⁴/8, point équilibré ≥ N^(4/3)/2** (lemme clé : amélioration de l'étape k^(3/2) de HWY en k²/2). Vérifiée exactement contre toutes les frontières SAT pour N ≤ 5. Ce n'est pas une borne pour le problème général (l'optimum à N = 5 sort du modèle). `code/gen22/g22b/notes.md` Valeurs exactes du modèle ST (G24d, SAT avec preuves DRAT) : 1, 3, 4, 7, 7, 8, 8 pour N = 2 à 8 ; encadrements jusqu'à N = 16. |
| **Extrémités de la courbe (G24g ; vérifié par Lean, G25-P6 : niveau L)** | Théorème 1 : f(N, 2) = f(N, 1) = N(N−1)/2 pour tout N, sur tout corps (deux produits ne font pas mieux qu'un seul) ; lemme de triangularisation T ; borne N²/(2(r−1)) pour r ≥ 2. La frontière est exactement plate sur deux pas près de chaque extrémité. Lean : `lean/PnpLean/Outils/Triang.lean`. `code/gen24/g24g/notes.md` |
| **Modèle monomial-C (G26-B5, G26-P3 ; relu par G26-P9 : niveau P+, avec réserves : matrices signées ±1 seulement, sinon Q = min Σ‖G_l‖·‖G′_l‖ ; caractéristique 0 seulement ; bornes vides pour N ≲ 20)** | Lemme clé : pour v, w dans ker(1 − Θ), ⟨v, Θ_s w⟩ = ⟨v, w⟩ (Θ_s partie symétrique). Théorème A : p ≥ (1/2 − ε)·N² − C(q, ε)·N ; donc p ≥ (1/2 − o(1))·N² si q = o(log N) ; certificats rationnels exacts (q ≤ 3 : 0,4987·N² − 720·N ; q ≤ 10 : N²/3 − 62·N). Il suffit que Σ‖G_l‖·‖G′_l‖ ≤ q : une construction générale à O(1) termes avec p ≤ (1/2 − c)·N² demande des facteurs de norme non bornée. Théorème B : p ≥ c·N²·ln N/q (le lemme MC2 donnait ln q). Théorème G : si tous les facteurs commutent avec une représentation d'un groupe fini, p ≥ (N² − Σ d_ρ·m_ρ²)/2 pour tout q. Théorème C (conditionnel) : une borne sur l'énergie des mots de longueur 4 donnerait un gain d'exposant. Obstruction (esquisse) : la méthode des polynômes ne peut pas dépasser N²·log N/q. Les 5 optimums exacts connus vérifient p ≥ N²/2 − q·N/3 (conjecture « lemme fort MC »). `code/gen26/p3/notes.md` |
| **Calculs exacts dans le modèle MC (G26-P4 ; lemmes et relaxation GF(2) relus par G26-P16 : P+ ; branches UNSAT r = 4, 5 sur kissat seul)** | Relaxation SAT sur GF(2) (rang_ℚ ≥ rang_GF2, donc UNSAT sur GF(2) ⇒ UNSAT dans MC). Lemmes implicites prouvés (tout corps) : Σ|E_l|·|C_l| ≥ N² − p ; rang ≥ N² − |dom ∩ im| ; lemme en escalier ; corollaire : si tous les rangs sont ≤ r (r ≥ 2), p ≥ N² − q·r² + 1. Résultats : N = 4 complet (p minimal = 6, 6, 5, 4, 3, 1, 1, 1, 1 pour q = 1 à 9 ; côtés UNSAT avec preuves DRAT vérifiées) ; **(9,3)@5 UNSAT** (dans MC, p minimal à q = 3, N = 5 est 10) ; **(6,4)@5 UNSAT** (p minimal à q = 4 est 7) ; branches r = 4, 5 sur kissat seul (sans DRAT) ; (5,5)@5 et (11,3)@6 : non décidés. `code/gen26/p4/` |
| **Profil spectral et nouveaux points exacts (G26-P11 ; les 9 points sont vérifiés par G26-P15 avec un vérificateur indépendant en arithmétique exacte : niveau P+ pour les points, P pour le reste)** | Spectre de Θ aux optimums : {1^(n−p), −1, −2} (plus −3 pour Q8) ; ‖Θ_s‖ = 2 ; Θ diagonal par blocs sur les N² cases (cases fixes + blocs-cliques de taille 2 à 4, chacun de coût 1 : p = nombre de blocs). Neuf nouveaux points exacts (vérifiés modulo 2, 3, 5, 7 et revérifiés) : (8,6)@6 ; (20,3), (17,4), (14,5), (12,6), (11,7)@7 (avec (8,8), la frontière N = 7 est complète) ; (22,5)@8, (20,6)@8 ; (22,6)@9. Q8 viole la règle empirique p ≥ N²/2 − qN/3 (constante ≥ 3/8). Dans le modèle MC à blocs-cliques : T ≥ (5 − 2θ)/3 (au-dessus de la droite pour θ > 1, mais dans le cercle). Aucun point au-delà de la droite. `code/gen26/p11/` |
| **Conséquences** | Aucun gain par termes intriqués sur un seul croisement. La courbe des paires réalisables de l'algorithme de Claude est au-dessus de la droite T = 2 − θ. P(3) = 4 sur F₂. |
| **Antériorité (G19e, puis G22a)** | L'objet est le problème 12 de Naslund (arXiv:1701.04475) = question 1.2 de Karam (arXiv:2310.11356) = largeur w(ID_N) de Hrubeš-Wigderson-Yehudayoff = « PT-rank » de Rossman-Zhu (arXiv:2512.01227). **G22a a lu neuf articles (trois en entier) : le théorème C n'apparaît dans aucun ; seul le cas k = 1 (p + q ≥ N) est connu ; aucun article ne traite p et q séparément, donc la question du compromis elle-même est nouvelle.** La preuve reste simple. Le théorème B (barrière) n'est pas trouvé ; il doit citer Efremenko-Garg-Oliveira-Wigderson (arXiv:1710.09502). Non lus : Siemaszko, Debnath (non accessibles). |
| **Avis honnête** | Le théorème C seul est un lemme. **Avec la barrière (G20d, relue), l'absence de témoin de bas degré (G21d), l'invariant de distance (G22d) et surtout la borne au-delà de la droite dans le modèle monomial (G22b-bis), l'ensemble forme maintenant une note cohérente sur le compromis (p, q) du problème 12 de Naslund : une question nouvelle (G22a), une borne générale vérifiée par Lean, la barrière qui lui correspond, et une borne plus forte dans un modèle naturel.** Condition restante : relecture humaine. (Lean : théorème C et borne du modèle monomial ST, faits.) Limite : rien de tout cela ne donne ρ > 1/2 dans le problème général. |
| **Ce qui manque** | L'antériorité est maintenant contrôlée sur la littérature principale (G22a). Reste : Siemaszko et Debnath (non accessibles), et une relecture humaine. |

### B3 — Réfutations algébriques sous-exponentielles et carte de la réfutation des 3-CNF aléatoires

| | |
|---|---|
| **Résultat** | (a) Toute k-CNF insatisfaisable à m clauses a une réfutation IPS sur ℚ de taille 2^(O(km/√log m)). (b) Carte temps contre densité : les algorithmes de 2026 sont les meilleurs seulement pour Δ ≤ c·log n, avec 2^(O(n/√(Δ·log n))) (proposition A). |
| **Niveau** | (a) E : une étape non contrôlée (la récursion est-elle un circuit arithmétique aveugle aux tables ?). (b) P pour la proposition A ; les références sont lues dans des résumés. |
| **Fichiers** | `code/gen18/g18i/notes.md` |
| **Limites** | Conditionnel à la prépublication de Claude. Nouveauté établie par une seule recherche. Rien pour Frege ni Extended Frege. |
| **Ce qui manque** | Contrôler l'étape (a). Une vraie recherche d'antériorité. Sans cela : une remarque, pas un article. |

### B5 — « Battre la largeur arborescente » après 2026, et ce que la conjecture du rang asymptotique ne donne pas

| | |
|---|---|
| **Résultat** | (T1) Conditionnel à la prépublication de Claude : #CSP binaire pondéré sur un graphe de degré max Δ à v variables et domaine s en temps s^(O(Δ·v/√log v)) ; sur les expanseurs de degré borné : s^(O(tw/√log tw)). (T2, T3) La capacité de liaison est O(k·log Δ/log k) ; le facteur log tw laissé ouvert par le théorème de Marx est une quantité de routage, et le fermer ne donnerait rien pour 3-SAT. (Théorème R) Tout développement du tenseur de croisement en termes de rang un sur les colonnes a au moins 2^(S/2) termes : la conjecture du rang asymptotique et ω = 2 changent seulement une constante. (Circuits) 3-SAT en 2^(O(V/a(V))) ⇒ E^NP sans circuits de taille c·n·a(n). |
| **Niveau** | T1 : E (esquisse, conditionnelle, constantes non contrôlées). T2, T3, R, circuits : P (non relus). |
| **Fichiers** | `code/gen21/g21c/notes.md` |
| **Limites** | Nouveauté non établie : T2 est sans doute du folklore ; l'article de Marx n'a pas été lu (hôte bloqué) ; plusieurs références citées de mémoire ou lues dans les résumés seulement. T3 couvre une seule famille de réductions. |
| **Avis honnête** | Matière réelle mais mince et non relue. Aujourd'hui : des remarques à mettre dans la section « travaux reliés » d'un autre article. Un article seul demanderait une relecture, la lecture des sources, et un contrôle de T1. |

### B6 — Version explicite et réoptimisée des deux algorithmes de 2026 (« rang inférieur »)

| | |
|---|---|
| **Résultat** | (Claude, G25-P2 puis G25-P9) Théorème A′ : #3-SAT en 2^((c1 + o(1))·W/√log₂ m) avec c1 = 10,6 (déterministe) ou 7,47 (Las Vegas), contre 42,4 dans le papier ; avec le lemme RS et des niveaux déséquilibrés (G25-P9) : 6,81 / 4,82 ; **puis (G26-B7, relu par G26-P2 : P+) c1 = 4,05 à 4,12 déterministe et Las Vegas (regroupement stratifié, lemme d'ordre aléatoire, dérandomisation), et 6,16 en mémoire polynomiale ; compromis temps-mémoire du balayage (1 − M)/3 ; seuil contre la force brute à densité 4,26 : log₂ m* = 2 237 (papier : 276 000)** ; **et une explication de structure : la loi « de Pythagore » (Prop. P) montre qu'aucun arbre de niveaux ne change la puissance du logarithme** ; lemme de balayage nouveau (base e0 ≤ 1/4 + 1/(2K0)) ; seuils explicites contre la force brute, le meilleur #3-SAT et PPSZ. (OpenAI, G25-P3 puis G26-P6) annexe D avec constantes explicites : E ≤ 75·k·m/√log₂ s (rigoureux si l'annexe D est juste), **puis 37,5·k·m/√log₂ s grâce à la formule de Smith des blocs de Cauchy, prouvée par G26-P6 et relue par G26-P7 (P+)**, ε légal 1/960 au lieu de 2⁻¹⁷, seuil 780 fois plus bas ; corrections de formulation de (D.5). **Conséquences explicites (G25-P10, sous une hypothèse non vérifiée) :** Max-Cut, spectre de Max-3-SAT, #CSP binaire, plus long cycle / chemin à degré borné (constante 5,3). |
| **Niveau** | P+ pour le théorème A′ de Claude (relu par G25-P8, code indépendant) ; P pour G25-P3 et G25-P9 ; conditionnel aux deux prépublications. Réserve : avec le balayage, la mémoire devient exponentielle (≈ 2^(0,95 × temps)). |
| **Fichiers** | `code/gen25/p2/`, `code/gen25/p3/`, `code/gen25/b4/` |
| **Troisième tour (G26-P5 ; relu par G26-P10 : R2F, R3″, c1 = 3,82 et seuils au niveau P+ ; S2 corrigé)** | Lemme de regroupement de Freedman (R2F) : ζ ≈ (2 + 0,85·√(L·L2))/t au lieu de (4 + √(8·L·L2))/t ; gain gratuit au niveau de base (R3″ : Hoeffding-Azuma sur intervalle, facteur de Serfling). c1 : 4,05 → 3,92 → **3,82** (déterministe = Las Vegas). Seuils à densité 4,26 : **1 943** contre 2^n (avant 2 237), 3 821 contre 1,6423^n, 13 410 contre 1,30704^n ; densité 1 : 59,3 / 93,0 / 355. Mémoire polynomiale : 6,16 → 5,88. Plancher de cette conception : regroupement aléatoire stratifié ≈ 3,73 (esquisse) ; tout lemme de regroupement uniforme ≈ 3,36 (comptage, asymptotique) ; toute conception de regroupement : α ≥ 1. Changement de structure S2 (cliques filles pondérées par arête ; esquisse) : P5 donnait 3,23, mais G26-P10 trouve une erreur d'un facteur 2 dans le coût de remplissage (les deux diagonales sont remplies) : **c1 ≈ 3,31** corrigé, seuils 1 514 / 2 911 / 10 020, mémoire polynomiale 4,85 ; seulement 0,05 sous le plancher 3,36. **G26-P13 (relu par G26-P14 : E, C, D confirmés ; B confirmé avec correction) :** S2 écrit en entier (le second ordre est réglé par un placement aléatoire : les plafonds se concentrent ; c1 = 3,312 entièrement compté) ; nouveau « découper au lieu de remplir » et borne de covariance (iii) prouvée : **c1 = 3,166**, seuils 1 352 / 2 683 / 9 269 à densité 4,26 (termes d'ordre inférieur négligés comme dans P2/P5/P10). **Correction de P14 :** le terme de concentration ε_pl = C·√(ln g/h) n'est pas négligeable aux bases (g, h ≈ 2^8 à 2^12) ; avec C = 1 à 3,7 : **c1 ≈ 3,18 à 3,21**, seuils ≈ 1 368 à 1 393 contre 2^n (P+ pour cette fourchette ; la constante C exacte n'est pas fixée). `code/gen26/p5/notes.md` |
| **Avis honnête** | Contribution modeste mais propre et vérifiable : elle corrige et précise les constantes, et montre que les algorithmes restent sans effet pratique (seuils astronomiques). Utile comme section d'un article ou comme note technique ; les auteurs des prépublications seraient les premiers intéressés. |

### B7 — Critère du « coin » : l'étape 2 dépend de la pente de la frontière près de θ = 2, pas du point équilibré

| | |
|---|---|
| **Résultat** | Dans le modèle de conception M de l'algorithme de Claude, avec des niveaux g × h déséquilibrés (lemme U) : α* = inf{α : h_F(μ) ≥ ‖(μ, 1 − μ)‖_{1/(1−α)} pour tout μ}, h_F = fonction support de la frontière de croisement. Donc α* ≥ max(log₂(1/ρ), γ) et **α* = 1 ⇔ γ = 1** (γ = exposant de coin de l'enveloppe concave). Une pente finie c au coin donne e(x) ≤ O(c·log x/x), soit un temps 2^(O(c·m·log log m/log m)). Exemple : la frontière p²q = N⁴ (ρ équilibré = 2/3) donne α* = 1. La loi de Pythagore de G25-P9 est le cas du cercle. |
| **Niveau** | P+ (G27-B2, relu par G27-R1 avec trois corrections : étape de preuve ajoutée, équivalence corrigée — « pente finie ⇒ étape 2 » seulement dans ce sens —, vitesse plus forte que B2). Niveau conception seulement ; dépend du lemme U (esquisse de G25-P9, confirmé par R1) et de la prépublication de Claude. |
| **Fichiers** | `code/gen27/b2/notes.md`, `code/gen27/r1/notes.md` |
| **Limites** | Aucune famille de pièces connue n'a de pente finie au coin (l'annexe B donne δ = η²/4, γ = 1/2 ; G27-T3, G28-C1 : toutes les familles testées ont un coin quadratique). Côté OpenAI : pas d'analogue (G28-O1). Si un coin de pente c ≤ 4 existait, il compterait dès log₂ m ≈ 100 à 900 (G28-O1, non relu). Ce n'est pas une borne inférieure pour les algorithmes. Gains seulement à des tailles astronomiques. |
| **Avis honnête** | Changement de cible important pour le projet : le problème ouvert devient « rigidité de rang de Kronecker du swap » au coin, et non plus le problème 12 de Naslund au point équilibré. Aujourd'hui : une section d'un article ou une note courte, seulement si une relecture humaine confirme le lemme U. |

### B4 — Bibliothèque Lean d'outils génériques

| | |
|---|---|
| **Résultat** | 23 énoncés génériques vérifiés et absents de Mathlib (dont le théorème C et le théorème A″) : la dimension du noyau divise le déterminant d'une matrice de polynômes ; anti-concentration avec constante explicite (1 − e⁻²)/4 ; rang de A(0) ≤ Σ rang E_k. |
| **Niveau** | L. |
| **Fichiers** | `lean/PnpLean/Outils.lean`, `lean/OUTILS.md` |
| **Nature** | Ce n'est pas un article de recherche. C'est une contribution possible à Mathlib, ou une courte note de formalisation. |
| **Ce qui manque** | Généraliser de ℚ à tout corps ; remplacer `uavg` par `Finset.expect` ; **une licence** (décision de l'auteur). |

---

## C. Ce qui n'est PAS assez pour un article (à ne pas présenter comme tel)

| Résultat | Pourquoi ce n'est pas un article | Fichier |
|---|---|---|
| Dictionnaire commun des deux algorithmes de 2026 ; loi « gain = 1/√D » ; six hybrides sans gain | Loi heuristique, pas un théorème. Calculs sur de très petits cas. Utile comme section d'exposé dans l'article 3, pas seul. | `code/gen18/g18c/notes.md` |
| Diagnostic du plancher 1/√2 (hypothèse (V), profil de croissance) | Aucun objet construit, aucun plancher nouveau prouvé. Le profil en escalier tel que défini est exclu (G19b). C'est une carte pour la recherche. | `code/gen18/g18h/notes.md`, `code/gen19/g19b/notes.md` |
| Reformulation du problème de Hurwitz-Naslund (droite contre cercle) | Aucune borne améliorée. La partie prouvée est dans B2. | `code/gen18/g18g/notes.md` |
| Sept familles de constructions sous le cercle (Pauli, Clifford, fenêtre d'Alman-Williams, Weyl, Hrubeš) | Résultat négatif : aucune ne passe sous le cercle. Le lien exact avec la non-rigidité (H = G + H·D·H) est propre, mais c'est une remarque. La borne « ρ ≥ 2/3 » est vraie seulement sous l'hypothèse de concavité (G21e) ; sans elle, les contraintes de dimension permettent ρ → 1/2. | `code/gen20/g20c/notes.md` |
| Pièces adaptées à l'instance ; mesures d'expansion des 3-CNF aléatoires | Résultat négatif attendu. La proposition R (le croisement restreint reste un échange complet) est une remarque utile pour la portée de l'article 3 ou de la note B2. | `code/gen20/g20e/notes.md` |
| Valeurs exactes de l'annexe B de Claude à N fini (triviale pour N ≤ 256 ; exposant < 1,6 seulement vers N ≈ 6·10¹¹) | Fait exact et utile comme remarque sur la portée pratique de l'algorithme de Claude ; pas un résultat seul. Peut entrer dans une section de la note B2. | `code/gen20/g20b/notes.md` |
| Théorème abstrait de l'annexe B pour les modules, et neuf familles de modules (Fock, tordus, Neveu-Schwarz, Ramond, symplectiques) | Résultat négatif : tous au-dessus du cercle. Le théorème abstrait (axiomes A1 à A4) et la pénalité d'espace fondamental sont propres, mais reposent sur des faits standards non refaits. Le contrôle exact « Ramond n = 8 = identité des huit carrés » est une remarque. | `code/gen21/g21a/notes.md` |
| Caractéristique p : la récursion de Claude sur F_{p^k} pour la décision ; annexe B p-restreinte ; Frobenius ; Dvir-Edelman | Résultat négatif : même cercle. La réparation de l'annexe B en caractéristique p (base duale modulo p, d = 1) est propre et vérifiée sur quatre cas, mais elle ne donne aucun gain : une remarque. | `code/gen21/g21b/notes.md` |
| Familles structurées pour le compromis (Cayley-Dickson jusqu'à N = 256, Clifford, groupes, Paley, tournois) | Résultat négatif : exposant 2. La mesure « le défaut des sédénions a un rang Θ(N²) » est un fait exact, une remarque. | `code/gen22/g22c/notes.md` |
| Barrière de Cardy (G22e) et règle de sélection de Voronoï (G23c) : √2 optimal pour toute algèbre de sommets fixe, y compris avec regraduation ; théorème abstrait corrigé ; carte des échappatoires | Résultat négatif partiel, propre, non relu. Le théorème A (élémentaire) et la remarque « pas de barrière combinatoire au-delà de la droite » pourraient former une section d'un article sur la construction de Claude, après relecture. Pas un article seul. | `code/gen22/g22e/notes.md` |
| Ansatz de groupe dans le modèle monomial (G23a) : forme close du point Q8 ; obstruction d'anisotropie (Chevalley-Warning) | Résultat négatif propre, non relu. La forme close du point Q8 (lien avec SU(3,2) et F₄) et l'obstruction pourraient compléter la note B2. Pas un article seul. | `code/gen23/g23a/notes.md` |
| Pièces de Hrubeš lues de façon non symétrique (G24a) : courbe exacte, limites du cadre de Clifford ; petite erreur de compte dans Hrubeš (704 contre 768 à N = 64) | Résultat négatif propre ; la courbe exacte et les limites (L1), (L2) pourraient être une section de la note B2 (« les constructions connues vues sur le compromis »). Pas un article seul. | `code/gen24/g24a/notes.md` |
| Question E1 (termes intriqués) | Non tranchée en général. La partie prouvée est dans B2. Le plus petit cas ouvert (N = 3) ne montre aucun gain, mais c'est un résultat numérique, pas une preuve. | `code/gen18/g18f/notes.md`, `code/gen19/g19c/notes.md` |
| Corollaire sur le problème de Hurwitz (borne N^√2) | **Déjà publié** par d'autres (Siemaszko, Debnath). | `JOURNAL.md`, résultat G4a |
| Conséquence pour l'arbre de Steiner | **Déjà revendiquée** par Wu et Wang (arXiv:2611.02548). | `code/gen9/g9e/` |
| Pistes fermées (algorithmes locaux, LP, SOS, anneaux Z_m, certificats FKO, etc.) | Résultats négatifs attendus ou connus. | `PISTES.md`, section 5 |

---

## D. Suivi

| Date | Changement |
|---|---|
| 24/11/2026 | Création du registre après la génération 18. Articles 1, 2, 3 ; candidats B1 à B4. |
| 25/11/2026 | G28-R1 : G27-T3 relu (P+) ; point (106,106)@32 ; pavage de Hurwitz-Radon (s, s) avec s = N·⌈N/ρ(N)⌉ : meilleur que N²/8 pour N ≥ 64. |
| 25/11/2026 | G27-T3 : problème du coin, aucune famille à coin linéaire ; théorème S (esquisse). |
| 25/11/2026 | G27-T2 : cinq nouveaux points exacts (B2), dont (25,5)@10. G27-T1 : annexe B, k borné. |
| 25/11/2026 | G27-B2 + G27-R1 : nouveau candidat B7 (critère du coin), P+. |
| 25/11/2026 | G27-B1 : point numérique (62,12)@16 (non exact). |
| 25/11/2026 | G26-P16 : lemmes de P4 et relaxation GF(2) confirmés (P+). |
| 25/11/2026 | G26-P4 : calculs SAT exacts dans le modèle MC ; (9,3)@5 et (6,4)@5 impossibles ; non relu. |
| 25/11/2026 | G26-P14 : P13 relu ; c1 ≈ 3,18 à 3,21 avec le terme de concentration gardé (P+) ; seuil ≈ 1 370 à 1 400. |
| 25/11/2026 | G26-P15 : les 9 points exacts de P11 sont confirmés (rang exact sur ℚ, deux méthodes) et non dominés. |
| 25/11/2026 | G26-P11 : profil spectral des optimums ; 9 nouveaux points exacts (B2) ; aucun point au-delà de la droite. |
| 25/11/2026 | G26-P13 : S2 rigoureux (3,312) ; découpage + covariance : c1 = 3,166, seuil 1 352 ; non relu. |
| 25/11/2026 | G26-P12 : théorèmes U et A2 relus (P+) ; (A) réfuté ; C1′ reste ouvert pour les blocs complets. |
| 25/11/2026 | G26-P1 : enveloppe de l'annexe B, k borné par 3 dans l'ansatz structuré (résultat négatif, pas de nouvel élément d'article). |
| 25/11/2026 | G26-P8 : seuils d'OpenAI (ligne H) recalculés : 2,25e5 ; théorèmes U et A2 vers C1′ (article 3), non relus. |
| 25/11/2026 | G26-P10 : P5 relu ; c1 = 3,82 (P+), seuil 1 943 ; S2 corrigé : 3,31, seuil 1 514. |
| 25/11/2026 | G26-P9 : bornes du modèle monomial-C (G26-P3) relues (P+), avec réserves. |
| 25/11/2026 | G26-P5 : constante de Claude 3,82 (Freedman) et 3,23 (structure S2, en partie esquissée) ; seuil 1 943 puis 1 419 ; non relu. |
| 25/11/2026 | G26-P3 : bornes inférieures dans le modèle monomial-C (B2), gain logarithmique seulement ; non relu. |
| 25/11/2026 | G26-P7 : formule de Smith relue (P+). |
| 25/11/2026 | G26-P6 : formule de Smith prouvée (article 3 et B6). |
| 25/11/2026 | G26-P2 : B7 relu (P+) ; seuil 2 237. |
| 25/11/2026 | G26-B7 : c1 = 4,12 et version à mémoire polynomiale ajoutées à B6 (non relues). |
| 25/11/2026 | G25-P10 : encadrement ajouté à l'article 3 ; conséquences explicites ajoutées à B6 ; génération 25 terminée. |
| 25/11/2026 | G25-P8 : théorème A′ relu (P+). |
| 25/11/2026 | G25-P9 : constante 4,82 / 6,81 et loi de Pythagore ajoutées au candidat B6. |
| 25/11/2026 | G25-P2 : nouveau candidat B6 (versions explicites des algorithmes de 2026). |
| 24/11/2026 | G25-P6 : f(N,2) = N(N−1)/2 vérifié par Lean ; bibliothèque B4 : 23 énoncés. |
| 24/11/2026 | G24g : théorème f(N,2) = N(N−1)/2 ajouté au candidat B2. |
| 24/11/2026 | G24d : valeurs exactes du modèle ST avec preuves DRAT, ajoutées au candidat B2. |
| 24/11/2026 | G24c : distance pour q ≥ 2 (sous-classe) ajoutée au candidat B2. |
| 24/11/2026 | G24a : pièces de Hrubeš (négatif), partie C. |
| 24/11/2026 | G23a : résultat négatif (ansatz de groupe), partie C. Génération 23 terminée. |
| 24/11/2026 | G23d : bornes analytiques indépendantes de N pour q = 1 ajoutées au candidat B2. |
| 24/11/2026 | G23e : la borne du modèle monomial (ST) est vérifiée par Lean ; bibliothèque B4 : 20 énoncés. |
| 24/11/2026 | G23c : règle de sélection de Voronoï ajoutée à la ligne « barrière de Cardy » (partie C). |
| 24/11/2026 | G23b : la borne du modèle monomial passe au niveau P+. |
| 24/11/2026 | G22b-bis : borne au-delà de la droite dans le modèle monomial ; le candidat B2 devient une note cohérente (sous condition de relecture). Génération 22 terminée. |
| 24/11/2026 | G22d : invariant de distance ajouté au candidat B2 (numérique ; une loi exacte pour q = 1). |
| 24/11/2026 | G22e : barrière de Cardy (algèbres fixes), partie C. |
| 24/11/2026 | G22a : antériorité du candidat B2 contrôlée sur neuf articles : le théorème C et la question du compromis semblent nouveaux. G22c : résultat négatif (familles structurées), partie C. |
| 24/11/2026 | G21d : un élément ajouté au candidat B2 (pas de témoin de bas degré à N = 3). Génération 21 terminée : aucun nouveau résultat positif. |
| 24/11/2026 | G21b : résultat négatif (caractéristique p), ajouté à la partie C. |
| 24/11/2026 | G21a : résultat négatif (modules), ajouté à la partie C. |
| 24/11/2026 | G21e : la barrière du candidat B2 passe au niveau P+ ; la borne 2/3 de G20c est limitée au cas concave. |
| 24/11/2026 | G21c : nouveau candidat B5 (largeur arborescente après 2026 ; limites de la conjecture du rang asymptotique), classé « remarques » pour le moment. |
| 24/11/2026 | G20b : valeurs exactes de l'annexe B à N fini, ajoutées à la partie C. |
| 24/11/2026 | G20c, G20e : résultats négatifs, ajoutés à la partie C. Aucun nouveau candidat. |
| 24/11/2026 | G20d : une barrière (théorème B, non relue) est ajoutée au candidat B2, qui devient une note possible sous conditions. |
| 24/11/2026 | G20g : le théorème C (candidat B2) est vérifié par Lean ; la bibliothèque B4 a 17 énoncés. |
| 24/11/2026 | G20f : les éléments de G19a (article 3) et de G19b (candidat B2) passent au niveau P+ ; le théorème 5 est corrigé. |
| 24/11/2026 | G19b : compléments au candidat B2 (théorèmes 3 à 5, table exacte). Génération 19 terminée : aucun nouveau candidat d'article. |
| 24/11/2026 | G19d : le théorème A″ (article 3) est vérifié par Lean. |
| 24/11/2026 | G19e : le candidat B2 devient le théorème C (P+). Avis : un lemme, pas un article seul. |
| 24/11/2026 | G19a : quatre éléments ajoutés à l'article 3 (origine des charges, inégalité de séparateur, théorème E, limite de l'aplatissement). |

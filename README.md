# Recherche P vs NP

Journal d'une recherche faite avec Claude Code (novembre 2026), dans la direction **P = NP** :
chercher un algorithme polynomial pour 3-SAT, et comprendre pourquoi chaque voie échoue.

**Résultat honnête : P vs NP n'est pas résolu.** Chaque voie testée échoue, et la cause de chaque
échec est mesurée ou expliquée par un théorème.

## À propos : pourquoi ce projet ?

**Ce projet est une forme de veille informatique sur les IA.** Je ne prétends pas résoudre P vs NP :
c'est l'un des plus grands problèmes ouverts des mathématiques et de l'informatique.

L'idée est de **confier à une IA une tâche a priori impossible, pour voir jusqu'où elle peut aller**.
Un problème ouvert est un bon terrain d'essai : la réponse n'existe nulle part, donc l'IA ne peut pas la recopier.
On voit alors ce qu'elle sait vraiment faire : explorer des pistes, écrire et tester du code, organiser des sous-agents,
rédiger un article, et surtout **reconnaître quand une piste échoue**.

Ce que j'observe en particulier :
- **la fiabilité** : chaque affirmation est vérifiée par des tests chiffrés ou par l'assistant de preuve **Lean 4**, pour ne pas croire l'IA sur parole ;
- **l'honnêteté des résultats** : les échecs sont documentés autant que les réussites ;
- **les limites** : ce que l'IA ne sait pas faire seule, et où un regard humain reste indispensable (les articles ne sont pas encore relus par des experts).

- **`resultats/article/Requena_determinant_bound_clause_product_schemes.pdf`** : l'article (version 2, 26 pages, en anglais ; source `main.tex` à côté). Preuves complètes, section Lean, section sur l'usage de l'IA. **Pas encore relu par des experts humains.**
- **`resultats/article/Requena_borne_determinant_schemas_produit_de_clauses.pdf`** : la version française de l'article (28 pages ; source `main_fr.tex`). La version anglaise fait foi.
- `resultats/brouillons/pre-v3/` : **brouillon** d'une future version 3 (pas l'article officiel), fait pour une vérification de la portée par une IA indépendante. Voir `LISEZ-MOI.md` dans ce dossier.
- `resultats/note_limits_s9.md` (version 1, remplacée par l'article) : note mathématique (en anglais) sur le seul résultat nouveau confirmé, une réponse négative à la question ouverte du §9 de Dias-Lima (arXiv:2611.00517). **Pas encore relue par des experts humains.**
- `lean/` : vérifications Lean 4 + Mathlib, **classées par article** (`lean/README.md`). `PnpLean/Article1/`, `Article2/`, `Article3/` : énoncés vérifiés, sans `sorry` (`cd lean && lake exe cache get && lake build`). Les dossiers `EnCours/` contiennent le travail non terminé.
- `lean/OUTILS.md` : **outils réutilisables** par d'autres chercheurs (23 énoncés génériques vérifiés par Lean, absents de Mathlib), rassemblés dans `lean/PnpLean/Outils.lean`.
- `ARTICLES.md` : **registre des articles possibles** : la matière de chaque article (1, 2, 3), les autres candidats, et ce qui n'est pas assez pour un article.
- `PISTES.md` : toutes les pistes de recherche (en cours, pas encore explorées, fermées). **À consulter en premier pour savoir où on en est.**
- `resultats/rapport_relecture_G8R.md` : rapport de relecture indépendant de la note (recommandation : révision majeure).
- `JOURNAL.md` : toutes les voies explorées, les résultats, les causes des échecs, l'actualité 2026 et les pistes restantes.
- `carnet.html` : la même chose sous forme de page.
- `code/` : programmes et données de chaque génération de sous-agents.
  - `code/bench/` : banc de test commun pour les algorithmes SAT (`python3 bench.py --help`).
  - `code/spec/requirements.md` : liste des familles d'algorithmes déjà éliminées par un théorème.

Les articles cités (Claude/Anthropic, OpenAI, O'Donnell) ne sont pas copiés ici : voir les liens dans `JOURNAL.md`.

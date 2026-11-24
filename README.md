# Recherche P vs NP

Journal d'une recherche faite avec Claude Code (novembre 2026), dans la direction **P = NP** :
chercher un algorithme polynomial pour 3-SAT, et comprendre pourquoi chaque voie échoue.

**Résultat honnête : P vs NP n'est pas résolu.** Chaque voie testée échoue, et la cause de chaque
échec est mesurée ou expliquée par un théorème.

- **`resultats/article/Requena_determinant_bound_clause_product_schemes.pdf`** : l'article (version 2, 26 pages, en anglais ; source `main.tex` à côté). Preuves complètes, section Lean, section sur l'usage de l'IA. **Pas encore relu par des experts humains.**
- **`resultats/article/Requena_borne_determinant_schemas_produit_de_clauses.pdf`** : la version française de l'article (28 pages ; source `main_fr.tex`). La version anglaise fait foi.
- `resultats/note_limits_s9.md` (version 1, remplacée par l'article) : note mathématique (en anglais) sur le seul résultat nouveau confirmé, une réponse négative à la question ouverte du §9 de Dias-Lima (arXiv:2611.00517). **Pas encore relue par des experts humains.**
- `lean/PnpLean/WitnessDet.lean` : preuve Lean 4 + Mathlib du lemme 1 de la note (déterminant témoin et borne de rang), sans `sorry`. Vérifier avec `cd lean && lake exe cache get && lake env lean PnpLean/WitnessDet.lean`.
- `PISTES.md` : toutes les pistes de recherche (en cours, pas encore explorées, fermées). **À consulter en premier pour savoir où on en est.**
- `resultats/rapport_relecture_G8R.md` : rapport de relecture indépendant de la note (recommandation : révision majeure).
- `JOURNAL.md` : toutes les voies explorées, les résultats, les causes des échecs, l'actualité 2026 et les pistes restantes.
- `carnet.html` : la même chose sous forme de page.
- `code/` : programmes et données de chaque génération de sous-agents.
  - `code/bench/` : banc de test commun pour les algorithmes SAT (`python3 bench.py --help`).
  - `code/spec/requirements.md` : liste des familles d'algorithmes déjà éliminées par un théorème.

Les articles cités (Claude/Anthropic, OpenAI, O'Donnell) ne sont pas copiés ici : voir les liens dans `JOURNAL.md`.

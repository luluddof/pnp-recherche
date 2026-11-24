# Pré-version 3 (brouillon)

**Ce document n'est pas l'article officiel.** L'article officiel est la version 2, dans `resultats/article/`. Elle est en relecture chez un expert humain. La version 3 sera faite après son retour.

Ce brouillon sert à faire vérifier la **portée** des affirmations par une IA indépendante.

## Fichiers

- `BROUILLON_pre-v3_Requena_determinant_bound.pdf` : le brouillon (41 pages, en anglais).
- `pre_v3.tex` : la source LaTeX.

## Différences avec la version 2

1. La section Lean est à jour : les deux théorèmes et tous leurs lemmes sont vérifiés par machine.
2. La comparaison avec le théorème 7.1(ii) de Dias et Lima est corrigée. Le texte dit que leur borne est meilleure pour les petits groupes.
3. Une nouvelle section 6 traite la formule F\* (refaite par un second agent IA, pas relue par un humain).
4. Une nouvelle section 11, « Scope of the claims », contient un tableau de portée et 8 questions pour le relecteur.

## Comment le faire vérifier par une autre IA

Donnez-lui **deux fichiers** :
1. ce brouillon (PDF) ;
2. l'article source : Dias et Lima, « Limits and Constants of OpenAI's Claimed ETH Refutation », arXiv:2611.00517.

Exemple de demande :

> Voici un brouillon d'article et l'article source auquel il répond. Ne résume pas le brouillon. Réponds aux 8 questions de la section 11 (« Questions for a scope review »). Pour chaque question, cite la phrase exacte du brouillon et la phrase exacte de l'article source. Cherche les affirmations trop fortes : dans le titre, le résumé, l'introduction et la conclusion. Dis aussi si la définition 5.1 et la question du §9 sont citées fidèlement. Classe chaque problème : erreur, affirmation trop forte, ou imprécision.

## Limites à connaître

- Les agents IA qui ont écrit et relu ce texte sont du même modèle : leurs relectures ne sont pas indépendantes.
- La section 6 suppose que F\* a « une copie par occurrence ». C'est notre lecture de la description de Dias et Lima, pas une citation.
- Les tests numériques portent sur de petites tailles : ils détectent des erreurs, mais ils ne prouvent rien pour les grandes tailles.

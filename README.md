
# NGS-IST Annotation — Streamlit (MVP)

Ce dépôt contient une **appli Streamlit** minimale pour annoter des résultats de métagénomique.

## Contenu

- `app.py` — l’application Streamlit
- `requirements.txt` — dépendances
- `cases_long_with_batch_ctrl_fixedpcr_noaccents.csv` — votre base (à placer à la racine du repo)

## Déployer sur Streamlit Community Cloud

1. **Créez un repo GitHub** et ajoutez `app.py`, `requirements.txt` et votre CSV.
2. Allez sur https://streamlit.io/cloud → New app → sélectionnez votre repo et `app.py`.
3. Lancez l’app.

> MVP : les réponses sont sauvegardées **localement** côté app (système de fichiers éphémère sur le cloud).  
> Chaque opérateur peut **télécharger** son CSV via le bouton “Télécharger les réponses”.
>
> **Étape 2** (à venir) : intégration **Google Sheets** pour une persistance centrale et la collaboration en temps réel.

## Utilisation

1. Saisissez votre **identifiant opérateur** (ex. `op1`).
2. L’app charge la base :
   - par défaut depuis `cases_long_with_batch_ctrl_fixedpcr_noaccents.csv` s’il est présent,
   - sinon via **Upload** (barre latérale).
3. L’app affiche :
   - Contexte (Matrice, Patient, Visite)
   - Cible
   - **Mesures** : RPM espèce/genre, contrôles, **ratio**
   - **Graphique 1** : barres combinées  
     (espèce = **vert**, genre = **rouge**) pour échantillon & contrôle
   - **Graphique 2** : **histogramme batch** bleu  
     + ligne **verte** (échantillon)  
     + lignes **pointillées noires** aux seuils **0**, **1** et **10** RPM.
4. Choisissez un label : `Negatif | Positif faible | Positif | Echec technique | Je ne sais pas`
5. Cliquez **“Enregistrer et suivant”**.
6. Téléchargez vos réponses via le bouton dédié (CSV).

## Champs sauvegardés

- `timestamp`, `operator_id`, `question_id`
- `filename`, `patient`, `visit`, `matrix`, `batch_group`, `pathogen`
- `rpm_species`, `rpm_genus`, `ctrl_env_rpm_species`, `ctrl_env_rpm_genus`, `ratio_species_genus`
- `label`, `notes`

## Notes

- Le **nom de fichier** et le **batch** ne sont **pas affichés** à l’opérateur, mais sont sauvegardés dans le CSV pour la traçabilité.
- La sélection de la “prochaine question” est **aléatoire** parmi celles non encore annotées. On pourra remplacer par un `assignments.csv` si souhaité.
- Pour un déploiement production multi-opérateurs, prévoir **Google Sheets** (ou Postgres/Supabase) pour la persistance.


# NGS-IST — Annotation operateur (Streamlit)

Application Streamlit pour l'annotation manuelle de resultats metagenomiques (NGS-IST) par des operateurs, avec enregistrement centralise dans Google Sheets.

## Contenu du depot
- `app.py` — application Streamlit
- `requirements.txt` — dependances Python
- `cases_questions_v3_1665_with_assignments.csv` — banque de questions (1665 lignes) avec la colonne `operator`
- `manifest_batch_stats.csv` — (option recommande) manifest agrege pour le graphe "Batch distribution"
  (Alternative : `cases_long_v2_noBLSE.csv` si vous preferez recalculer dynamiquement)
- `README.md` — ce document

## Pre-requis
- Un compte de service Google Cloud (JSON) avec acces a Google Sheets.
- Une Google Sheet partagee avec l'email du service account (role Editeur).
  - L'onglet par defaut s'appelle `responses` (l'app creera l'entete si besoin).
- Deploiement sur Streamlit Cloud (ou execution locale).

## Secrets Streamlit (obligatoire)
Dans Settings -> Secrets, collez :
```
[gcp_service_account]
type = "service_account"
project_id = "VOTRE_PROJECT_ID"
private_key_id = "VOTRE_PRIVATE_KEY_ID"
private_key = """-----BEGIN PRIVATE KEY-----
VOTRE_CLE_PRIVEE_MULTILIGNE
-----END PRIVATE KEY-----"""
client_email = "service-account@votre-projet.iam.gserviceaccount.com"
client_id = "VOTRE_CLIENT_ID"
token_uri = "https://oauth2.googleapis.com/token"

# Google Sheet
SHEET_ID = "VOTRE_SHEET_ID"
SHEET_TAB = "responses"

# Utilisateurs
OPERATORS = "op1,op2,op3,op4,op5"
ADMIN_OPERATORS = "admin"
OPERATOR_PINS_JSON = "{\"op1\":\"1234\",\"op2\":\"2345\",\"op3\":\"3456\",\"op4\":\"4567\",\"op5\":\"5678\",\"admin\":\"9999\"}"

# Banque de questions
CASES_PATH = "cases_questions_v3_1665_with_assignments.csv"

# Graphe 2 (choisissez UNE des deux voies)
# 1) Manifest agrege (leger, recommande)
MANIFEST_PATH = "manifest_batch_stats.csv"
LONG_BASE_PATH = ""

# 2) Recalcul dynamique depuis la base longue (identique au comportement valide precedemment)
# MANIFEST_PATH = ""
# LONG_BASE_PATH = "cases_long_v2_noBLSE.csv"
```
Important : le champ private_key doit rester multiligne entre triple guillemets.

### Partager la Google Sheet
1. Ouvrez la feuille (https://sheets.google.com), onglet `responses`.
2. Bouton Partager -> ajoutez l'email du service account (champ `client_email` du JSON) en Editeur.
3. Recuperez l'ID de la feuille : la partie entre `/d/` et `/edit` dans l'URL.

## Deployer sur Streamlit Cloud
1. Creez un repo GitHub avec les fichiers ci-dessus.
2. Sur https://streamlit.io/cloud -> New app -> selectionnez le repo et `app.py` -> Deploy.
3. Dans la sidebar de l'app, admin peut tester la connexion Sheets via le bouton dedie.

## Utilisation (operateurs)
- Choisir Identifiant (ex. `op1`) et entrer le PIN si configure.
- Lire le contexte (matrice/patient/visite/cible/batch).
- Examiner les 2 graphes :
  1. Sample vs Control (Species/Genus) — barres vertes (espece) / rouges (genre), transparence pour les controles.
  2. Batch distribution — histogramme agrege par classes de RPM (0–1, 1–10, 10–100, 100–1000, >=1000), barre verte = valeur de l'echantillon, lignes pointillees aux seuils 0/1/10.
- Selectionner un label : `negatif`, `positif_faible`, `positif`, `echec_technique`, `je_ne_sais_pas`.
- Cliquer Enregistrer & Question suivante.

## Conseils & depannage
- 403 / permission denied : la Google Sheet n'est pas partagee avec l'email du service account (role Editeur).
- Invalid private key : respectez les retours a la ligne entre triple guillemets.
- Worksheet not found : verifiez `SHEET_TAB = "responses"` (ou renommez l'onglet dans la feuille).
- Pas de graphe 2 : renseignez MANIFEST_PATH ou LONG_BASE_PATH et verifiez l'existence du fichier dans le repo.
- Progression : un compteur affiche repondu / assigne pour l'operateur courant.

## Notes techniques
- RPM utilise pour les decisions/graphes :
  - Espece pour les cibles exactes (ex. Neisseria gonorrhoeae).
  - Genre pour les `spp` et Papillomavirus (somme des genres contenant "papillomavirus").

Bonnes annotations !

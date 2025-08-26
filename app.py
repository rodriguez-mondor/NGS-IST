import streamlit as st
import gspread
from google.oauth2.service_account import Credentials
from datetime import datetime, timezone

st.set_page_config(page_title="GSheet test", page_icon="📝", layout="centered")
st.title("Test de connexion Google Sheets (app fictive)")

# Auth via service account JSON placé dans st.secrets["gcp_service_account"]
if "gcp_service_account" not in st.secrets:
    st.error("❌ st.secrets['gcp_service_account'] manquant. Voir README pour la configuration.")
    st.stop()

SCOPES = ["https://www.googleapis.com/auth/spreadsheets"]
creds_info = st.secrets["gcp_service_account"]
sa_email = creds_info.get("client_email", "(inconnu)")
credentials = Credentials.from_service_account_info(creds_info, scopes=SCOPES)
client = gspread.authorize(credentials)

st.caption(f"Service account : **{sa_email}** — partagez la Google Sheet avec cet email (rôle Éditeur).")

# Paramètres (peuvent aussi venir de st.secrets)
SHEET_ID = st.secrets.get("gsheet_id", "")
sheet_id = st.text_input("Google Sheet ID", value=SHEET_ID, placeholder="ex: 1A2B3C...")
ws_name = st.text_input("Worksheet (onglet)", value=st.secrets.get("worksheet_name", "Sheet1"))

st.subheader("Message à écrire")
msg = st.text_area("Votre message", placeholder="Tapez un message de test…", height=120)
operator = st.selectbox("Opérateur", ["op1","op2","op3","op4","op5","admin"], index=0)
write_btn = st.button("Écrire la ligne")

def sanitize_cell(s: str) -> str:
    # Empêche l'interprétation en formule dans Sheets
    if s and s[0] in ("=", "+", "-", "@", "\t"):
        return "'" + s
    return s

if write_btn:
    if not sheet_id.strip():
        st.error("Veuillez renseigner le Google Sheet ID.")
        st.stop()
    try:
        sh = client.open_by_key(sheet_id.strip())
        try:
            ws = sh.worksheet(ws_name)
        except gspread.WorksheetNotFound:
            ws = sh.add_worksheet(title=ws_name, rows=100, cols=10)

        ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
        row = [ts, operator, sanitize_cell(msg)]
        ws.append_row(row, value_input_option="USER_ENTERED")

        st.success(f"✅ Ligne écrite dans '{ws_name}' du document {sheet_id[:6]}…")
        st.write({"timestamp_utc": ts, "operator": operator, "message_preview": (msg[:60] + "…") if len(msg) > 60 else msg})
    except gspread.SpreadsheetNotFound:
        st.error("❌ Feuille introuvable. Avez-vous partagé la feuille avec le service account ?")
        st.info(f"Partagez la feuille avec: **{sa_email}** (rôle: Éditeur)")
    except Exception as e:
        st.exception(e)

st.divider()
with st.expander("Aide rapide"):
    st.markdown("""
    - **Sheet ID** : la chaîne dans l'URL entre `/d/` et `/edit`.
    - **Partage** : ajoutez l'email du service account affiché ci-dessus en **Éditeur**.
    - **Secrets** (facultatif) : `gsheet_id` et `worksheet_name` peuvent être définis dans `st.secrets`.
    - Colonnes ajoutées : **timestamp_utc**, **operator**, **message**.
    """)

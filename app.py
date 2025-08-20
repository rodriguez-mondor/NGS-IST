
import os, json
from datetime import datetime

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import gspread
from google.oauth2.service_account import Credentials


# =========================
# Config depuis Secrets
# =========================
OPERATORS = [x.strip() for x in st.secrets.get("OPERATORS", "op1,op2").split(",") if x.strip()]
ADMIN_OPERATORS = [x.strip() for x in st.secrets.get("ADMIN_OPERATORS", "admin").split(",") if x.strip()]
try:
    OPERATOR_PINS = json.loads(st.secrets.get("OPERATOR_PINS_JSON", "{}"))
except Exception:
    OPERATOR_PINS = {}

CASES_PATH = st.secrets.get("CASES_PATH", "cases_questions_v3_1665_with_assignments.csv")
LONG_BASE_PATH = st.secrets.get("LONG_BASE_PATH", "")  # optionnel pour histogramme batch

SHEET_ID = st.secrets.get("SHEET_ID", "")
SHEET_TAB = st.secrets.get("SHEET_TAB", "responses")

RESPONSE_COLUMNS = [
    "timestamp", "operator", "question_id", "label", "notes",
    "filename", "patient", "visit", "matrix", "pathogen",
    "used_rpm", "used_ctrl_rpm", "batch_group",
    "batch_n_total", "batch_n_positive", "batch_median_rpm_pos",
    "ratio_species_genus", "rpm_band", "pcr_status"
]


# =========================
# Helpers
# =========================
@st.cache_resource
def get_gsheet_ws():
    """Ouvre/initialise la feuille Google Sheets (onglet responses)."""
    scopes = ["https://www.googleapis.com/auth/spreadsheets"]
    creds = Credentials.from_service_account_info(dict(st.secrets["gcp_service_account"]), scopes=scopes)
    gc = gspread.authorize(creds)
    sh = gc.open_by_key(SHEET_ID)
    try:
        ws = sh.worksheet(SHEET_TAB)
    except gspread.WorksheetNotFound:
        ws = sh.add_worksheet(title=SHEET_TAB, rows=1, cols=len(RESPONSE_COLUMNS))
    # En-tête si vide
    if len(ws.get_all_values()) == 0:
        ws.append_row(RESPONSE_COLUMNS, value_input_option="RAW")
    return ws


@st.cache_data
def load_cases():
    """Charge la banque de questions (1665 lignes)."""
    df = pd.read_csv(CASES_PATH)
    need = {
        "question_id", "operator", "filename", "patient", "visit", "matrix", "pathogen",
        "rpm_species", "rpm_genus", "ctrl_env_rpm_species", "ctrl_env_rpm_genus",
        "used_rpm", "used_ctrl_rpm", "ratio_species_genus", "species_or_spp",
        "batch_group", "batch_n_total", "batch_n_positive", "batch_median_rpm_pos",
        "pcr_status", "rpm_band"
    }
    missing = need - set(df.columns)
    if missing:
        st.error(f"Colonnes manquantes dans {CASES_PATH}: {missing}")
    return df


@st.cache_data
def load_long_base():
    """Charge la base longue (optionnelle) pour le graph de distribution batch."""
    if not LONG_BASE_PATH or not os.path.exists(LONG_BASE_PATH):
        return None
    try:
        df = pd.read_csv(LONG_BASE_PATH)
        need = {"batch_group", "matrix", "pathogen", "filename", "rpm_species", "rpm_genus"}
        if not need.issubset(df.columns):
            return None
        return df
    except Exception:
        return None


def used_vals(row):
    """Retourne (used_rpm, used_ctrl_rpm) en fonction de species_vs_spp/papilloma."""
    if str(row.get("species_or_spp", "species")) == "species":
        return float(row.get("rpm_species", 0.0)), float(row.get("ctrl_env_rpm_species", 0.0))
    return float(row.get("rpm_genus", 0.0)), float(row.get("ctrl_env_rpm_genus", 0.0))


def save_to_gsheet(d):
    ws = get_gsheet_ws()
    ws.append_row([d.get(c, "") for c in RESPONSE_COLUMNS], value_input_option="RAW")


def save_to_local(d):
    path = "/tmp/responses_PROD.csv"
    if not os.path.exists(path):
        pd.DataFrame(columns=RESPONSE_COLUMNS).to_csv(path, index=False)
    df = pd.read_csv(path)
    df = pd.concat([df, pd.DataFrame([d], columns=RESPONSE_COLUMNS)], ignore_index=True)
    df.to_csv(path, index=False)


@st.cache_data(ttl=120)
def answered_ids_for(op):
    """Recharge la progression depuis Google Sheets (toutes les 2 minutes)."""
    ws = get_gsheet_ws()
    recs = ws.get_all_records()
    if not recs:
        return set()
    df = pd.DataFrame(recs)
    df = df[df["operator"] == op]
    return set(df["question_id"].astype(str))


def next_pending(my_cases, answered_ids):
    """Donne la prochaine question en attente pour l’opérateur."""
    pending = my_cases[~my_cases["question_id"].astype(str).isin(answered_ids)].copy()
    if pending.empty:
        return None
    pending = pending.sort_values(["pathogen", "matrix", "rpm_band", "filename", "sample_id"])
    return pending.iloc[0].to_dict()


# =========================
# UI
# =========================
st.set_page_config(page_title="NGS-IST Annotation", layout="wide")
st.title("NGS-IST — Annotation opérateur")

cases = load_cases()
long_base = load_long_base()

with st.sidebar:
    st.header("Connexion")
    op = st.selectbox("Identifiant", options=OPERATORS + ADMIN_OPERATORS, index=0)
    pin_needed = OPERATOR_PINS.get(op)
    if pin_needed:
        pin = st.text_input("PIN", type="password", max_chars=12)
        if pin != pin_needed:
            st.warning("Entrez votre PIN pour continuer.")
            st.stop()

    # Mode TEST visible uniquement pour admin (écrit en local /tmp)
    TEST_MODE = st.checkbox("Mode TEST (local)", value=False) if op in ADMIN_OPERATORS else False

    # Bouton test Google Sheets (admin)
    if op in ADMIN_OPERATORS and st.button("Tester connexion Google Sheets"):
        try:
            ws = get_gsheet_ws()
            test = {c: "" for c in RESPONSE_COLUMNS}
            test.update({
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "operator": op,
                "question_id": "TEST",
                "label": "negatif"
            })
            ws.append_row([test.get(c, "") for c in RESPONSE_COLUMNS], value_input_option="RAW")
            st.success("OK : ligne test écrite dans le Google Sheet ✅")
        except Exception as e:
            st.error(f"Erreur : {e}")

# Progression
answered_ids = answered_ids_for(op) if not TEST_MODE else set()
my_cases = cases[cases["operator"] == op].copy()
total_assigned = len(my_cases)
done = len([x for x in answered_ids if x in set(my_cases["question_id"].astype(str))])
st.info(f"**Progression {op} : {done} / {total_assigned}**")

row = next_pending(my_cases, answered_ids)
if row is None:
    st.success("🎉 Vous avez terminé toutes vos questions.")
    st.stop()

# =========================
# Question en cours
# =========================
c1, c2 = st.columns([1, 1])

with c1:
    st.subheader("🧪 Contexte")
    st.markdown(f"- **Matrice** : {row['matrix']}")
    st.markdown(f"- **Patient** : {row['patient']} — **Visite** : {row['visit']}")
    st.markdown(f"- **Batch** : {row['batch_group']}")
    st.markdown(
        f"- **Cible** : <span style='color:green;font-weight:700;font-size:1.05rem'>{row['pathogen']}</span>",
        unsafe_allow_html=True
    )
    used_rpm, used_ctrl = used_vals(row)

with c2:
    st.subheader("Graphiques")

    # Graph 1 : Sample vs Control (species/genre)
    fig1, ax1 = plt.subplots()
    labels = ["Species (sample)", "Species (ctrl)", "Genus (sample)", "Genus (ctrl)"]
    vals = [
        float(row.get("rpm_species", 0.0)),
        float(row.get("ctrl_env_rpm_species", 0.0)),
        float(row.get("rpm_genus", 0.0)),
        float(row.get("ctrl_env_rpm_genus", 0.0)),
    ]
    colors = ["#2ca02c", "#2ca02c", "#d62728", "#d62728"]  # vert espèce, rouge genre
    alphas = [1.0, 0.4, 1.0, 0.4]
    xs = np.arange(len(labels))
    for i, v in enumerate(vals):
        ax1.bar(xs[i], v, alpha=alphas[i], color=colors[i])
    ax1.set_xticks(xs, labels, rotation=15, ha="right")
    ax1.set_ylabel("RPM")
    ax1.set_title("Sample vs Control (Species/Genus)")
    st.pyplot(fig1, clear_figure=True)

    # Graph 2 : Batch distribution (si base longue dispo)
    if long_base is not None:
        dfB = long_base[
            (long_base["batch_group"] == row["batch_group"]) &
            (long_base["matrix"] == row["matrix"]) &
            (long_base["pathogen"] == row['pathogen'])
        ].copy()
        if not dfB.empty:
            def _u(r):
                return float(r["rpm_species"]) if str(row.get("species_or_spp", "species")) == "species" else float(r["rpm_genus"])
            dfB["u"] = dfB.apply(_u, axis=1)

            fig2, ax2 = plt.subplots()
            ax2.hist(dfB["u"].values, bins=[0, 1, 10, 100, 1000], edgecolor="black")
            ax2.axvline(used_rpm, color="#2ca02c", linewidth=3)  # barre verte échantillon
            for thr in [0, 1, 10]:
                ax2.axvline(thr, linestyle="--", color="black", linewidth=1)
            n_total = int(row.get("batch_n_total", len(dfB)))
            ax2.set_title(f"Batch distribution ({row['matrix']}, {row['pathogen']}, N={n_total})")
            ax2.set_xlabel("RPM (used)")
            ax2.set_ylabel("Count")
            st.pyplot(fig2, clear_figure=True)
        else:
            st.caption("Batch distribution indisponible (base longue absente pour ce batch).")
    else:
        st.caption("Batch distribution indisponible (LONG_BASE_PATH non fourni).")

st.markdown("---")
st.subheader("Votre interprétation")

LABELS = ["negatif", "positif_faible", "positif", "echec_technique", "je_ne_sais_pas"]
label = st.radio("Choisissez un label :", options=LABELS, horizontal=True, index=None)
notes = st.text_area("Notes (optionnel)")

btn = st.button("Enregistrer & Question suivante", type="primary", use_container_width=True)

if btn:
    if not label:
        st.warning("Veuillez sélectionner un label.")
        st.stop()

    out = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "operator": op,
        "question_id": str(row["question_id"]),
        "label": label,
        "notes": notes or "",
        "filename": row["filename"],
        "patient": row["patient"],
        "visit": row["visit"],
        "matrix": row["matrix"],
        "pathogen": row["pathogen"],
        "used_rpm": float(used_rpm),
        "used_ctrl_rpm": float(used_ctrl),
        "batch_group": row["batch_group"],
        "batch_n_total": int(row.get("batch_n_total", 0)) if str(row.get("batch_n_total", "")).isdigit() else 0,
        "batch_n_positive": int(row.get("batch_n_positive", 0)) if str(row.get("batch_n_positive", "")).isdigit() else 0,
        "batch_median_rpm_pos": float(row.get("batch_median_rpm_pos", 0.0)),
        "ratio_species_genus": float(row.get("ratio_species_genus", 0.0)),
        "rpm_band": row.get("rpm_band", ""),
        "pcr_status": row.get("pcr_status", "NA"),
    }

    try:
        if TEST_MODE:
            save_to_local(out)
        else:
            if not SHEET_ID:
                st.error("SHEET_ID manquant dans les secrets.")
                st.stop()
            save_to_gsheet(out)
        st.success("Réponse enregistrée ✅")
        st.cache_data.clear()  # rafraîchir progression
        st.experimental_rerun()
    except Exception as e:
        st.error(f"Erreur d’enregistrement : {e}")

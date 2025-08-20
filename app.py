
import os
import io
import json
from datetime import datetime

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ----------------------
# Config & defaults
# ----------------------
DEFAULT_OPERATORS = ["op1","op2","op3","op4","op5"]
DEFAULT_ADMINS = ["admin"]

def get_cfg():
    # Lis depuis Secrets si dispo, sinon valeurs par défaut
    operators = [x.strip() for x in st.secrets.get("OPERATORS", ",".join(DEFAULT_OPERATORS)).split(",") if x.strip()]
    admins = [x.strip() for x in st.secrets.get("ADMIN_OPERATORS", ",".join(DEFAULT_ADMINS)).split(",") if x.strip()]
    pins = {}
    try:
        pins = json.loads(st.secrets.get("OPERATOR_PINS_JSON", "{}"))
    except Exception:
        pins = {}
    cases_path = st.secrets.get("CASES_PATH", "cases_questions_v3_1665_with_assignments.csv")
    manifest_path = st.secrets.get("MANIFEST_PATH", "")  # optionnel
    long_base_path = st.secrets.get("LONG_BASE_PATH", "")  # optionnel
    return {
        "operators": operators,
        "admins": admins,
        "pins": pins,
        "cases_path": cases_path,
        "manifest_path": manifest_path,
        "long_base_path": long_base_path,
    }


RESPONSE_COLUMNS = [
    "timestamp", "operator", "question_id", "label", "notes",
    "filename", "patient", "visit", "matrix", "pathogen",
    "used_rpm", "used_ctrl_rpm", "batch_group",
    "batch_n_total", "batch_n_positive", "batch_median_rpm_pos",
    "ratio_species_genus", "rpm_band", "pcr_status"
]


# ----------------------
# Data loaders
# ----------------------
@st.cache_data
def load_cases(cases_path):
    if not os.path.exists(cases_path):
        st.error(f"Fichier de questions introuvable : {cases_path}. Placez '{cases_path}' à la racine du repo.")
        return pd.DataFrame()
    df = pd.read_csv(cases_path)
    need = {
        "question_id", "operator", "filename", "patient", "visit", "matrix", "pathogen",
        "rpm_species", "rpm_genus", "ctrl_env_rpm_species", "ctrl_env_rpm_genus",
        "used_rpm", "used_ctrl_rpm", "ratio_species_genus", "species_or_spp",
        "batch_group", "batch_n_total", "batch_n_positive", "batch_median_rpm_pos",
        "pcr_status", "rpm_band"
    }
    missing = need - set(df.columns)
    if missing:
        st.error(f"Colonnes manquantes dans {cases_path}: {missing}")
    return df


@st.cache_data
def load_manifest(manifest_path):
    if not manifest_path or not os.path.exists(manifest_path):
        return None
    try:
        df = pd.read_csv(manifest_path)
        need = ["batch_group","matrix","pathogen",
                "count_0_1","count_1_10","count_10_100","count_100_1000","count_ge_1000",
                "n_total","n_positive","median_pos"]
        missing = [c for c in need if c not in df.columns]
        if missing:
            st.warning(f"Colonnes manquantes dans MANIFEST_PATH: {missing}")
        return df
    except Exception as e:
        st.warning(f"Lecture MANIFEST_PATH impossible: {e}")
        return None


@st.cache_data
def load_long_base(long_base_path):
    if not long_base_path or not os.path.exists(long_base_path):
        return None
    try:
        df = pd.read_csv(long_base_path)
        need = {"batch_group", "matrix", "pathogen", "filename", "rpm_species", "rpm_genus"}
        if not need.issubset(df.columns):
            st.warning("Colonnes manquantes dans LONG_BASE_PATH; graph 2 (calcul) desactive.")
            return None
        return df
    except Exception as e:
        st.warning(f"Lecture LONG_BASE_PATH impossible: {e}")
        return None


# ----------------------
# Helpers
# ----------------------
def used_vals(row):
    if str(row.get("species_or_spp", "species")) == "species":
        return float(row.get("rpm_species", 0.0)), float(row.get("ctrl_env_rpm_species", 0.0))
    return float(row.get("rpm_genus", 0.0)), float(row.get("ctrl_env_rpm_genus", 0.0))


def local_test_path():
    return "/tmp/responses_TEST.csv"


def load_local_responses():
    path = local_test_path()
    if not os.path.exists(path):
        return pd.DataFrame(columns=RESPONSE_COLUMNS)
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame(columns=RESPONSE_COLUMNS)


def save_to_local(d):
    path = local_test_path()
    if not os.path.exists(path):
        pd.DataFrame(columns=RESPONSE_COLUMNS).to_csv(path, index=False)
    df = pd.read_csv(path)
    df = pd.concat([df, pd.DataFrame([d], columns=RESPONSE_COLUMNS)], ignore_index=True)
    df.to_csv(path, index=False)


def answered_ids_local(op):
    df = load_local_responses()
    if df.empty:
        return set()
    df = df[df["operator"] == op]
    return set(df["question_id"].astype(str))


def next_pending(my_cases, answered_ids):
    pending = my_cases[~my_cases["question_id"].astype(str).isin(answered_ids)].copy()
    if pending.empty:
        return None
    sort_cols = [c for c in ["pathogen", "matrix", "rpm_band", "filename", "patient"] if c in pending.columns]
    if sort_cols:
        pending = pending.sort_values(sort_cols)
    return pending.iloc[0].to_dict()


# ----------------------
# UI
# ----------------------
st.set_page_config(page_title="NGS-IST — TEST local", layout="wide")
st.title("NGS-IST — Mode TEST (local, sans Google Sheets)")

cfg = get_cfg()
cases = load_cases(cfg["cases_path"])
manifest = load_manifest(cfg["manifest_path"])
long_base = load_long_base(cfg["long_base_path"])

with st.sidebar:
    st.header("Connexion")
    # En mode test, tout le monde peut utiliser la case "Mode TEST"
    TEST_MODE = st.checkbox("Mode TEST (local)", value=True)

    # Liste d'opérateurs toujours avec 5 op par défaut + admin
    all_ops = cfg["operators"] if cfg["operators"] else DEFAULT_OPERATORS
    all_admins = cfg["admins"] if cfg["admins"] else DEFAULT_ADMINS
    op = st.selectbox("Identifiant", options=all_ops + all_admins, index=0)

    # En TEST, si admin choisi, proposer de simuler un opérateur
    sim_op = None
    if TEST_MODE and op in all_admins:
        st.info("Vous êtes en **admin** — sélectionnez un opérateur à simuler pour tester le flux.")
        sim_op = st.selectbox("Simuler l'opérateur", options=all_ops, index=0)

    # Boutons utilitaires pour le mode TEST
    st.markdown("### Outils TEST")
    if st.button("Réinitialiser réponses locales (TEST)"):
        try:
            path = local_test_path()
            if os.path.exists(path):
                os.remove(path)
            st.success("Fichier de réponses locales réinitialisé.")
        except Exception as e:
            st.error(f"Impossible de réinitialiser: {e!r}")

    # Bouton de téléchargement du CSV local
    df_local = load_local_responses()
    if not df_local.empty:
        buf = io.StringIO()
        df_local.to_csv(buf, index=False)
        st.download_button("Télécharger responses_TEST.csv", buf.getvalue(), file_name="responses_TEST.csv", mime="text/csv")

# Détermination de l'opérateur effectif pour filtrer les questions
effective_op = sim_op if (TEST_MODE and op in cfg["admins"] and sim_op) else op

# Filtrer les questions
if cases.empty:
    st.stop()

my_cases = cases.copy()
if effective_op not in cfg["admins"]:
    # En test, on respecte l'assignation si elle existe, sinon on propose tout
    if "operator" in my_cases.columns and effective_op in my_cases["operator"].unique():
        my_cases = my_cases[my_cases["operator"] == effective_op].copy()
    else:
        st.warning("Aucune assignation trouvée pour cet opérateur dans la banque. Affichage de toutes les questions (TEST).")

# Progression locale basée sur responses_TEST.csv
answered_ids = answered_ids_local(effective_op)
total_assigned = len(my_cases)
done = len([x for x in answered_ids if x in set(my_cases["question_id"].astype(str))])
st.info(f"**Progression {effective_op} (TEST) : {done} / {total_assigned}**")

row = next_pending(my_cases, answered_ids)
if row is None:
    st.success("🎉 Vous avez terminé toutes les questions (TEST) pour cet identifiant.")
    st.stop()

# ----------------------
# Question en cours
# ----------------------
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
    colors = ["#2ca02c", "#2ca02c", "#d62728", "#d62728"]
    alphas = [1.0, 0.4, 1.0, 0.4]
    xs = np.arange(len(labels))
    for i, v in enumerate(vals):
        ax1.bar(xs[i], v, alpha=alphas[i], color=colors[i])
    ax1.set_xticks(xs, labels, rotation=15, ha="right")
    ax1.set_ylabel("RPM")
    ax1.set_title("Sample vs Control (Species/Genus)")
    st.pyplot(fig1, clear_figure=True)

    # Graph 2 : Batch distribution
    bins_labels = ["[0,1)", "[1,10)", "[10,100)", "[100,1000)", "≥1000"]
    have_graph2 = False

    manifest = load_manifest(cfg["manifest_path"])
    if manifest is not None:
        row_m = manifest[
            (manifest["batch_group"] == row["batch_group"]) &
            (manifest["matrix"] == row["matrix"]) &
            (manifest["pathogen"] == row["pathogen"])
        ]
        if not row_m.empty:
            r = row_m.iloc[0]
            counts = [
                int(r.get("count_0_1", 0)),
                int(r.get("count_1_10", 0)),
                int(r.get("count_10_100", 0)),
                int(r.get("count_100_1000", 0)),
                int(r.get("count_ge_1000", 0)),
            ]
            fig2, ax2 = plt.subplots()
            xs2 = np.arange(len(bins_labels))
            ax2.bar(xs2, counts, edgecolor="black")
            ax2t = ax2.twiny()
            ax2t.set_xlim(0, 1000)
            ax2t.axvline(used_rpm, color="#2ca02c", linewidth=3)
            for thr in [0, 1, 10]:
                ax2t.axvline(thr, linestyle="--", color="black", linewidth=1)
            n_total = int(r.get("n_total", 0))
            ax2.set_xticks(xs2, bins_labels)
            ax2.set_ylabel("Count")
            ax2.set_title(f"Batch distribution ({row['matrix']}, {row['pathogen']}, N={n_total})")
            st.pyplot(fig2, clear_figure=True)
            have_graph2 = True

    if not have_graph2 and cfg["long_base_path"]:
        dfB = load_long_base(cfg["long_base_path"])
        if dfB is not None:
            dfB = dfB[
                (dfB["batch_group"] == row["batch_group"]) &
                (dfB["matrix"] == row["matrix"]) &
                (dfB["pathogen"] == row["pathogen"])
            ].copy()
            if not dfB.empty:
                def _u(rrow):
                    p = str(rrow["pathogen"])
                    if p.endswith(" spp") or p == "Papillomavirus":
                        return float(rrow["rpm_genus"])
                    return float(rrow["rpm_species"])
                u = dfB.apply(_u, axis=1).values
                counts, _ = np.histogram(u, bins=[0,1,10,100,1000])
                ge1000 = int((u >= 1000).sum())
                counts = list(counts) + [ge1000]
                fig2, ax2 = plt.subplots()
                xs2 = np.arange(len(bins_labels))
                ax2.bar(xs2, counts, edgecolor="black")
                ax2t = ax2.twiny()
                ax2t.set_xlim(0, 1000)
                ax2t.axvline(used_rpm, color="#2ca02c", linewidth=3)
                for thr in [0, 1, 10]:
                    ax2t.axvline(thr, linestyle="--", color="black", linewidth=1)
                n_total = len(u)
                ax2.set_xticks(xs2, bins_labels)
                ax2.set_ylabel("Count")
                ax2.set_title(f"Batch distribution ({row['matrix']}, {row['pathogen']}, N={n_total})")
                st.pyplot(fig2, clear_figure=True)
                have_graph2 = True

    if not have_graph2:
        st.caption("Batch distribution indisponible (ni MANIFEST_PATH, ni LONG_BASE_PATH utilisable).")

st.markdown("---")
st.subheader("Votre interpretation (TEST)")

LABELS = ["negatif", "positif_faible", "positif", "echec_technique", "je_ne_sais_pas"]
label = st.radio("Choisissez un label :", options=LABELS, horizontal=True, index=None)
notes = st.text_area("Notes (optionnel)")

btn = st.button("Enregistrer & passer à la suivante (TEST)", type="primary", use_container_width=True)

if btn:
    if not label:
        st.warning("Veuillez sélectionner un label.")
        st.stop()

    out = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "operator": effective_op,
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
        save_to_local(out)
        st.success("Réponse enregistrée en local ✅ (/tmp/responses_TEST.csv)")
        st.cache_data.clear()
        st.experimental_rerun()
    except Exception as e:
        st.error(f"Erreur d’enregistrement local : {e!r}")

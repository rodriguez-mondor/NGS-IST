# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np

# ==========================
# CONFIG
# ==========================
st.set_page_config(page_title="NGS-IST Dashboard", layout="wide")

# URL GitHub RAW pour manifest_batch_stats
MANIFEST_BATCH_STATS_URL = "https://raw.githubusercontent.com/rodriguez-mondor/NGS-IST/main/manifest_batch_stats.csv"

# ==========================
# CHARGEMENT DES DONNÉES
# ==========================
@st.cache_data(show_spinner=False)
def load_cases_from_uploader(uploaded_file):
    return pd.read_csv(uploaded_file)

@st.cache_data(show_spinner=False)
def load_cases_from_default(path: str):
    return pd.read_csv(path)

@st.cache_data(show_spinner=False)
def load_manifest_batch_stats():
    return pd.read_csv(MANIFEST_BATCH_STATS_URL)

# ==========================
# UI CHARGEMENT
# ==========================
st.sidebar.header("🔧 Chargement des données")
cases_file = st.sidebar.file_uploader("Fichier cases (CSV)", type=["csv"])
default_cases_path = "cases_long_with_batch_ctrl_fixedpcr_noaccents.csv"

cases = None
if cases_file is not None:
    cases = load_cases_from_uploader(cases_file)
else:
    # Essaye de charger depuis le chemin par défaut (comme dans la v4)
    try:
        cases = load_cases_from_default(default_cases_path)
        st.sidebar.info(f"Chargé depuis le fichier par défaut : {default_cases_path}")
    except Exception:
        st.warning("Charge un CSV dans la barre latérale (ou place le fichier par défaut dans le repo).")
        st.stop()

# Chargement manifest_batch_stats depuis GitHub Raw
try:
    manifest_stats = load_manifest_batch_stats()
except Exception as e:
    st.error(f"Impossible de charger manifest_batch_stats depuis GitHub : {e}")
    st.stop()

st.sidebar.success("Données chargées ✔️")

# ==========================
# FONCTIONS
# ==========================
def plot_batch_hist(row: pd.Series, manifest_stats: pd.DataFrame):
    """Construit l'histogramme batch à partir de manifest_batch_stats.

    Utilise les classes pré-agrégées: count_0_1, count_1_10, count_10_100, count_100_1000, count_ge_1000.

    """
    batch = str(row["batch_group"]) if "batch_group" in row else None
    matrix = str(row["matrix"]) if "matrix" in row else None
    pathogen = str(row["pathogen"]) if "pathogen" in row else None

    if not all([batch, matrix, pathogen]):
        st.warning("Colonnes manquantes pour construire l'histogramme (batch_group/matrix/pathogen).")
        return

    subset = manifest_stats[        (manifest_stats["batch_group"] == batch)        & (manifest_stats["matrix"] == matrix)        & (manifest_stats["pathogen"] == pathogen)    ]

    if subset.empty:
        st.warning(f"Aucune donnée batch trouvée dans manifest_batch_stats pour {batch} - {matrix} - {pathogen}")
        return

    stats = subset.iloc[0]

    # Distribution par classes
    counts = {
        "0-1": stats.get("count_0_1", 0),
        "1-10": stats.get("count_1_10", 0),
        "10-100": stats.get("count_10_100", 0),
        "100-1000": stats.get("count_100_1000", 0),
        "≥1000": stats.get("count_ge_1000", 0),
    }
    df_plot = pd.DataFrame(list(counts.items()), columns=["Classe RPM", "Nombre échantillons"])

    fig = px.bar(
        df_plot,
        x="Classe RPM",
        y="Nombre échantillons",
        title=f"Distribution RPM dans le batch ({batch}, {matrix})",
        text_auto=True,
    )
    st.plotly_chart(fig, use_container_width=True)

    # Résumé
    n_total = int(stats.get("n_total", 0))
    n_positive = int(stats.get("n_positive", 0))
    median_pos = float(stats.get("median_pos", 0.0))

    st.markdown(
        f"""

        **Résumé du batch :**

        - Nombre total d'échantillons : {n_total}

        - Nombre positifs : {n_positive}

        - Médiane des positifs : {median_pos:.2f}

        """

    )

# ==========================
# AFFICHAGE SIMPLE (DEMO)
# ==========================
st.title("🧬 NGS-IST — Annotation (v4-1, batch stats)")

# Pour démo: affiche la première ligne et son histogramme batch
try:
    row = cases.iloc[0]
except Exception:
    st.error("Le fichier cases est vide.")
    st.stop()

st.subheader("Échantillon courant (démo)")
st.write({k: row[k] for k in ["patient", "visit", "matrix", "batch_group", "pathogen"] if k in cases.columns})

plot_batch_hist(row, manifest_stats)

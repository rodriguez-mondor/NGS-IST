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

# Chargement principal des données
@st.cache_data
def load_cases(path: str):
    return pd.read_csv(path)

@st.cache_data
def load_manifest_batch_stats():
    return pd.read_csv(MANIFEST_BATCH_STATS_URL)

# ==========================
# LOAD DATA
# ==========================
st.sidebar.header("🔧 Chargement")
cases_file = st.sidebar.file_uploader("Fichier cases (CSV)", type=["csv"])
if cases_file is None:
    st.stop()

cases = load_cases(cases_file)
manifest_stats = load_manifest_batch_stats()

st.sidebar.success("Données chargées")

# ==========================
# UTIL FUNCTIONS
# ==========================
def plot_batch_hist(row, manifest_stats):
    """Construit l'histogramme batch à partir de manifest_batch_stats"""
    batch = row["batch_group"]
    matrix = row["matrix"]
    pathogen = row["pathogen"]

    subset = manifest_stats[        (manifest_stats["batch_group"] == batch)        & (manifest_stats["matrix"] == matrix)        & (manifest_stats["pathogen"] == pathogen)    ]

    if subset.empty:
        st.warning(f"Aucune donnée batch trouvée dans manifest_batch_stats pour {batch} - {matrix} - {pathogen}")
        return

    stats = subset.iloc[0]

    # Distribution par classes
    counts = {
        "0-1": stats["count_0_1"],
        "1-10": stats["count_1_10"],
        "10-100": stats["count_10_100"],
        "100-1000": stats["count_100_1000"],
        "≥1000": stats["count_ge_1000"],
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
    st.markdown(
        f"""
        **Résumé du batch :**
        - Nombre total d'échantillons : {int(stats['n_total'])}
        - Nombre positifs : {int(stats['n_positive'])}
        - Médiane des positifs : {stats['median_pos']:.2f}
        """
    )

# ==========================
# MAIN DISPLAY
# ==========================
st.title("🧬 NGS-IST — Annotation")

# Exemple : on parcourt les lignes du tableau de cas
# Ici j’affiche juste une ligne pour montrer l’histogramme
row = cases.iloc[0]  # en pratique, tu auras ta logique d’itération/annotation
st.write("Échantillon courant :", row.to_dict())

plot_batch_hist(row, manifest_stats)

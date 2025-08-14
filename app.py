
import hashlib
from datetime import datetime

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="NGS-IST Annotation", layout="wide")

# ---------------------------
# Helpers
# ---------------------------

def question_id(row: pd.Series) -> str:
    base = f"{row.get('filename','')}|{row.get('pathogen','')}"
    return hashlib.sha1(base.encode("utf-8")).hexdigest()

def is_species_level(patho: str) -> bool:
    p = str(patho).lower()
    return (" spp" not in p) and ("papillomavirus" not in p)

def rpm_used_col(patho: str) -> str:
    return "rpm_species" if is_species_level(patho) else "rpm_genus"

def pick_next_question(df, answered_ids, operator_filter=None):
    # prefer items not yet labeled by anyone
    remaining = df[~df["question_id"].isin(answered_ids)].copy()
    if operator_filter is not None:
        pass  # kept for future per-operator assignment logic
    if remaining.empty:
        return None
    # simple strategy: shuffle and return first
    return remaining.sample(1, random_state=42).iloc[0]

def plot_four_bars(sample_rpm_species, ctrl_rpm_species, sample_rpm_genus, ctrl_rpm_genus, title):
    labels = ["Sample\nspecies", "Ctrl env.\nspecies", "Sample\ngenus", "Ctrl env.\ngenus"]
    values = [sample_rpm_species, ctrl_rpm_species, sample_rpm_genus, ctrl_rpm_genus]
    colors = ["green", "green", "red", "red"]
    fig, ax = plt.subplots(figsize=(5.8, 3.0))
    ax.bar(labels, values, color=colors)
    ax.set_ylabel("RPM")
    ax.set_title(title)
    st.pyplot(fig)

def plot_batch_hist(batch_vals, sample_val, species_mode=True, title_suffix=""):
    # batch histogram (blue), sample line green, thresholds 0/1/10 as dotted black
    fig, ax = plt.subplots(figsize=(6.0, 3.0))
    ax.hist(batch_vals, bins=20)  # default blue
    ax.axvline(sample_val, color="green", linestyle="-", linewidth=2, label="Sample")
    for thr in [0, 1, 10]:
        ax.axvline(thr, color="black", linestyle=":", linewidth=1)
    ax.set_xlabel(f"RPM ({'species' if species_mode else 'genus'})")
    ax.set_ylabel("Count")
    ax.set_title(f"Batch distribution {title_suffix}")
    st.pyplot(fig)

# ---------------------------
# Sidebar: inputs & setup
# ---------------------------

st.sidebar.header("Configuration")

operator_id = st.sidebar.text_input("Votre identifiant opérateur", value="", help="Utilisez un ID court (ex. op1, op2...).")
base_csv_source = st.sidebar.radio("Source de la base", ["Fichier dans le repo", "Upload manuel"], index=0)

default_path = "cases_long_with_batch_ctrl_fixedpcr_noaccents.csv"
uploaded = None
if base_csv_source == "Upload manuel":
    uploaded = st.sidebar.file_uploader("Charger la base CSV", type=["csv"])

responses_file = st.sidebar.text_input("Nom du fichier de réponses (local)", value="responses.csv",
                                       help="Les réponses sont aussi proposées en téléchargement.")

st.sidebar.markdown("---")
st.sidebar.caption("Version MVP : stockage local & téléchargement. Une intégration Google Sheets pourra être ajoutée ensuite pour la collaboration en temps réel.")

# ---------------------------
# Load data
# ---------------------------

@st.cache_data(show_spinner=False)
def load_cases(path_or_buffer):
    df = pd.read_csv(path_or_buffer)
    # basic sanity
    needed = {"sample_id","filename","patient","visit","matrix","batch_group","pathogen",
              "rpm_species","rpm_genus","ctrl_env_rpm_species","ctrl_env_rpm_genus"}
    missing = needed - set(df.columns)
    if missing:
        st.error(f"Colonnes manquantes dans la base: {missing}")
        return None
    # compute helper IDs
    df = df.copy()
    df["question_id"] = df.apply(question_id, axis=1)
    return df

if uploaded is not None:
    cases = load_cases(uploaded)
else:
    try:
        cases = load_cases(default_path)
    except Exception as e:
        cases = None
        st.error(f"Impossible de charger {default_path}. Uploadez la base dans la barre latérale.")

if cases is None:
    st.stop()

# ---------------------------
# Load/save responses (local CSV)
# ---------------------------

def load_responses(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
        if "question_id" not in df.columns:
            # backward compatibility
            df["question_id"] = df.apply(lambda r: hashlib.sha1(f"{r.get('filename','')}|{r.get('pathogen','')}".encode("utf-8")).hexdigest(), axis=1)
        return df
    except FileNotFoundError:
        return pd.DataFrame(columns=[
            "timestamp","operator_id","question_id","filename","patient","visit","matrix","batch_group","pathogen",
            "rpm_species","rpm_genus","ctrl_env_rpm_species","ctrl_env_rpm_genus","ratio_species_genus","label","notes"
        ])

def append_response(path: str, row: dict):
    # naive append; in cloud the filesystem is ephemeral, but ok for MVP
    df = load_responses(path)
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    df.to_csv(path, index=False)
    return df

responses = load_responses(responses_file)

# ---------------------------
# Header / progress
# ---------------------------

st.title("NGS-IST — Annotation opérateur (MVP)")

if not operator_id:
    st.info("Entrez votre **identifiant opérateur** dans la barre latérale pour commencer.")
    st.stop()

already = set(responses["question_id"].unique())
total = cases["question_id"].nunique()
done = len(already)
st.write(f"**Progression globale** : {done}/{total} questions annotées.")

# ---------------------------
# Next question
# ---------------------------

row = pick_next_question(cases, already, operator_filter=operator_id)
if row is None:
    st.success("Il n'y a plus de questions à annoter. Merci !")
    st.download_button("Télécharger vos réponses", data=load_responses(responses_file).to_csv(index=False).encode("utf-8"),
                       file_name=responses_file, mime="text/csv")
    st.stop()

# ---------------------------
# Display context (no filename/batch to operator)
# ---------------------------
colA, colB = st.columns([1,1])
with colA:
    st.subheader("Contexte")
    st.markdown(f"- **Matrice** : {row['matrix']}")
    st.markdown(f"- **Patient** : {row['patient']}")
    st.markdown(f"- **Visite** : {row['visit']}")
with colB:
    st.subheader("Cible")
    st.markdown(f"**{row['pathogen']}**")

# KPIs
ratio = (row["rpm_species"] / row["rpm_genus"]) if row["rpm_genus"] > 0 else 0.0
k1, k2, k3, k4, k5 = st.columns(5)
k1.metric("RPM espèce (échantillon)", f"{row['rpm_species']:.3f}")
k2.metric("RPM espèce (ctrl env.)", f"{row['ctrl_env_rpm_species']:.3f}")
k3.metric("RPM genre (échantillon)", f"{row['rpm_genus']:.3f}")
k4.metric("RPM genre (ctrl env.)", f"{row['ctrl_env_rpm_genus']:.3f}")
k5.metric("Ratio espèce/genre", f"{ratio:.3f}")

st.markdown("---")

# ---------------------------
# Graphs
# ---------------------------
g1, g2 = st.columns([1,1])

with g1:
    plot_four_bars(
        row["rpm_species"], row["ctrl_env_rpm_species"],
        row["rpm_genus"], row["ctrl_env_rpm_genus"],
        title="RPMs — échantillon vs contrôle"
    )

with g2:
    # Prepare batch series from current dataframe (same batch_group, matrix & pathogen)
    species_mode = is_species_level(row["pathogen"])
    used_col = "rpm_species" if species_mode else "rpm_genus"
    batch_vals = cases[
        (cases["batch_group"] == row["batch_group"]) &
        (cases["matrix"] == row["matrix"]) &
        (cases["pathogen"] == row["pathogen"])
    ][used_col].astype(float).values
    plot_batch_hist(batch_vals, row[used_col], species_mode=species_mode,
                    title_suffix=f"({row['matrix']}, {row['pathogen']})")

st.markdown("---")

# ---------------------------
# Answer form
# ---------------------------
label = st.radio("Votre évaluation", ["Negatif","Positif faible","Positif","Echec technique","Je ne sais pas"], horizontal=True)
notes = st.text_input("Commentaire (optionnel)")

c1, c2, c3 = st.columns([1,1,1])
with c1:
    if st.button("⏭️ Passer (sans enregistrer)"):
        st.experimental_rerun()

with c2:
    if st.button("💾 Enregistrer et suivant"):
        record = {
            "timestamp": datetime.utcnow().isoformat(),
            "operator_id": operator_id,
            "question_id": row["question_id"],
            "filename": row["filename"],
            "patient": row["patient"],
            "visit": row["visit"],
            "matrix": row["matrix"],
            "batch_group": row["batch_group"],
            "pathogen": row["pathogen"],
            "rpm_species": row["rpm_species"],
            "rpm_genus": row["rpm_genus"],
            "ctrl_env_rpm_species": row["ctrl_env_rpm_species"],
            "ctrl_env_rpm_genus": row["ctrl_env_rpm_genus"],
            "ratio_species_genus": ratio,
            "label": label,
            "notes": notes
        }
        append_response(responses_file, record)
        st.experimental_rerun()

with c3:
    st.download_button("⬇️ Télécharger les réponses", data=load_responses(responses_file).to_csv(index=False).encode("utf-8"),
                       file_name=responses_file, mime="text/csv")

# Debug / details (optional)
with st.expander("Détails techniques (debug)"):
    st.write({k: row[k] for k in ["filename","batch_group"]})
    st.dataframe(load_responses(responses_file).tail(10))

# -*- coding: utf-8 -*-
import os
import hashlib
import tempfile
from datetime import datetime

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title='NGS-IST Annotation', layout='wide')

# ---------------------------
# Helpers
# ---------------------------

def question_id(row: pd.Series) -> str:
    base = f"{row.get('filename','')}|{row.get('pathogen','')}"
    return hashlib.sha1(base.encode('utf-8')).hexdigest()

def is_species_level(patho: str) -> bool:
    p = str(patho).lower()
    return (' spp' not in p) and ('papillomavirus' not in p)

def rpm_used_col(patho: str) -> str:
    return 'rpm_species' if is_species_level(patho) else 'rpm_genus'

def pick_next_question(df, answered_ids):
    remaining = df[~df['question_id'].isin(answered_ids)].copy()
    if remaining.empty:
        return None
    # simple random pick
    return remaining.sample(1).iloc[0]

def get_current_row(cases: pd.DataFrame, answered_ids: set):
    """Keep the same question across reruns until user saves or skips."""
    qid = st.session_state.get('current_qid')
    if qid and qid not in answered_ids:
        row = cases.loc[cases['question_id'] == qid]
        if not row.empty:
            return row.iloc[0]
    # Need to select a new question
    row = pick_next_question(cases, answered_ids)
    if row is not None:
        st.session_state['current_qid'] = row['question_id']
    return row

def clear_current_row():
    if 'current_qid' in st.session_state:
        del st.session_state['current_qid']

def plot_four_bars(sample_rpm_species, ctrl_rpm_species, sample_rpm_genus, ctrl_rpm_genus, title, fs=9):
    labels = ['Sample\nspecies', 'Ctrl env.\nspecies', 'Sample\ngenus', 'Ctrl env.\ngenus']
    values = [sample_rpm_species, ctrl_rpm_species, sample_rpm_genus, ctrl_rpm_genus]
    colors = ['green', 'green', 'red', 'red']
    fig, ax = plt.subplots(figsize=(5.8, 3.0))
    ax.bar(labels, values, color=colors)
    ax.set_ylabel('RPM', fontsize=fs)
    ax.set_title(title, fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    st.pyplot(fig)

def plot_batch_hist(batch_vals, sample_val, species_mode=True, title_suffix='', fs=9):
    # Histogram in default blue; sample line green; thresholds 0/1/10 dotted black
    fig, ax = plt.subplots(figsize=(6.0, 3.0))
    ax.hist(batch_vals, bins=20)
    ax.axvline(sample_val, color='green', linestyle='-', linewidth=2, label='Sample')
    for thr in [0, 1, 10]:
        ax.axvline(thr, color='black', linestyle=':', linewidth=1)
    ax.set_xlabel(f"RPM ({'species' if species_mode else 'genus'})", fontsize=fs)
    ax.set_ylabel('Count', fontsize=fs)
    ax.set_title(f'Batch distribution {title_suffix}', fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    st.pyplot(fig)

# ---------------------------
# Sidebar: inputs & setup
# ---------------------------

st.sidebar.header('Configuration')

operator_id = st.sidebar.text_input('Votre identifiant operateur', value='', help='Utilisez un ID court (ex. op1, op2...).')
base_csv_source = st.sidebar.radio('Source de la base', ['Fichier dans le repo', 'Upload manuel'], index=0)

default_path = 'cases_long_with_batch_ctrl_fixedpcr_noaccents.csv'
uploaded = None
if base_csv_source == 'Upload manuel':
    uploaded = st.sidebar.file_uploader('Charger la base CSV', type=['csv'])

# Mode (Test/Prod) to isolate files
mode = st.sidebar.selectbox('Mode', ['Test', 'Production'], index=0)
default_resp_name = 'responses_TEST.csv' if mode == 'Test' else 'responses.csv'

responses_name = st.sidebar.text_input('Nom du fichier de reponses', value=default_resp_name,
                                       help='Sauvegarde dans un repertoire temporaire du serveur.')
# Always write to a guaranteed-writable temp directory
RESPONSES_PATH = os.path.join(tempfile.gettempdir(), responses_name)

st.sidebar.markdown(f"**Chemin d'ecriture** : `{RESPONSES_PATH}`")
st.sidebar.caption('La zone temporaire est volatile. Telechargez vos reponses pour les conserver localement.')

# Optional: Reset button for TEST mode
if mode == 'Test':
    if st.sidebar.button('🧹 Reset fichier de test'):
        try:
            if os.path.exists(RESPONSES_PATH):
                os.remove(RESPONSES_PATH)
            clear_current_row()
            st.sidebar.success('Fichier de reponses de test supprime.')
        except Exception as e:
            st.sidebar.error(f'Echec du reset: {e}')

# ---------------------------
# Load data
# ---------------------------

@st.cache_data(show_spinner=False)
def load_cases(path_or_buffer):
    df = pd.read_csv(path_or_buffer)
    needed = {'sample_id','filename','patient','visit','matrix','batch_group','pathogen',
              'rpm_species','rpm_genus','ctrl_env_rpm_species','ctrl_env_rpm_genus'}
    missing = needed - set(df.columns)
    if missing:
        st.error(f'Colonnes manquantes dans la base: {missing}')
        return None
    df = df.copy()
    df['question_id'] = df.apply(question_id, axis=1)
    return df

if uploaded is not None:
    cases = load_cases(uploaded)
else:
    try:
        cases = load_cases(default_path)
    except Exception as e:
        cases = None
        st.error(f'Impossible de charger {default_path}. Uploadez la base dans la barre laterale.')

if cases is None:
    st.stop()

# ---------------------------
# Load/save responses (temp dir)
# ---------------------------

def load_responses(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
        if 'question_id' not in df.columns:
            df['question_id'] = df.apply(lambda r: hashlib.sha1(f"{r.get('filename','')}|{r.get('pathogen','')}".encode('utf-8')).hexdigest(), axis=1)
        return df
    except FileNotFoundError:
        return pd.DataFrame(columns=[
            'timestamp','operator_id','question_id','filename','patient','visit','matrix','batch_group','pathogen',
            'rpm_species','rpm_genus','ctrl_env_rpm_species','ctrl_env_rpm_genus','ratio_species_genus','label','notes'
        ])

def append_response(path: str, row: dict):
    try:
        df = load_responses(path)
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
        # Ensure dir exists
        os.makedirs(os.path.dirname(path), exist_ok=True)
        df.to_csv(path, index=False)
        return df, None
    except Exception as e:
        return None, str(e)

responses = load_responses(RESPONSES_PATH)

# ---------------------------
# Header / progress
# ---------------------------

st.title('NGS-IST — Annotation operateur (MVP)')

if not operator_id:
    st.info("Entrez votre **identifiant operateur** dans la barre laterale pour commencer.")
    st.stop()

already = set(responses['question_id'].unique())
total = cases['question_id'].nunique()
done = len(already)

# Operator-specific progression
responses_by_op = responses[responses['operator_id'] == operator_id]
done_by_op = responses_by_op['question_id'].nunique()
remaining_global = total - done

pc1, pc2, pc3 = st.columns(3)
pc1.metric('Progression globale', f'{done}/{total}')
pc2.metric('Vos reponses', f'{done_by_op}/{total}')
pc3.metric('Restantes (global)', f'{remaining_global}')

# ---------------------------
# Current question (locked)
# ---------------------------

row = get_current_row(cases, already)
if row is None:
    st.success("Il n'y a plus de questions a annoter. Merci !")
    st.download_button('Telecharger vos reponses', data=load_responses(RESPONSES_PATH).to_csv(index=False).encode('utf-8'),
                       file_name=os.path.basename(RESPONSES_PATH), mime='text/csv')
    st.stop()

# ---------------------------
# Display context (no filename/batch to operator)
# ---------------------------
colA, colB = st.columns([1,1])
with colA:
    st.subheader('Contexte')
    st.markdown(f"- **Matrice** : {row['matrix']}")
    st.markdown(f"- **Patient** : {row['patient']}")
    st.markdown(f"- **Visite** : {row['visit']}")
with colB:
    st.subheader('Cible')
    st.markdown(f"<h3 style='color:green; font-weight:700; margin:0'>{row['pathogen']}</h3>", unsafe_allow_html=True)

# KPIs
ratio = (row['rpm_species'] / row['rpm_genus']) if row['rpm_genus'] > 0 else 0.0
k1, k2, k3, k4, k5 = st.columns(5)
k1.metric('RPM espece (echantillon)', f"{row['rpm_species']:.3f}")
k2.metric('RPM espece (ctrl env.)', f"{row['ctrl_env_rpm_species']:.3f}")
k3.metric('RPM genre (echantillon)', f"{row['rpm_genus']:.3f}")
k4.metric('RPM genre (ctrl env.)', f"{row['ctrl_env_rpm_genus']:.3f}")
k5.metric('Ratio espece/genre', f'{ratio:.3f}')

st.markdown('---')

# ---------------------------
# Graphs
# ---------------------------
g1, g2 = st.columns([1,1])

with g1:
    plot_four_bars(
        row['rpm_species'], row['ctrl_env_rpm_species'],
        row['rpm_genus'], row['ctrl_env_rpm_genus'],
        title='RPMs — echantillon vs controle',
        fs=9
    )

with g2:
    species_mode = is_species_level(row['pathogen'])
    used_col = 'rpm_species' if species_mode else 'rpm_genus'

    # 1) Tous les échantillons du batch pour la même matrice
    batch_samples = cases[
        (cases['batch_group'].astype(str) == str(row['batch_group'])) &
        (cases['matrix'].astype(str) == str(row['matrix']))
    ]['sample_id'].unique()

    # 2) RPM de la cible pour les échantillons où elle est présente
    sub = cases[
        (cases['batch_group'].astype(str) == str(row['batch_group'])) &
        (cases['matrix'].astype(str) == str(row['matrix'])) &
        (cases['pathogen'].astype(str) == str(row['pathogen']))
    ][['sample_id', used_col]]

    rpm_map = dict(zip(sub['sample_id'], sub[used_col].astype(float)))

    # 3) Série complète incluant les zéros manquants
    import numpy as np
    batch_vals = np.array([float(rpm_map.get(sid, 0.0)) for sid in batch_samples], dtype=float)

    title_sfx = f"({row['matrix']}, {row['pathogen']}, N={n_batch})"
    plot_batch_hist(batch_vals, float(row[used_col]), species_mode=species_mode, title_suffix=title_sfx)


st.markdown('---')

# ---------------------------
# Answer form
# ---------------------------
radio_key = f"radio_{row['question_id']}"
label = st.radio('Votre evaluation', ['Negatif','Positif faible','Positif','Echec technique','Je ne sais pas'],
                 horizontal=True, key=radio_key)
notes = st.text_input('Commentaire (optionnel)', key=f"notes_{row['question_id']}")

c1, c2, c3 = st.columns([1,1,1])
with c1:
    if st.button('⏭️ Passer (sans enregistrer)'):
        clear_current_row()
        st.rerun()

with c2:
    if st.button('💾 Enregistrer et suivant'):
        record = {
            'timestamp': datetime.utcnow().isoformat(),
            'operator_id': operator_id,
            'question_id': row['question_id'],
            'filename': row['filename'],
            'patient': row['patient'],
            'visit': row['visit'],
            'matrix': row['matrix'],
            'batch_group': row['batch_group'],
            'pathogen': row['pathogen'],
            'rpm_species': row['rpm_species'],
            'rpm_genus': row['rpm_genus'],
            'ctrl_env_rpm_species': row['ctrl_env_rpm_species'],
            'ctrl_env_rpm_genus': row['ctrl_env_rpm_genus'],
            'ratio_species_genus': ratio,
            'label': label,
            'notes': notes
        }
        _, err = append_response(RESPONSES_PATH, record)
        if err:
            st.error(f"Erreur d'ecriture du fichier de reponses : {err}")
        clear_current_row()
        st.rerun()

with c3:
    st.download_button('⬇️ Telecharger les reponses',
                       data=load_responses(RESPONSES_PATH).to_csv(index=False).encode('utf-8'),
                       file_name=os.path.basename(RESPONSES_PATH), mime='text/csv')

# Debug / details (optional)
with st.expander('Details techniques (debug)'):
    st.write({k: row[k] for k in ['filename','batch_group']})
    st.write(f'Fichier de reponses : {RESPONSES_PATH}')
    st.dataframe(load_responses(RESPONSES_PATH).tail(10))

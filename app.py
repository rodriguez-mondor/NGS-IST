# -*- coding: utf-8 -*-
import os
import hashlib
import tempfile
from datetime import datetime

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import gspread
from google.oauth2.service_account import Credentials
from matplotlib.ticker import FixedLocator, FixedFormatter, MaxNLocator



# ===== Added: manifest_batch_stats loader and helpers =====
LOCAL_MANIFEST = "manifest_batch_pathogen_stats.csv"  # nouveau manifest complet (batch×site×pathogène)
MANIFEST_BATCH_STATS_URL = "https://raw.githubusercontent.com/rodriguez-mondor/NGS-IST/main/manifest_batch_stats.csv"  # fallback


def load_manifest_batch_stats():
    """Charge le manifest depuis le fichier local fixé en dur.
    Normalise les colonnes pour compatibilité avec l'app.
    """
    import pandas as pd
    MANIFEST_LOCAL = "manifest_batch_pathogen_stats.csv"
    try:
        df = pd.read_csv(MANIFEST_LOCAL)
    except Exception as e:
        st.error(f"Impossible de charger le manifest local '{MANIFEST_LOCAL}': {e}")
        return pd.DataFrame()
    # Normalisation colonnes
    if 'matrix' not in df.columns and 'site' in df.columns:
        df['matrix'] = df['site']
    if 'batch_group' not in df.columns and {'batch_id','batch_matrix'}.issubset(df.columns):
        df['batch_group'] = df['batch_id'].astype(str) + df['batch_matrix'].astype(str)
    for c in ['batch_number','n_total','n_positive']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

def expand_counts_to_values(stats_row):
    # Map des bins -> valeur représentative pour histogramme
    mids = {
        "count_0_1": 0.5,
        "count_1_10": 3.0,
        "count_10_100": 30.0,
        "count_100_1000": 300.0,
        "count_ge_1000": 3000.0,
    }
    vals = []
    for k, m in mids.items():
        c = int(stats_row.get(k, 0))
        if c > 0:
            vals.extend([m] * c)
    return np.array(vals, dtype=float)
# ===== End added =====
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
    # Histogramme (esthétique): axe X en symlog 0→1000, ticks fixes sans '0' doublé
    edges = [0, 1, 10, 100, 1000, 10000]  # bins inchangés (source des données)
    fig, ax = plt.subplots(figsize=(6.2, 3.2))
    # On récupère les comptes pour ajuster l'axe Y proprement
    n, _bins, _patches = ax.hist(batch_vals, bins=edges)
    # Ligne de l'échantillon: vert si species, rouge si genus/spp
    ax.axvline(sample_val, color=('green' if species_mode else 'red'), linestyle='-', linewidth=2)
    # Repères visuels utiles
    for thr in [1, 10]:
        ax.axvline(thr, color='black', linestyle=':', linewidth=1)

    # Axe X en symlog + limites
    ax.set_xscale('symlog', linthresh=0.1)
    ax.set_xlim(0, 1000)

    # Ticks majeurs fixes (sans 0), labels explicites pour éviter le doublon "0"
    major_ticks = [0.1, 1, 10, 100, 1000]
    major_labels = ['0.1', '1', '10', '100', '1000']
    ax.xaxis.set_major_locator(FixedLocator(major_ticks))
    ax.xaxis.set_major_formatter(FixedFormatter(major_labels))
    # Tick mineur non étiqueté à 0 (permet d'afficher 0 sans label)
    ax.xaxis.set_minor_locator(FixedLocator([0]))

    # Axe Y: min à 25, sinon adapté au max des comptes; ticks entiers
    y_max = int((max(n) if len(n) else 0))
    # arrondi au multiple de 5 supérieur
    import math
    y_target = max(25, int(math.ceil(y_max / 5.0) * 5))
    ax.set_ylim(0, y_target if y_target > 0 else 25)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Libellés
    ax.set_xlabel(f"log RPM ({'species' if species_mode else 'genus'})", fontsize=fs)
    ax.set_ylabel('Count', fontsize=fs)
    ax.set_title(f'Batch distribution {title_suffix}', fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    st.pyplot(fig)
# ---------------------------

# Sidebar: inputs & setup
# ---------------------------

st.sidebar.header('Configuration')

operator_choices = ['— choisir —','op1','op2','op3','op4','op5','admin']
operator_id = st.sidebar.selectbox("Votre identifiant opérateur", operator_choices, index=0)
if operator_id == '— choisir —':
    st.sidebar.warning("Merci de choisir votre identifiant opérateur pour commencer.")
    st.stop()
if operator_id == 'admin':
    admin_mode = st.sidebar.radio('Mode admin', ['Test','Supervision'], index=0)
else:
    _mode = st.sidebar.radio('Choisir le mode :', ['Production','Test'], index=1,
                              help="En mode Test, rien n'est enregistré. En Production, écriture Google Sheet.")
    mode = _mode

if operator_id == 'admin':
    base_csv_source = st.sidebar.radio('Source de la base', ['Fichier dans le repo', 'Upload manuel'], index=0)
    default_path = 'cases_questions_v3_1665_with_assignments_withPCR_v3.csv'
    uploaded = None
    if base_csv_source == 'Upload manuel':
        uploaded = st.sidebar.file_uploader('Charger la base CSV', type=['csv'])
else:
    base_csv_source = 'Fichier dans le repo'
    default_path = 'cases_questions_v3_1665_with_assignments_withPCR_v3.csv'
    uploaded = None

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

# Optional: Manifest upload (for complete batch N including zeros)
st.sidebar.markdown('---')
st.sidebar.subheader('Manifest (optionnel)')
st.sidebar.caption("CSV attendu : sample_id,batch_group,matrix — pour compter tous les echantillons du batch, y compris ceux a 0.")
manifest = None
manifest_file = st.sidebar.file_uploader('Charger manifest.csv', type=['csv'])
if manifest_file is None:
    # try to load manifest.csv from repo if present
    if os.path.exists('manifest.csv'):
        try:
            manifest = pd.read_csv('manifest.csv', dtype=str)
        except Exception:
            manifest = None
else:
    try:
        manifest = pd.read_csv(manifest_file, dtype=str)
    except Exception as e:
        st.sidebar.error(f"Manifest illisible : {e}")
        manifest = None

if manifest is not None:
    for col in ['sample_id','batch_group','matrix']:
        if col not in manifest.columns:
            st.sidebar.error(f"Manifest invalide : colonne manquante {col}")
            manifest = None
            break
    if manifest is not None:
        manifest['sample_id'] = manifest['sample_id'].astype(str)
        manifest['batch_group'] = manifest['batch_group'].astype(str)
        manifest['matrix'] = manifest['matrix'].astype(str)

# ---------------------------
# Load data

# ---------------------------



def load_cases(path_or_buffer=None):
    """Charge le fichier des questions depuis un fichier local fixé en dur.
    Normalise les colonnes (site, operator_id, types) pour compatibilité avec l'app.
    """
    import pandas as pd
    CASES_LOCAL = "cases_questions_v3_1665_with_assignments_withPCR_v3.csv"
    try:
        df = pd.read_csv(CASES_LOCAL)
    except Exception as e:
        st.error(f"Impossible de charger le fichier cases local '{CASES_LOCAL}': {e}")
        return pd.DataFrame()
    # Normalisations colonnes
    if "site" not in df.columns and "matrix" in df.columns:
        df["site"] = df["matrix"]
    if "operator_id" not in df.columns and "operator" in df.columns:
        df["operator_id"] = df["operator"]
    for c in ["question_id","sample_id","filename","patient","visit","matrix","site","batch_group","pathogen"]:
        if c in df.columns:
            df[c] = df[c].astype(str)
    for c in ["rpm_species","rpm_genus","ctrl_env_rpm_species","ctrl_env_rpm_genus"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    return df


# === Chargement des 'cases' (questions) ===
if uploaded is not None:
    cases = load_cases(uploaded)
else:
    cases = load_cases(default_path)

if cases is None or cases.empty:
    st.error("Aucune donnée 'cases' chargée. Vérifiez le CSV.")
    st.stop()

# Normalisation minimale
if 'question_id' not in cases.columns:
    cases['question_id'] = cases.apply(question_id, axis=1)
for c in ['question_id','sample_id','filename','patient','visit','matrix','site','batch_group','pathogen']:
    if c in cases.columns:
        cases[c] = cases[c].astype(str)

# === Chargement du manifest statistiques (batch × site × pathogène) ===
try:
    manifest_stats = load_manifest_batch_stats()
    if manifest_stats is not None and not manifest_stats.empty:
        # Harmonise 'matrix' si nécessaire (certains manifests utilisent 'site')
        if 'matrix' not in manifest_stats.columns and 'site' in manifest_stats.columns:
            manifest_stats['matrix'] = manifest_stats['site']
        # Cast types utiles
        for c in ['batch_number','n_total','n_positive']:
            if c in manifest_stats.columns:
                manifest_stats[c] = pd.to_numeric(manifest_stats[c], errors='coerce')
        # Assure str sur clés
        for c in ['batch_group','matrix','pathogen']:
            if c in manifest_stats.columns:
                manifest_stats[c] = manifest_stats[c].astype(str)
    else:
        manifest_stats = None
except Exception as _e:
    manifest_stats = None





# === Google Sheet helpers ===
@st.cache_data(ttl=60)
def load_responses_gsheet():
    try:
        creds_info = st.secrets.get("gcp_service_account", None)
        if not creds_info:
            return pd.DataFrame()
        scopes = ["https://www.googleapis.com/auth/spreadsheets.readonly"]
        credentials = Credentials.from_service_account_info(creds_info, scopes=scopes)
        client = gspread.authorize(credentials)
        sheet_id = (st.secrets.get("gsheet_id", "") or st.secrets.get("SHEET_ID", "")).strip()
        if not sheet_id:
            return pd.DataFrame()
        ws_name = st.secrets.get("worksheet_name", st.secrets.get("SHEET_TAB", "responses"))
        sh = client.open_by_key(sheet_id)
        ws = sh.worksheet(ws_name)
        rows = ws.get_all_records()
        if not rows:
            return pd.DataFrame()
        df = pd.DataFrame(rows)
        return df
    except Exception as e:
        st.warning(f"Lecture Google Sheet impossible: {e}")
        return pd.DataFrame()

def append_response_gsheet(record: dict):
    creds_info = st.secrets.get("gcp_service_account", None)
    if not creds_info:
        raise RuntimeError("st.secrets['gcp_service_account'] manquant.")
    scopes = ["https://www.googleapis.com/auth/spreadsheets"]
    credentials = Credentials.from_service_account_info(creds_info, scopes=scopes)
    client = gspread.authorize(credentials)
    sheet_id = (st.secrets.get("gsheet_id", "") or st.secrets.get("SHEET_ID", "")).strip()
    if not sheet_id:
        raise RuntimeError("ID Google Sheet absent (gsheet_id / SHEET_ID).")
    ws_name = st.secrets.get("worksheet_name", st.secrets.get("SHEET_TAB", "responses"))
    sh = client.open_by_key(sheet_id)
    try:
        ws = sh.worksheet(ws_name)
    except gspread.WorksheetNotFound:
        ws = sh.add_worksheet(title=ws_name, rows=1000, cols=20)
    def _sanitize_cell(s):
        try:
            s = "" if s is None else str(s)
        except Exception:
            s = str(s)
        if s and s[:1] in ("=", "+", "-", "@", "\t"):
            return "'" + s
        return s
    ordered_keys = [
        'timestamp','operator_id','question_id','filename','patient','visit','matrix','batch_group','pathogen',
        'rpm_species','rpm_genus','ctrl_env_rpm_species','ctrl_env_rpm_genus','ratio_species_genus','label','notes',
        'status','mode','questions_version','app_version'
    ]
    row = [_sanitize_cell(record.get(k, "")) for k in ordered_keys]
    ws.append_row(row, value_input_option="USER_ENTERED")
    return True

def detect_assignment_column(df: pd.DataFrame):
    for c in ['assigned_operator','operator','assigned_to','assignment','assignee','operator_id']:
        if c in df.columns:
            return c
    return None

def assigned_set_for_operator(df: pd.DataFrame, assign_col: str, op_id: str) -> set:
    if not assign_col or assign_col not in df.columns:
        return set()
    ser = df[assign_col].astype(str).fillna("")
    def has_op_cell(v: str) -> bool:
        toks = [t.strip().lower() for t in v.split(',') if t.strip()]
        return op_id.lower() in toks
    mask = ser.apply(has_op_cell)
    return set(df.loc[mask, 'question_id'].astype(str).unique().tolist())

def compute_questions_version(df: pd.DataFrame) -> str:
    try:
        payload = df.to_csv(index=False).encode('utf-8')
        import hashlib
        return hashlib.sha1(payload).hexdigest()[:12]
    except Exception:
        return datetime.utcnow().strftime('%Y%m%d%H%M%S')
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

if operator_id == 'admin':
    responses = load_responses_gsheet() if 'admin_mode' in locals() and admin_mode=='Supervision' else pd.DataFrame()
else:
    responses = load_responses_gsheet() if ('_mode' in locals() and _mode=='Production') else pd.DataFrame()


# ---------------------------
# Header / progress
# ---------------------------

st.title('NGS-IST — Annotation operateur (MVP)')

# === Progress & dashboards ===
questions_version = compute_questions_version(cases)
assign_col = detect_assignment_column(cases)

if operator_id == 'admin':
    if 'admin_mode' in locals() and admin_mode == 'Supervision':
        ops = ['op1','op2','op3','op4','op5']
        def assigned_for(op):
            aset = assigned_set_for_operator(cases, assign_col, op)
            return len(aset)
        assigned_map = {op: assigned_for(op) for op in ops}
        def answered_for(op):
            if responses.empty: return 0
            df = responses.copy()
            df = df[df.get('mode','Production')=='Production']
            df = df[df.get('operator_id','').astype(str).str.lower()==op.lower()]
            if 'timestamp' in df.columns:
                df = df.sort_values('timestamp').drop_duplicates(subset=['question_id'], keep='last')
            df = df[df.get('status','answered')=='answered']
            return int(df['question_id'].astype(str).nunique()) if 'question_id' in df.columns else 0
        answered_map = {op: answered_for(op) for op in ops}
        total_assigned_global = int(sum(assigned_map.values()))
        total_answered_global = int(sum(answered_map.values()))
        assigned_unique = set()
        for op in ops:
            assigned_unique |= set() if assign_col is None else assigned_set_for_operator(cases, assign_col, op)
        covered_unique = set()
        if not responses.empty and 'question_id' in responses.columns:
            dfcov = responses.copy()
            dfcov = dfcov[dfcov.get('mode','Production')=='Production']
            if 'timestamp' in dfcov.columns:
                dfcov = dfcov.sort_values('timestamp').drop_duplicates(subset=['question_id'], keep='last')
            dfcov = dfcov[dfcov.get('status','answered')=='answered']
            covered_unique = set(dfcov['question_id'].astype(str).tolist())
        st.subheader('Supervision (Admin)')
        cA, cB = st.columns(2)
        cA.metric('Couverture unique (≥1 réponse)', f"{len(covered_unique)}/{len(assigned_unique) if assigned_unique else 0}")
        cB.metric('Charge complétée (avec redondance)', f"{total_answered_global}/{total_assigned_global}")
        import pandas as _pd
        df_ops = _pd.DataFrame({'op': ops,
                                'assignés': [assigned_map[o] for o in ops],
                                'répondus': [answered_map[o] for o in ops]})
        df_ops['%'] = df_ops.apply(lambda r: f"{(r['répondus']/r['assignés']*100):.0f}%" if r['assignés'] else '0%', axis=1)
        st.dataframe(df_ops, use_container_width=True)
        st.stop()
else:
    assigned_set = assigned_set_for_operator(cases, assign_col, operator_id)
    total_for_op = int(len(assigned_set))
    if ('_mode' in locals() and _mode=='Production') and not responses.empty:
        df_op = responses.copy()
        df_op = df_op[df_op.get('mode','Production')=='Production']
        df_op = df_op[df_op.get('operator_id','').astype(str).str.lower()==operator_id.lower()]
        if 'timestamp' in df_op.columns:
            df_op = df_op.sort_values('timestamp').drop_duplicates(subset=['question_id'], keep='last')
        df_op = df_op[df_op.get('status','answered')=='answered']
        answered_for_op = int(df_op['question_id'].astype(str).nunique()) if 'question_id' in df_op.columns else 0
    else:
        answered_for_op = 0
    remaining_for_op = max(0, total_for_op - answered_for_op)
    pc1, pc2, pc3 = st.columns(3)
    pc1.metric('Total pour vous', f'{total_for_op}')
    pc2.metric('Vos réponses', f'{answered_for_op}/{total_for_op}')
    pc3.metric('Restantes', f'{remaining_for_op}')

if not operator_id:
    st.info("Entrez votre **identifiant operateur** dans la barre laterale pour commencer.")
    st.stop()

if operator_id!='admin' and ('_mode' in locals() and _mode=='Production') and not responses.empty:
    _df = responses.copy()
    _df = _df[_df.get('mode','Production')=='Production']
    _df = _df[_df.get('operator_id','').astype(str).str.lower()==operator_id.lower()]
    if 'timestamp' in _df.columns:
        _df = _df.sort_values('timestamp').drop_duplicates(subset=['question_id'], keep='last')
    _df = _df[_df.get('status','answered')=='answered']
    already = set(_df['question_id'].astype(str).tolist())
else:
    already = set()
done_by_op = 0['question_id'].nunique()
remaining_global = total - done

pc1, pc2, pc3 = st.columns(3)
pc1.metric('Progression globale', f'{done}/{total}')
pc2.metric('Vos reponses', f'{done_by_op}/{total}')
pc3.metric('Restantes (global)', f'{remaining_global}')

# ---------------------------
# Current question (locked)
# ---------------------------

cases_selection = cases if operator_id=='admin' else cases[cases['question_id'].astype(str).isin(list(assigned_set))]
row = get_current_row(cases_selection, already)
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
    # green, same size as subheader
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

    if manifest_stats is None:
        st.warning("manifest_batch_stats indisponible — fallback sur la méthode v4.")
        # ==== Fallback to original code (kept minimal): compute from cases/manifest if present ====
        batch_mask = (
            (cases['batch_group'].astype(str) == str(row['batch_group'])) &
            (cases['matrix'].astype(str) == str(row['matrix']))
        )
        target_mask = batch_mask & (cases['pathogen'].astype(str) == str(row['pathogen']))
        sub = cases.loc[target_mask, ['sample_id', used_col]].copy()
        sub['sample_id'] = sub['sample_id'].astype(str)
        rpm_map = dict(zip(sub['sample_id'], sub[used_col].astype(float)))
        batch_samples = cases.loc[batch_mask, 'sample_id'].astype(str).unique()
        batch_vals = np.array([rpm_map.get(sid, 0.0) for sid in batch_samples], dtype=float)
        n_batch = len(batch_samples)
        _arr = np.asarray(batch_vals, dtype=float)
        n_pos = int((_arr > 0).sum())
        title_sfx = f"({row['matrix']}, {row['pathogen']}, N={n_pos}/{n_batch})"
        plot_batch_hist(batch_vals, float(row[used_col]), species_mode=species_mode, title_suffix=title_sfx, fs=9)
    else:
        subset = manifest_stats[
            (manifest_stats["batch_group"] == str(row["batch_group"])) &
            (manifest_stats["matrix"] == str(row["matrix"])) &
            (manifest_stats["pathogen"] == str(row["pathogen"]))
        ]
        if subset.empty:
            st.warning("Pas de données pré-calculées dans manifest_batch_stats pour ce batch — fallback méthode v4.")
            # Fallback minimal comme ci-dessus
            batch_mask = (
                (cases['batch_group'].astype(str) == str(row['batch_group'])) &
                (cases['matrix'].astype(str) == str(row['matrix']))
            )
            target_mask = batch_mask & (cases['pathogen'].astype(str) == str(row['pathogen']))
            sub = cases.loc[target_mask, ['sample_id', used_col]].copy()
            sub['sample_id'] = sub['sample_id'].astype(str)
            rpm_map = dict(zip(sub['sample_id'], sub[used_col].astype(float)))
            batch_samples = cases.loc[batch_mask, 'sample_id'].astype(str).unique()
            batch_vals = np.array([rpm_map.get(sid, 0.0) for sid in batch_samples], dtype=float)
            n_batch = len(batch_samples)
            _arr = np.asarray(batch_vals, dtype=float)
            n_pos = int((_arr > 0).sum())
            title_sfx = f"({row['matrix']}, {row['pathogen']}, N={n_pos}/{n_batch})"
            plot_batch_hist(batch_vals, float(row[used_col]), species_mode=species_mode, title_suffix=title_sfx, fs=9)
        else:
            stats_row = subset.iloc[0]
            batch_vals = expand_counts_to_values(stats_row)
            n_batch = int(stats_row.get("n_total", len(batch_vals)))
            try:
                n_pos = int(stats_row.get('n_positive', 0))
                if n_pos == 0:
                    n_pos = (
                        int(stats_row.get('count_0_1', 0))
                        + int(stats_row.get('count_1_10', 0))
                        + int(stats_row.get('count_10_100', 0))
                        + int(stats_row.get('count_100_1000', 0))
                        + int(stats_row.get('count_ge_1000', 0))
                    )
            except Exception:
                n_pos = (
                    int(stats_row.get('count_0_1', 0))
                    + int(stats_row.get('count_1_10', 0))
                    + int(stats_row.get('count_10_100', 0))
                    + int(stats_row.get('count_100_1000', 0))
                    + int(stats_row.get('count_ge_1000', 0))
                )
            title_sfx = f"({row['matrix']}, {row['pathogen']}, N={n_pos}/{n_batch})"
            plot_batch_hist(batch_vals, float(row[used_col]), species_mode=species_mode, title_suffix=title_sfx, fs=9)

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
            'notes': notes,
            'status': 'answered',
            'mode': (_mode if operator_id!='admin' else 'Production'),
            'questions_version': questions_version,
            'app_version': 'v5.5'
        }
        if operator_id!='admin' and ('_mode' in locals() and _mode=='Production'):
            try:
                append_response_gsheet(record)
                st.success("✅ Enregistré (Production) dans la Google Sheet.")
            except Exception as e:
                st.error(f"Erreur d'écriture Google Sheet : {e}")
                st.stop()
        else:
            st.info("🧪 Mode Test ou Admin : aucune écriture.")
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

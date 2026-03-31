from __future__ import annotations

import statistics
import time
from datetime import datetime
from typing import List

import pandas as pd
import streamlit as st

from src.alignment import SequenceRecord, build_multiple_alignment, format_fasta, parse_fasta
from src.ncbi_client import fetch_fasta_by_ids, parse_accession_input, search_nucleotide
from src.phylo import (
    average_similarity,
    build_neighbor_joining_tree,
    distance_matrix_from_similarity,
    newick_from_tree,
    similarity_matrix,
)
from src.sample_projects import SAMPLE_PROJECTS
from src.visualization import render_heatmap_html, render_tree_svg


st.set_page_config(
    page_title="Genome Atlas",
    layout="wide",
    initial_sidebar_state="collapsed",
)


DEMO_ACCESSIONS = SAMPLE_PROJECTS[0].accessions_text
DEMO_QUERY = SAMPLE_PROJECTS[1].search_query
DEMO_FASTA = SAMPLE_PROJECTS[2].fasta_text


CUSTOM_CSS = """
<style>
    :root {
        --ink: #ecf7ff;
        --muted: #98a9bd;
        --panel: rgba(11, 19, 33, 0.82);
        --panel-strong: #101a2f;
        --line: rgba(125, 163, 201, 0.16);
        --accent: #29d3b2;
        --accent-2: #ff8a3d;
        --accent-soft: rgba(41, 211, 178, 0.16);
        --warm: rgba(255, 138, 61, 0.10);
    }
    .stApp {
        background:
            radial-gradient(circle at 0% 0%, rgba(255, 138, 61, 0.14), transparent 26%),
            radial-gradient(circle at 100% 0%, rgba(41, 211, 178, 0.12), transparent 30%),
            radial-gradient(circle at 50% 100%, rgba(30, 87, 153, 0.18), transparent 35%),
            linear-gradient(180deg, #050914 0%, #091121 45%, #0c1830 100%);
        color: var(--ink);
        overflow-x: hidden;
    }
    header[data-testid="stHeader"] {
        background: linear-gradient(90deg, rgba(8,15,28,0.92) 0%, rgba(10,22,40,0.92) 100%) !important;
        border-bottom: 1px solid rgba(125, 163, 201, 0.12);
        backdrop-filter: blur(10px);
    }
    [data-testid="stToolbar"] {
        display: none !important;
    }
    #MainMenu {
        visibility: hidden;
    }
    [data-testid="stToolbar"] {
        right: 1rem;
        top: 0.55rem;
    }
    [data-testid="stToolbar"] button,
    [data-testid="stToolbar"] a,
    [data-testid="stDecoration"] {
        color: #73f0d6 !important;
    }
    [data-testid="stStatusWidget"] {
        color: #d9ecff !important;
    }
    .block-container {
        max-width: 1300px;
        padding-top: 2.3rem;
        padding-bottom: 2.5rem;
        overflow-x: clip;
    }
    h1, h2, h3 {
        letter-spacing: -0.02em;
        color: #f2f8ff;
    }
    .shell {
        background: var(--panel);
        border: 1px solid var(--line);
        box-shadow: 0 24px 60px rgba(0, 0, 0, 0.26);
        border-radius: 28px;
        backdrop-filter: blur(8px);
        padding: 1.25rem;
        margin-bottom: 1rem;
    }
    .hero {
        display: grid;
        grid-template-columns: 1.25fr 0.85fr;
        gap: 1.2rem;
        align-items: stretch;
        width: 100%;
        overflow: hidden;
    }
    .hero-panel {
        background:
            radial-gradient(circle at top right, rgba(255,138,61,0.16), transparent 28%),
            linear-gradient(145deg, #0f172a 0%, #112744 46%, #0e7490 100%);
        color: white;
        border-radius: 26px;
        padding: 1.8rem;
        min-height: 230px;
        box-shadow: 0 24px 60px rgba(0, 0, 0, 0.32);
        min-width: 0;
    }
    .hero-panel h1 {
        color: white;
        font-size: 3rem;
        line-height: 0.96;
        margin: 0 0 0.9rem 0;
    }
    .hero-panel p {
        margin: 0;
        font-size: 1.05rem;
        line-height: 1.65;
        max-width: 52rem;
        opacity: 0.96;
    }
    .hero-points {
        margin: 0;
        padding-left: 1.2rem;
        font-size: 1.02rem;
        line-height: 1.8;
        opacity: 0.96;
    }
    .hero-points li {
        margin-bottom: 0.2rem;
    }
    .hero-note {
        display: inline-block;
        margin-bottom: 0.8rem;
        padding: 0.35rem 0.7rem;
        border-radius: 999px;
        background: rgba(115, 240, 214, 0.10);
        border: 1px solid rgba(115, 240, 214, 0.18);
        font-size: 0.85rem;
        font-weight: 700;
        letter-spacing: 0.04em;
        text-transform: uppercase;
    }
    .hero-side {
        background: linear-gradient(180deg, rgba(13,22,38,0.94) 0%, rgba(11,19,33,0.98) 100%);
        border-radius: 26px;
        padding: 1.3rem;
        border: 1px solid var(--line);
        box-shadow: inset 0 1px 0 rgba(255,255,255,0.04);
        min-width: 0;
    }
    .mini-title {
        font-size: 0.82rem;
        font-weight: 800;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #ffb073;
        margin-bottom: 0.55rem;
    }
    .hero-list {
        margin: 0;
        padding-left: 1.15rem;
        color: #d7e7f8;
        line-height: 1.65;
    }
    .feature-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 0.9rem;
        margin-top: 0.9rem;
        width: 100%;
    }
    .feature-card {
        background: linear-gradient(180deg, rgba(13,22,38,0.86) 0%, rgba(9,17,33,0.92) 100%);
        border-radius: 22px;
        padding: 1rem 1.05rem;
        border: 1px solid var(--line);
        box-shadow: 0 12px 28px rgba(0, 0, 0, 0.22);
    }
    .feature-card h3 {
        margin: 0 0 0.35rem 0;
        font-size: 1.06rem;
    }
    .feature-card p {
        margin: 0;
        color: #a8bad0;
        line-height: 1.6;
        font-size: 0.96rem;
    }
    .panel-title {
        font-size: 1.22rem;
        font-weight: 800;
        color: var(--ink);
        margin-bottom: 0.25rem;
    }
    .panel-copy {
        color: var(--muted);
        line-height: 1.65;
        margin-bottom: 1rem;
    }
    .callout {
        border-left: 4px solid #ff8a3d;
        background: rgba(255, 138, 61, 0.10);
        color: #ffd6bb;
        padding: 0.9rem 1rem;
        border-radius: 14px;
        margin: 0.9rem 0 0.35rem 0;
    }
    .result-chip {
        display: inline-block;
        padding: 0.45rem 0.75rem;
        border-radius: 999px;
        background: linear-gradient(180deg, rgba(41,211,178,0.14) 0%, rgba(30,87,153,0.18) 100%);
        border: 1px solid rgba(115, 240, 214, 0.20);
        color: #9cf6e2;
        font-weight: 700;
        margin-right: 0.45rem;
        margin-bottom: 0.45rem;
    }
    [data-testid="stMetric"] {
        background: linear-gradient(180deg, rgba(14,24,42,0.96) 0%, rgba(10,18,32,1) 100%);
        border-radius: 18px;
        border: 1px solid rgba(125, 163, 201, 0.12);
        padding: 0.95rem 1rem;
        box-shadow: 0 14px 32px rgba(0, 0, 0, 0.22);
    }
    .metric-caption {
        color: #8fa4bc;
        font-size: 0.92rem;
        margin-top: -0.25rem;
    }
    .panel-subcopy {
        color: var(--muted);
        line-height: 1.6;
        margin-bottom: 0.9rem;
    }
    .project-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 1rem;
        margin: 0.9rem 0 1rem 0;
        width: 100%;
    }
    .project-card {
        background: linear-gradient(180deg, rgba(16,27,48,0.98) 0%, rgba(10,18,32,0.98) 100%);
        border: 1px solid rgba(115, 240, 214, 0.12);
        border-radius: 24px;
        padding: 1rem;
        box-shadow: 0 14px 34px rgba(0, 0, 0, 0.22);
        min-height: 190px;
    }
    .project-tag {
        display: inline-block;
        border-radius: 999px;
        background: rgba(255, 138, 61, 0.14);
        color: #ffb073;
        padding: 0.35rem 0.65rem;
        font-size: 0.78rem;
        font-weight: 800;
        letter-spacing: 0.04em;
        text-transform: uppercase;
        margin-bottom: 0.7rem;
    }
    .project-card h4 {
        margin: 0 0 0.35rem 0;
        color: var(--ink);
        font-size: 1.15rem;
    }
    .project-card p {
        margin: 0;
        color: var(--muted);
        line-height: 1.6;
        font-size: 0.95rem;
    }
    .download-strip {
        background: linear-gradient(135deg, rgba(15,24,43,0.92) 0%, rgba(10,20,37,0.96) 100%);
        border: 1px solid rgba(125, 163, 201, 0.12);
        border-radius: 22px;
        padding: 1rem;
        margin-bottom: 1rem;
    }
    .stButton > button,
    .stDownloadButton > button,
    div[data-testid="stFormSubmitButton"] > button {
        background: linear-gradient(135deg, #0ea5a4 0%, #1d4ed8 100%) !important;
        color: white !important;
        border: 0 !important;
        border-radius: 16px !important;
        min-height: 3rem !important;
        font-weight: 700 !important;
        box-shadow: 0 14px 32px rgba(14, 165, 164, 0.22) !important;
        transition: transform 0.18s ease, box-shadow 0.18s ease, filter 0.18s ease !important;
    }
    .stButton > button:hover,
    .stDownloadButton > button:hover,
    div[data-testid="stFormSubmitButton"] > button:hover {
        transform: translateY(-1px);
        filter: brightness(1.03);
        box-shadow: 0 18px 34px rgba(14, 165, 164, 0.26) !important;
    }
    .stButton > button:focus,
    .stDownloadButton > button:focus,
    div[data-testid="stFormSubmitButton"] > button:focus {
        box-shadow: 0 0 0 0.2rem rgba(34, 211, 238, 0.22) !important;
    }
    .stTextInput > div > div > input,
    .stTextArea textarea,
    div[data-baseweb="select"] > div,
    .stNumberInput input {
        background: rgba(9, 17, 30, 0.96) !important;
        color: var(--ink) !important;
        border: 1px solid rgba(115, 240, 214, 0.16) !important;
        border-radius: 18px !important;
        box-shadow: inset 0 1px 0 rgba(255,255,255,0.02), 0 8px 20px rgba(0, 0, 0, 0.16) !important;
    }
    .stTextInput > div > div > input::placeholder,
    .stTextArea textarea::placeholder {
        color: #6f859d !important;
    }
    .stTextInput > label,
    .stTextArea > label,
    .stSlider > label,
    .stRadio > label {
        color: var(--ink) !important;
        font-weight: 700 !important;
    }
    .stTextInput small,
    .stTextArea small {
        color: var(--muted) !important;
    }
    div[data-testid="stRadio"] > div {
        gap: 0.8rem;
    }
    div[data-testid="stRadio"] label {
        background: rgba(12, 21, 38, 0.88);
        border: 1px solid rgba(115, 240, 214, 0.14);
        border-radius: 999px;
        padding: 0.5rem 0.9rem;
        box-shadow: 0 8px 18px rgba(0, 0, 0, 0.14);
    }
    div[data-testid="stRadio"] label p {
        color: var(--ink) !important;
        font-weight: 600 !important;
    }
    div[data-testid="stRadio"] label:has(input:checked) {
        background: linear-gradient(135deg, rgba(14,165,164,0.18) 0%, rgba(29,78,216,0.18) 100%);
        border-color: rgba(115, 240, 214, 0.35);
    }
    div[data-testid="stRadio"] input[type="radio"] {
        accent-color: #29d3b2;
    }
    .stSlider [data-baseweb="slider"] div[role="slider"] {
        background: #ff8a3d !important;
        border-color: #ff8a3d !important;
        box-shadow: 0 0 0 6px rgba(255, 138, 61, 0.16);
    }
    .stSlider [data-baseweb="slider"] > div > div {
        background: linear-gradient(90deg, #29d3b2 0%, #1d4ed8 60%, #ff8a3d 100%) !important;
    }
    .stAlert {
        border-radius: 18px !important;
        border: 1px solid var(--line) !important;
    }
    [data-testid="stDataFrame"] {
        border-radius: 18px;
        overflow: hidden;
        border: 1px solid var(--line);
    }
    .stCodeBlock, pre {
        border-radius: 18px !important;
    }
    .stMarkdown, .stText, label, p, li, div {
        color: inherit;
    }
    [data-testid="stMetricLabel"], [data-testid="stMetricValue"] {
        color: #f5fbff !important;
    }
    [data-testid="stTabs"] button {
        color: #b8cce4 !important;
    }
    [data-testid="stTabs"] button[aria-selected="true"] {
        color: #f4fbff !important;
    }
    @media (max-width: 1200px) {
        .hero {
            grid-template-columns: 1fr;
        }
    }
    @media (max-width: 1100px) {
        .hero, .project-grid {
            grid-template-columns: 1fr;
        }
    }
    @media (max-width: 768px) {
        .block-container {
            padding-top: 1.4rem;
            padding-left: 0.75rem;
            padding-right: 0.75rem;
        }
        .shell, .hero-panel, .hero-side, .feature-card, .project-card {
            border-radius: 20px;
        }
        .hero-panel {
            padding: 1.25rem;
        }
        .hero-panel h1 {
            font-size: 2.25rem;
        }
        .hero-points {
            font-size: 0.96rem;
            line-height: 1.65;
        }
    }
</style>
"""


def init_state() -> None:
    defaults = {
        "input_mode": "NCBI Search",
        "accessions_text": "",
        "search_query": DEMO_QUERY,
        "search_limit": 5,
        "fasta_text": "",
        "email": "",
        "api_key": "",
        "analysis_error": "",
        "records": None,
        "results": None,
        "run_history": [],
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def load_sample_project(project_key: str, auto_run: bool = True) -> None:
    project = next(item for item in SAMPLE_PROJECTS if item.key == project_key)
    st.session_state.input_mode = project.input_mode
    st.session_state.search_query = project.search_query
    st.session_state.accessions_text = project.accessions_text
    st.session_state.fasta_text = project.fasta_text
    st.session_state.search_limit = project.search_limit
    if auto_run:
        run_analysis()


@st.cache_data(show_spinner=False)
def load_records(
    input_mode: str,
    accessions_text: str,
    search_query: str,
    search_limit: int,
    fasta_text: str,
    email: str,
    api_key: str,
) -> List[SequenceRecord]:
    if input_mode == "Accession IDs":
        accessions = parse_accession_input(accessions_text)
        if not accessions:
            raise ValueError("Paste at least 2 accession IDs, one per line or comma-separated.")
        return fetch_fasta_by_ids(accessions, email=email, api_key=api_key)

    if input_mode == "NCBI Search":
        if not search_query.strip():
            raise ValueError("Enter an NCBI nucleotide search query.")
        ids = search_nucleotide(
            term=search_query.strip(),
            retmax=search_limit,
            email=email,
            api_key=api_key,
        )
        return fetch_fasta_by_ids(ids, email=email, api_key=api_key)

    if not fasta_text.strip():
        raise ValueError("Paste at least 2 FASTA records to analyze.")
    records = parse_fasta(fasta_text)
    if len(records) < 2:
        raise ValueError("At least 2 FASTA records are required.")
    return records


@st.cache_data(show_spinner=False)
def run_pipeline(records: List[SequenceRecord]) -> dict:
    alignment = build_multiple_alignment(records)
    sim = similarity_matrix(alignment)
    dist = distance_matrix_from_similarity(sim)
    tree = build_neighbor_joining_tree(alignment, dist)
    return {
        "alignment": alignment,
        "similarity": sim,
        "distance": dist,
        "tree": tree,
        "newick": newick_from_tree(tree) + ";",
    }


def matrix_frame(matrix: List[List[float]], labels: List[str]) -> pd.DataFrame:
    return pd.DataFrame(matrix, index=labels, columns=labels)


def sequence_summary_frame(records: List[SequenceRecord]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Label": record.label,
                "Identifier": record.identifier,
                "Length": len(record.sequence),
                "Description": record.description,
            }
            for record in records
        ]
    )


def similarity_summary_text(labels: List[str], similarity: List[List[float]]) -> str:
    comparisons = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            comparisons.append((similarity[i][j], labels[i], labels[j]))

    if not comparisons:
        return "Only one sequence is available, so there are no pairwise similarity comparisons to summarize."

    highest = max(comparisons, key=lambda item: item[0])
    lowest = min(comparisons, key=lambda item: item[0])
    average_value = average_similarity(similarity) * 100
    return (
        f"The closest pair is {highest[1]} and {highest[2]} at {highest[0] * 100:.2f}% identity, "
        f"while the most distant pair is {lowest[1]} and {lowest[2]} at {lowest[0] * 100:.2f}% identity. "
        f"The overall average pairwise identity across this dataset is {average_value:.2f}%."
    )


def metric_items(records: List[SequenceRecord], alignment: List[SequenceRecord], similarity: List[List[float]]) -> List[tuple[str, str, str]]:
    lengths = [len(record.sequence) for record in records]
    return [
        ("Sequences", str(len(records)), "input records"),
        ("Mean Length", f"{statistics.mean(lengths):.1f} bp", "raw DNA length"),
        ("Alignment", f"{len(alignment[0].sequence)} bp", "aligned columns"),
        ("Identity", f"{average_similarity(similarity) * 100:.2f}%", "average pairwise similarity"),
    ]


def render_project_picker() -> None:
    st.markdown('<div class="panel-title">One-Click Example Projects</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="panel-subcopy">Open a ready-made dataset instantly, or skip this and enter your own data below.</div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="project-grid">'
        + "".join(
            (
                '<div class="project-card">'
                f'<div class="project-tag">{project.subtitle}</div>'
                f'<h4>{project.title}</h4>'
                f'<p>{project.description}</p>'
                '</div>'
            )
            for project in SAMPLE_PROJECTS
        )
        + '</div>',
        unsafe_allow_html=True,
    )

    for col, project in zip(st.columns(len(SAMPLE_PROJECTS)), SAMPLE_PROJECTS):
        with col:
            if st.button(f"Open {project.title}", key=f"project_{project.key}", use_container_width=True):
                load_sample_project(project.key, auto_run=True)


def render_input_mode_fields() -> None:
    st.radio(
        "Input mode",
        ["NCBI Search", "Accession IDs", "Paste FASTA"],
        horizontal=True,
        key="input_mode",
    )

    input_mode = st.session_state.input_mode
    if input_mode == "NCBI Search":
        left, right = st.columns([2.3, 1])
        with left:
            st.text_input(
                "NCBI nucleotide query",
                key="search_query",
                help="Example: COX1[Gene] AND mitochondrion[filter] AND primates[Organism] AND 500:3000[SLEN]",
            )
        with right:
            st.slider("Result count", 3, 15, key="search_limit")
        return

    if input_mode == "Accession IDs":
        st.text_area(
            "Accession IDs",
            key="accessions_text",
            height=170,
            placeholder="Example:\nPV276824.1\nPQ495639.1\nPQ495638.1",
            help="Paste one ID per line or use commas.",
        )
        return

    st.text_area(
        "FASTA input",
        key="fasta_text",
        height=220,
        placeholder=">Seq1\nATGCTAGCTAGC\n>Seq2\nATGCTAGATAGC",
        help="Paste DNA FASTA records directly.",
    )


def render_ncbi_credentials() -> None:
    left, right = st.columns(2)
    with left:
        st.text_input("Email for NCBI requests", key="email")
    with right:
        st.text_input("NCBI API key", key="api_key", type="password")


def render_metrics(records: List[SequenceRecord], alignment: List[SequenceRecord], similarity: List[List[float]]) -> None:
    for column, (label, value, caption) in zip(
        st.columns(4),
        metric_items(records, alignment, similarity),
    ):
        with column:
            st.metric(label, value)
            st.markdown(f'<div class="metric-caption">{caption}</div>', unsafe_allow_html=True)


def run_analysis() -> None:
    start_time = time.perf_counter()
    try:
        with st.spinner("Collecting DNA sequences and building the tree..."):
            records = load_records(
                st.session_state.input_mode,
                st.session_state.accessions_text,
                st.session_state.search_query,
                st.session_state.search_limit,
                st.session_state.fasta_text,
                st.session_state.email,
                st.session_state.api_key,
            )
            if len(records) < 2:
                raise ValueError("At least 2 sequences are required to continue.")
            results = run_pipeline(records)
    except Exception as exc:
        st.session_state.analysis_error = str(exc)
        st.session_state.records = None
        st.session_state.results = None
        return

    st.session_state.analysis_error = ""
    st.session_state.records = records
    st.session_state.results = results
    elapsed_seconds = round(time.perf_counter() - start_time, 3)
    alignment = results["alignment"]
    similarity = results["similarity"]
    lengths = [len(record.sequence) for record in records]
    history_item = {
        "run_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "input_mode": st.session_state.input_mode,
        "sequence_count": len(records),
        "mean_length": round(statistics.mean(lengths), 2),
        "alignment_length": len(alignment[0].sequence),
        "avg_identity_pct": round(average_similarity(similarity) * 100, 3),
        "runtime_seconds": elapsed_seconds,
    }
    st.session_state.run_history = (st.session_state.run_history + [history_item])[-30:]


def render_header() -> None:
    st.markdown(CUSTOM_CSS, unsafe_allow_html=True)
    st.markdown(
        """
        <div class="hero">
            <div class="hero-panel">
                <div class="hero-note">Comparative Genomics Workspace</div>
                <h1>Genome Atlas</h1>
                <ul class="hero-points">
                    <li>Fetch DNA sequences from NCBI or paste FASTA directly</li>
                    <li>Align sequences and calculate pairwise similarity</li>
                    <li>Build a distance matrix for comparative analysis</li>
                    <li>Generate a Neighbor Joining phylogenetic tree</li>
                </ul>
            </div>
            <div class="hero-side">
                <div class="mini-title">What This App Does</div>
                <ul class="hero-list">
                    <li>Fetch nucleotide FASTA records directly from NCBI</li>
                    <li>Accept user-pasted FASTA for offline testing</li>
                    <li>Build a lightweight multiple alignment</li>
                    <li>Compute pairwise sequence similarity and distance</li>
                    <li>Render a phylogenetic tree and export the results</li>
                </ul>
                <div class="callout">
                    Use the demo loaders below if you want a guaranteed working example before trying your own data.
                </div>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_input_panel() -> None:
    st.markdown('<div class="shell">', unsafe_allow_html=True)
    st.markdown('<div class="panel-title">Input Studio</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="panel-copy">Choose how you want to provide DNA sequences, then run the analysis.</div>',
        unsafe_allow_html=True,
    )
    render_project_picker()
    render_input_mode_fields()
    render_ncbi_credentials()

    submitted = st.button("Analyze DNA Sequences", type="primary", use_container_width=True)

    if submitted:
        run_analysis()

    if st.session_state.analysis_error:
        st.error(st.session_state.analysis_error)

    st.markdown("</div>", unsafe_allow_html=True)


def render_results() -> None:
    records = st.session_state.records
    results = st.session_state.results

    if not records or not results:
        st.info("Run an analysis to see sequences, alignment, similarity, distance, and phylogenetic tree output.")
        return

    alignment = results["alignment"]
    similarity = results["similarity"]
    distance = results["distance"]
    tree = results["tree"]
    newick = results["newick"]
    labels = [record.label for record in records]
    fasta_text = format_fasta(records)
    alignment_text = format_fasta(alignment)
    summary_text = similarity_summary_text(labels, similarity)

    render_metrics(records, alignment, similarity)

    overview_tab, sequence_tab, matrix_tab, tree_tab, monitor_tab = st.tabs(
        ["Overview", "Sequences", "Matrices", "Tree", "Monitoring"]
    )

    with overview_tab:
        st.markdown('<div class="shell">', unsafe_allow_html=True)
        st.markdown('<div class="panel-title">Sequence Summary</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="panel-subcopy">A line-up view of the fetched DNA records so the identifiers and descriptions stay readable.</div>',
            unsafe_allow_html=True,
        )
        st.dataframe(sequence_summary_frame(records), use_container_width=True, hide_index=True)
        st.markdown("</div>", unsafe_allow_html=True)

        st.markdown('<div class="shell">', unsafe_allow_html=True)
        st.markdown('<div class="panel-title">Tree Preview</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="panel-subcopy">The tree is shown in a taller, full-width panel with branch-length hints and clearer tip labels.</div>',
            unsafe_allow_html=True,
        )
        st.markdown(
            render_tree_svg(tree, width=1160, height=max(420, 170 + 96 * len(records))),
            unsafe_allow_html=True,
        )
        st.markdown("</div>", unsafe_allow_html=True)

        st.markdown('<div class="shell">', unsafe_allow_html=True)
        st.markdown('<div class="panel-title">Similarity Snapshot</div>', unsafe_allow_html=True)
        st.markdown(
            f'<div class="panel-subcopy">{summary_text}</div>',
            unsafe_allow_html=True,
        )
        st.markdown(render_heatmap_html(similarity, labels), unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)

    with sequence_tab:
        left, right = st.columns(2)
        with left:
            st.markdown('<div class="shell">', unsafe_allow_html=True)
            st.markdown('<div class="panel-title">DNA FASTA</div>', unsafe_allow_html=True)
            st.code(fasta_text, language="text")
            st.markdown("</div>", unsafe_allow_html=True)
        with right:
            st.markdown('<div class="shell">', unsafe_allow_html=True)
            st.markdown('<div class="panel-title">Multiple Alignment</div>', unsafe_allow_html=True)
            st.code(alignment_text, language="text")
            st.markdown("</div>", unsafe_allow_html=True)

    with matrix_tab:
        sim_df = matrix_frame(similarity, labels)
        dist_df = matrix_frame(distance, labels)
        left, right = st.columns(2)
        with left:
            st.markdown('<div class="shell">', unsafe_allow_html=True)
            st.markdown('<div class="panel-title">Similarity Matrix</div>', unsafe_allow_html=True)
            st.markdown(render_heatmap_html(similarity, labels), unsafe_allow_html=True)
            st.markdown("</div>", unsafe_allow_html=True)
        with right:
            st.markdown('<div class="shell">', unsafe_allow_html=True)
            st.markdown('<div class="panel-title">Distance Matrix</div>', unsafe_allow_html=True)
            st.dataframe(dist_df, use_container_width=True)
            st.markdown("</div>", unsafe_allow_html=True)

    with tree_tab:
        st.markdown('<div class="shell">', unsafe_allow_html=True)
        st.markdown('<div class="panel-title">Neighbor Joining Tree</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="panel-subcopy">Expanded tree canvas with improved branch styling, labels, and branch-length annotations.</div>',
            unsafe_allow_html=True,
        )
        st.markdown(
            render_tree_svg(tree, width=1200, height=max(480, 190 + 100 * len(records))),
            unsafe_allow_html=True,
        )
        st.markdown('<div class="panel-title">Newick Export</div>', unsafe_allow_html=True)
        st.code(newick, language="text")
        st.markdown("</div>", unsafe_allow_html=True)

    with monitor_tab:
        st.markdown('<div class="shell">', unsafe_allow_html=True)
        st.markdown('<div class="panel-title">Run Monitoring</div>', unsafe_allow_html=True)
        st.markdown(
            '<div class="panel-subcopy">Track runtime and analysis characteristics across recent runs in this session.</div>',
            unsafe_allow_html=True,
        )
        history = st.session_state.run_history
        if not history:
            st.info("No monitoring data yet. Run an analysis to generate monitoring graphs.")
        else:
            history_df = pd.DataFrame(history)
            st.dataframe(history_df, use_container_width=True, hide_index=True)

            st.markdown('<div class="panel-title">Runtime Trend (seconds)</div>', unsafe_allow_html=True)
            st.line_chart(
                history_df.set_index("run_time")[["runtime_seconds"]],
                use_container_width=True,
            )

            col1, col2 = st.columns(2)
            with col1:
                st.markdown('<div class="panel-title">Sequence Count Trend</div>', unsafe_allow_html=True)
                st.line_chart(
                    history_df.set_index("run_time")[["sequence_count"]],
                    use_container_width=True,
                )
            with col2:
                st.markdown('<div class="panel-title">Average Identity Trend (%)</div>', unsafe_allow_html=True)
                st.line_chart(
                    history_df.set_index("run_time")[["avg_identity_pct"]],
                    use_container_width=True,
                )

            st.markdown('<div class="panel-title">Alignment Length Trend</div>', unsafe_allow_html=True)
            st.bar_chart(
                history_df.set_index("run_time")[["alignment_length"]],
                use_container_width=True,
            )
        st.markdown("</div>", unsafe_allow_html=True)


def main() -> None:
    init_state()
    render_header()
    render_input_panel()
    render_results()


if __name__ == "__main__":
    main()

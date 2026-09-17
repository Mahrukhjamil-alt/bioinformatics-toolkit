import streamlit as st
from sections import single_seq, psa, msa, phylo, orf_finder, restriction_sites, motif_search

st.set_page_config(page_title="Bioinformatics Toolkit", layout="wide")

st.markdown(
    """
    <style>
    .top-nav {
        position: sticky;
        top: 0;
        z-index: 999;
        background: linear-gradient(90deg, #1f2937, #111827);
        padding: 14px 22px;
        border-radius: 0 0 14px 14px;
        box-shadow: 0 4px 14px rgba(0,0,0,0.35);
        margin-bottom: 30px;
        display: flex;
        flex-wrap: wrap;
        gap: 10px;
        align-items: center;
    }
    .top-nav .brand {
        color: #ffffff;
        font-weight: 800;
        font-size: 17px;
        margin-right: 18px;
        white-space: nowrap;
    }
    .top-nav a {
        color: #e5e7eb;
        text-decoration: none;
        font-weight: 600;
        font-size: 14px;
        padding: 8px 16px;
        border-radius: 999px;
        background-color: rgba(255,255,255,0.06);
        border: 1px solid rgba(255,255,255,0.08);
        transition: all 0.2s ease;
        white-space: nowrap;
    }
    .top-nav a:hover {
        background-color: #ff4b4b;
        border-color: #ff4b4b;
        color: white;
        transform: translateY(-2px);
    }
    </style>
    <div class="top-nav">
        <span class="brand">🧬 BioToolkit</span>
        <a href="#single-seq">Single Sequence</a>
        <a href="#psa">Pairwise Alignment</a>
        <a href="#msa">MSA</a>
        <a href="#phylo">Phylogenetic Tree</a>
        <a href="#orf">ORF Finder</a>
        <a href="#restriction">Restriction Sites</a>
        <a href="#motif">Motif Search</a>
    </div>
    """,
    unsafe_allow_html=True,
)

st.title("Bioinformatics Toolkit")

single_seq.render()
st.markdown("---")
psa.render()
st.markdown("---")
msa.render()
st.markdown("---")
phylo.render()
st.markdown("---")
orf_finder.render()
st.markdown("---")
restriction_sites.render()
st.markdown("---")
motif_search.render()
import streamlit as st
from Bio import Align
from Bio.Align import substitution_matrices
from sample_data import SAMPLE_MSA_FASTA, SAMPLE_MSA_FASTA_DNA
from utils import parse_fasta, load_uploaded_once
from msa_utils import center_star_msa


def render():
    st.header("Multiple Sequence Alignment (MSA)", anchor="msa")
    st.write("Paste 3+ sequences in FASTA format.")

    if "msa_fasta" not in st.session_state:
        st.session_state.msa_fasta = ""

    c1, c2, c3 = st.columns([1, 1, 4])
    with c1:
        if st.button("Load Protein Sample", key="load_msa_sample_protein"):
            st.session_state.msa_fasta = SAMPLE_MSA_FASTA
    with c2:
        if st.button("Load DNA Sample", key="load_msa_sample_dna"):
            st.session_state.msa_fasta = SAMPLE_MSA_FASTA_DNA
    with c3:
        if st.button("Clear", key="clear_msa"):
            st.session_state.msa_fasta = ""

    uploaded = st.file_uploader("Or upload a multi-FASTA file", type=["fasta", "fa", "txt"], key="msa_upload")
    load_uploaded_once(uploaded, "msa_fasta", "_msa_upload_id", multi=True)

    st.text_area("Sequences (FASTA format):", height=200, key="msa_fasta")
    msa_matrix = st.selectbox("Substitution matrix:", ["BLOSUM62", "BLOSUM50", "PAM250", "NUC.4.4"], key="msa_matrix")

    if st.button("Run MSA", key="run_msa", type="primary"):
        seqs = parse_fasta(st.session_state.msa_fasta)

        if len(seqs) < 2:
            st.warning("Please provide at least 2 sequences.")
        else:
            DNA_CHARS = set("ACGTUN")
            all_chars = set("".join(seqs.values()))
            looks_like_dna = all_chars <= DNA_CHARS
            is_dna_matrix = msa_matrix == "NUC.4.4"

            if is_dna_matrix and not looks_like_dna:
                st.warning(
                    f"'{msa_matrix}' is a nucleotide matrix, but your sequences contain "
                    f"protein letters ({sorted(all_chars - DNA_CHARS)}). "
                    "Select BLOSUM62, BLOSUM50, or PAM250 instead."
                )
            elif not is_dna_matrix and looks_like_dna:
                st.warning(
                    f"Your sequences look like DNA, but '{msa_matrix}' is a protein matrix. "
                    "Select NUC.4.4 instead, or confirm this is really a protein sequence."
                )
            else:
                try:
                    aligner = Align.PairwiseAligner()
                    aligner.mode = "global"
                    try:
                        aligner.substitution_matrix = substitution_matrices.load(msa_matrix)
                    except Exception as e:
                        st.warning(f"Could not load matrix '{msa_matrix}' ({e}). Using default scoring.")
                    aligner.open_gap_score = -10
                    aligner.extend_gap_score = -0.5

                    with st.spinner("Aligning..."):
                        aligned = center_star_msa(seqs, aligner)

                    st.divider()
                    st.success(f"Aligned {len(aligned)} sequences.")
                    max_len = max(len(n) for n in aligned)
                    lines = [f"{name.ljust(max_len)}  {seq}" for name, seq in aligned.items()]
                    st.code("\n".join(lines))

                    fasta_output = "\n".join(f">{name}\n{seq}" for name, seq in aligned.items())
                    st.download_button("Download Alignment (FASTA)", fasta_output, file_name="msa_alignment.fasta", mime="text/plain")
                except Exception as e:
                    st.error(f"MSA failed: {e}")
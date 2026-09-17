import streamlit as st
from Bio import Align
from Bio.Align import substitution_matrices
from sample_data import SAMPLE_PAIRS
from utils import clean_sequence, parse_fasta


def render():
    st.header("Pairwise Sequence Alignment (PSA)", anchor="psa")
    st.write("Align two sequences using a substitution matrix.")

    if "psa_seq1" not in st.session_state:
        st.session_state.psa_seq1 = ""
    if "psa_seq2" not in st.session_state:
        st.session_state.psa_seq2 = ""

    psa_matrix = st.selectbox("Substitution matrix:", list(SAMPLE_PAIRS.keys()), key="psa_matrix")

    c1, c2 = st.columns([1, 5])
    with c1:
        label, s1_sample, s2_sample = SAMPLE_PAIRS[psa_matrix]
        if st.button(f"Load Sample ({label})", key="load_psa_sample"):
            st.session_state.psa_seq1 = s1_sample
            st.session_state.psa_seq2 = s2_sample
    with c2:
        if st.button("Clear", key="clear_psa"):
            st.session_state.psa_seq1 = ""
            st.session_state.psa_seq2 = ""

    uploaded = st.file_uploader("Or upload a FASTA file with 2 sequences", type=["fasta", "fa", "txt"], key="psa_upload")
    if uploaded is not None:
        identifier = f"{uploaded.name}-{uploaded.size}"
        if st.session_state.get("_psa_upload_id") != identifier:
            content = uploaded.read().decode("utf-8", errors="ignore")
            values = list(parse_fasta(content).values())
            if len(values) >= 2:
                st.session_state.psa_seq1, st.session_state.psa_seq2 = values[0], values[1]
            elif len(values) == 1:
                st.session_state.psa_seq1 = values[0]
            st.session_state._psa_upload_id = identifier

    col1, col2 = st.columns(2)
    with col1:
        st.text_area("Sequence 1:", height=120, key="psa_seq1")
    with col2:
        st.text_area("Sequence 2:", height=120, key="psa_seq2")

    psa_mode = st.radio("Alignment type:", ["global", "local"], horizontal=True, key="psa_mode")

    if st.button("Run Alignment", key="run_psa", type="primary"):
        s1 = clean_sequence(st.session_state.psa_seq1)
        s2 = clean_sequence(st.session_state.psa_seq2)

        if not s1 or not s2:
            st.warning("Please provide both sequences.")
        else:
            try:
                aligner = Align.PairwiseAligner()
                aligner.mode = psa_mode
                try:
                    aligner.substitution_matrix = substitution_matrices.load(psa_matrix)
                except Exception as e:
                    st.warning(f"Could not load matrix '{psa_matrix}' ({e}). Using default scoring.")
                aligner.open_gap_score = -10
                aligner.extend_gap_score = -0.5

                best = aligner.align(s1, s2)[0]

                st.divider()
                st.success(f"Score: {best.score:.1f}")
                alignment_text = str(best)
                st.code(alignment_text)

                c1, c2, c3 = st.columns(3)
                c1.metric("Seq 1 Length", len(s1))
                c2.metric("Seq 2 Length", len(s2))
                c3.metric("Matrix", psa_matrix)

                report = f"Score: {best.score:.1f}\nMatrix: {psa_matrix}\nMode: {psa_mode}\n\n{alignment_text}"
                st.download_button("Download Alignment (TXT)", report, file_name="psa_alignment.txt", mime="text/plain")
            except Exception as e:
                st.error(f"Alignment failed: {e}")
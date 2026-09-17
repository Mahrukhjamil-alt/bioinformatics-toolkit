import io
import streamlit as st
from Bio import Align, Phylo
from Bio.Align import substitution_matrices, MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sample_data import SAMPLE_MSA_FASTA, SAMPLE_MSA_FASTA_DNA
from utils import parse_fasta, load_uploaded_once
from msa_utils import center_star_msa


def render():
    st.header("Phylogenetic Tree", anchor="phylo")
    st.write("Paste 3+ sequences in FASTA format. They'll be aligned, then a tree built.")

    if "phylo_fasta" not in st.session_state:
        st.session_state.phylo_fasta = ""

    c1, c2, c3 = st.columns([1, 1, 4])
    with c1:
        if st.button("Load Protein Sample", key="load_phylo_sample_protein"):
            st.session_state.phylo_fasta = SAMPLE_MSA_FASTA
    with c2:
        if st.button("Load DNA Sample", key="load_phylo_sample_dna"):
            st.session_state.phylo_fasta = SAMPLE_MSA_FASTA_DNA
    with c3:
        if st.button("Clear", key="clear_phylo"):
            st.session_state.phylo_fasta = ""

    uploaded = st.file_uploader("Or upload a multi-FASTA file", type=["fasta", "fa", "txt"], key="phylo_upload")
    load_uploaded_once(uploaded, "phylo_fasta", "_phylo_upload_id", multi=True)

    st.text_area("Sequences (FASTA format):", height=200, key="phylo_fasta")
    method = st.radio("Tree method:", ["UPGMA", "Neighbor-Joining (NJ)"], horizontal=True, key="phylo_method")

    if st.button("Build Tree", key="run_phylo", type="primary"):
        seqs = parse_fasta(st.session_state.phylo_fasta)

        if len(seqs) < 3:
            st.warning("Please provide at least 3 sequences.")
        else:
            DNA_CHARS = set("ACGTUN")
            all_chars = set("".join(seqs.values()))
            is_dna = all_chars <= DNA_CHARS
            matrix_name = "NUC.4.4" if is_dna else "BLOSUM62"

            try:
                with st.spinner("Aligning..."):
                    aligner = Align.PairwiseAligner()
                    aligner.mode = "global"
                    try:
                        aligner.substitution_matrix = substitution_matrices.load(matrix_name)
                    except Exception:
                        pass
                    aligner.open_gap_score = -10
                    aligner.extend_gap_score = -0.5
                    aligned = center_star_msa(seqs, aligner)

                with st.spinner("Building tree..."):
                    records = [SeqRecord(Seq(seq), id=name) for name, seq in aligned.items()]
                    msa_obj = MultipleSeqAlignment(records)
                    dm = DistanceCalculator("identity").get_distance(msa_obj)
                    constructor = DistanceTreeConstructor()
                    tree = constructor.upgma(dm) if method == "UPGMA" else constructor.nj(dm)

                st.divider()
                st.caption(f"Alignment used {matrix_name} scoring (auto-detected {'DNA' if is_dna else 'protein'}).")
                st.subheader("Tree")
                buf = io.StringIO()
                Phylo.draw_ascii(tree, file=buf)
                ascii_tree = buf.getvalue()
                st.code(ascii_tree)

                newick_buf = io.StringIO()
                Phylo.write(tree, newick_buf, "newick")
                newick_text = newick_buf.getvalue()

                with st.expander("Newick format (for iTOL / FigTree)"):
                    st.code(newick_text)

                dl1, dl2 = st.columns(2)
                with dl1:
                    st.download_button("Download Tree (Newick)", newick_text, file_name="tree.nwk", mime="text/plain")
                with dl2:
                    st.download_button("Download Tree (ASCII TXT)", ascii_tree, file_name="tree_ascii.txt", mime="text/plain")
            except Exception as e:
                st.error(f"Tree construction failed: {e}")
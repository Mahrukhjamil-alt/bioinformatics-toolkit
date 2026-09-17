import streamlit as st
import pandas as pd
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from sample_data import SAMPLE_SINGLE_PROTEIN, SAMPLE_SINGLE_DNA
from utils import clean_sequence, load_uploaded_once

DNA_CHARS = set("ACGTUN")


def render():
    st.header("Single Sequence Analysis", anchor="single-seq")
    st.write("Paste one sequence (DNA or protein). All relevant calculations run together.")

    if "single_seq" not in st.session_state:
        st.session_state.single_seq = ""

    c1, c2, c3 = st.columns([1, 1, 4])
    with c1:
        if st.button("Load Protein Sample", key="load_single_protein"):
            st.session_state.single_seq = SAMPLE_SINGLE_PROTEIN
    with c2:
        if st.button("Load DNA Sample", key="load_single_dna"):
            st.session_state.single_seq = SAMPLE_SINGLE_DNA
    with c3:
        if st.button("Clear", key="clear_single"):
            st.session_state.single_seq = ""

    uploaded = st.file_uploader("Or upload a FASTA/text file", type=["fasta", "fa", "txt"], key="single_upload")
    load_uploaded_once(uploaded, "single_seq", "_single_upload_id")

    st.text_area("Sequence:", height=140, key="single_seq")

    if st.button("Analyze", key="run_single", type="primary"):
        seq = clean_sequence(st.session_state.single_seq)

        if not seq:
            st.warning("Please provide a sequence.")
            return

        is_dna = set(seq) <= DNA_CHARS
        report_lines = [f"Length: {len(seq)}", f"Detected Type: {'DNA/RNA' if is_dna else 'Protein'}"]

        st.divider()
        st.subheader("Basic Info")
        col1, col2 = st.columns(2)
        col1.metric("Length", len(seq))
        col2.metric("Detected Type", "DNA/RNA" if is_dna else "Protein")

        if is_dna:
            st.divider()
            st.subheader("Nucleotide Composition")
            try:
                gc = gc_fraction(seq) * 100
                mw = molecular_weight(seq, seq_type="DNA")
                c1, c2, c3 = st.columns(3)
                c1.metric("GC Content (%)", round(gc, 2))
                c2.metric("AT Content (%)", round(100 - gc, 2))
                c3.metric("Molecular Weight (Da)", round(mw, 2))

                base_counts = {b: seq.count(b) for b in "ACGTU" if seq.count(b)}
                st.write("  ".join(f"{b}: {c}" for b, c in base_counts.items()))
                st.bar_chart(pd.DataFrame({"Count": list(base_counts.values())}, index=list(base_counts.keys())))

                report_lines += [
                    f"GC Content: {round(gc, 2)}%",
                    f"AT Content: {round(100 - gc, 2)}%",
                    f"Molecular Weight: {round(mw, 2)} Da",
                    f"Base Counts: {base_counts}",
                ]
            except Exception as e:
                st.error(f"Could not compute nucleotide stats: {e}")

            st.divider()
            st.subheader("Transcription & Translation")
            try:
                bio_seq = Seq(seq)
                mrna = bio_seq.transcribe()
                protein = bio_seq.translate(to_stop=True)
                st.write("**mRNA (5'→3'):**")
                st.code(str(mrna))
                st.write("**Protein translation (frame 1, stops at first stop codon):**")
                st.code(str(protein))
                report_lines += [f"mRNA: {mrna}", f"Protein (frame 1): {protein}"]
            except Exception as e:
                st.error(f"Could not transcribe/translate: {e}")

            st.divider()
            st.subheader("Reverse Complement")
            try:
                comp = bio_seq.complement()
                rev = bio_seq[::-1]
                revcomp = bio_seq.reverse_complement()
                st.write("**Complement:**")
                st.code(str(comp))
                st.write("**Reverse:**")
                st.code(str(rev))
                st.write("**Reverse Complement (5'→3'):**")
                st.code(str(revcomp))
                report_lines += [f"Complement: {comp}", f"Reverse: {rev}", f"Reverse Complement: {revcomp}"]
            except Exception as e:
                st.error(f"Could not compute reverse complement: {e}")

        else:
            st.divider()
            st.subheader("Physicochemical Properties")
            try:
                analysis = ProteinAnalysis(seq)
                mw, pi = analysis.molecular_weight(), analysis.isoelectric_point()
                arom, instab = analysis.aromaticity(), analysis.instability_index()

                c1, c2 = st.columns(2)
                c1.metric("Molecular Weight (Da)", round(mw, 2))
                c1.metric("Isoelectric Point (pI)", round(pi, 2))
                c2.metric("Aromaticity", round(arom, 3))
                c2.metric("Instability Index", round(instab, 2))

                st.subheader("Amino Acid Composition")
                comp = {aa: pct for aa, pct in analysis.amino_acids_percent.items() if pct > 0}
                st.write("  ".join(f"{aa}: {pct:.1f}%" for aa, pct in sorted(comp.items(), key=lambda x: -x[1])))
                chart_df = pd.DataFrame({"Percent": list(comp.values())}, index=list(comp.keys())).sort_values("Percent", ascending=False)
                st.bar_chart(chart_df)

                report_lines += [
                    f"Molecular Weight: {round(mw, 2)} Da",
                    f"Isoelectric Point: {round(pi, 2)}",
                    f"Aromaticity: {round(arom, 3)}",
                    f"Instability Index: {round(instab, 2)}",
                    f"Amino Acid Composition: {comp}",
                ]
            except Exception as e:
                st.error(f"Could not compute protein stats: {e}")

        st.divider()
        report_text = "Sequence Analysis Report\n" + "=" * 30 + "\n" + "\n".join(report_lines)
        st.download_button("Download Report (TXT)", report_text, file_name="sequence_report.txt", mime="text/plain")
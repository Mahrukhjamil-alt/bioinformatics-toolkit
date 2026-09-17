import streamlit as st
import pandas as pd
from collections import Counter
from Bio.Seq import Seq
from sample_data import SAMPLE_SINGLE_DNA
from utils import clean_sequence, load_uploaded_once


def find_orfs(seq: str, min_length: int = 30):
    orfs = []
    strands = [("+", seq), ("-", str(Seq(seq).reverse_complement()))]
    stops = {"TAA", "TAG", "TGA"}

    for strand_label, s in strands:
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i:i + 3] == "ATG":
                    j = i
                    found_stop = False
                    while j < len(s) - 2:
                        if s[j:j + 3] in stops:
                            found_stop = True
                            break
                        j += 3
                    if found_stop:
                        orf_len = j + 3 - i
                        if orf_len >= min_length:
                            protein = str(Seq(s[i:j + 3]).translate(to_stop=True))
                            orfs.append({
                                "Strand": strand_label, "Frame": frame + 1,
                                "Start": i + 1, "End": j + 3,
                                "Length (nt)": orf_len, "Protein": protein,
                            })
                        i = j + 3
                    else:
                        break
                else:
                    i += 3
    return sorted(orfs, key=lambda x: -x["Length (nt)"])


def codon_usage(seq: str):
    seq = seq[: len(seq) - len(seq) % 3]
    return Counter(seq[i:i + 3] for i in range(0, len(seq), 3))


def render():
    st.header("ORF Finder & Codon Usage", anchor="orf")
    st.write("Find Open Reading Frames (ATG → stop codon) in all 6 reading frames, and view codon usage.")

    if "orf_seq" not in st.session_state:
        st.session_state.orf_seq = ""

    c1, c2 = st.columns([1, 5])
    with c1:
        if st.button("Load DNA Sample", key="load_orf_sample"):
            st.session_state.orf_seq = SAMPLE_SINGLE_DNA
    with c2:
        if st.button("Clear", key="clear_orf"):
            st.session_state.orf_seq = ""

    uploaded = st.file_uploader("Or upload a DNA FASTA/text file", type=["fasta", "fa", "txt"], key="orf_upload")
    load_uploaded_once(uploaded, "orf_seq", "_orf_upload_id")

    st.text_area("DNA Sequence:", height=140, key="orf_seq")
    min_len = st.slider("Minimum ORF length (nucleotides):", 15, 300, 30, step=3)

    if st.button("Find ORFs", key="run_orf", type="primary"):
        seq = clean_sequence(st.session_state.orf_seq)

        if not seq:
            st.warning("Please provide a DNA sequence.")
        elif not set(seq) <= set("ACGTUN"):
            st.warning("This doesn't look like a DNA sequence.")
        else:
            try:
                orfs = find_orfs(seq, min_length=min_len)

                st.divider()
                st.subheader(f"Open Reading Frames Found: {len(orfs)}")
                if orfs:
                    df = pd.DataFrame(orfs)
                    st.dataframe(df, use_container_width=True)

                    top = df.head(10).copy()
                    top["Label"] = top.apply(lambda r: f"{r['Strand']}F{r['Frame']}:{r['Start']}", axis=1)
                    st.subheader("Top ORF Lengths")
                    st.bar_chart(top.set_index("Label")[["Length (nt)"]])

                    st.download_button("Download ORFs (CSV)", df.to_csv(index=False), file_name="orfs.csv", mime="text/csv")
                else:
                    st.info("No ORFs found with the current minimum length. Try lowering it.")

                st.divider()
                st.subheader("Codon Usage (frame 1, forward strand)")
                usage = codon_usage(seq)
                usage_df = pd.DataFrame({"Count": list(usage.values())}, index=list(usage.keys())).sort_values("Count", ascending=False)
                st.bar_chart(usage_df)
            except Exception as e:
                st.error(f"ORF finding failed: {e}")
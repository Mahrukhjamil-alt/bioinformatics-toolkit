import streamlit as st
import pandas as pd
import re
from sample_data import SAMPLE_SINGLE_PROTEIN
from utils import clean_sequence, load_uploaded_once


def render():
    st.header("Motif / Pattern Search", anchor="motif")
    st.write("Search for a substring or regex pattern within a sequence.")

    if "motif_seq" not in st.session_state:
        st.session_state.motif_seq = ""

    c1, c2 = st.columns([1, 5])
    with c1:
        if st.button("Load Sample", key="load_motif_sample"):
            st.session_state.motif_seq = SAMPLE_SINGLE_PROTEIN
    with c2:
        if st.button("Clear", key="clear_motif"):
            st.session_state.motif_seq = ""

    uploaded = st.file_uploader("Or upload a FASTA/text file", type=["fasta", "fa", "txt"], key="motif_upload")
    load_uploaded_once(uploaded, "motif_seq", "_motif_upload_id")

    st.text_area("Sequence:", height=140, key="motif_seq")
    pattern = st.text_input("Pattern to search for:", placeholder="e.g. GHG or C..C (regex)")
    use_regex = st.checkbox("Treat pattern as regex", value=False)

    if st.button("Search", key="run_motif", type="primary"):
        seq = clean_sequence(st.session_state.motif_seq)

        if not seq or not pattern:
            st.warning("Please provide both a sequence and a pattern.")
        else:
            try:
                positions = []
                if use_regex:
                    for m in re.finditer(pattern, seq, re.IGNORECASE):
                        positions.append((m.start() + 1, m.end(), m.group()))
                else:
                    p = pattern.upper()
                    start = 0
                    while True:
                        idx = seq.find(p, start)
                        if idx == -1:
                            break
                        positions.append((idx + 1, idx + len(p), p))
                        start = idx + 1

                st.divider()
                st.metric("Matches Found", len(positions))

                if positions:
                    df = pd.DataFrame(positions, columns=["Start", "End", "Match"])
                    st.dataframe(df, use_container_width=True)
                    st.download_button("Download Matches (CSV)", df.to_csv(index=False), file_name="motif_matches.csv", mime="text/csv")
                else:
                    st.info("No matches found.")
            except re.error as e:
                st.error(f"Invalid regex pattern: {e}")
            except Exception as e:
                st.error(f"Search failed: {e}")
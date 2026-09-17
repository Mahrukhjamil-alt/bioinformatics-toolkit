import streamlit as st
import pandas as pd
from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Seq import Seq
from sample_data import SAMPLE_SINGLE_DNA
from utils import clean_sequence, load_uploaded_once

COMMON_ENZYMES = ["EcoRI", "BamHI", "HindIII", "PstI", "NotI", "XhoI", "SalI", "SmaI", "KpnI", "SacI"]


def render():
    st.header("Restriction Enzyme Sites", anchor="restriction")
    st.write("Find common restriction enzyme cut sites in a DNA sequence.")

    if "restriction_seq" not in st.session_state:
        st.session_state.restriction_seq = ""

    c1, c2 = st.columns([1, 5])
    with c1:
        if st.button("Load DNA Sample", key="load_restriction_sample"):
            st.session_state.restriction_seq = SAMPLE_SINGLE_DNA
    with c2:
        if st.button("Clear", key="clear_restriction"):
            st.session_state.restriction_seq = ""

    uploaded = st.file_uploader("Or upload a DNA FASTA/text file", type=["fasta", "fa", "txt"], key="restriction_upload")
    load_uploaded_once(uploaded, "restriction_seq", "_restriction_upload_id")

    st.text_area("DNA Sequence:", height=140, key="restriction_seq")
    selected_enzymes = st.multiselect("Enzymes to check:", COMMON_ENZYMES, default=["EcoRI", "BamHI", "HindIII", "PstI"])

    if st.button("Find Sites", key="run_restriction", type="primary"):
        seq = clean_sequence(st.session_state.restriction_seq)

        if not seq:
            st.warning("Please provide a DNA sequence.")
        elif not selected_enzymes:
            st.warning("Please select at least one enzyme.")
        elif not set(seq) <= set("ACGTUN"):
            st.warning("This doesn't look like a DNA sequence.")
        else:
            try:
                rb = RestrictionBatch(selected_enzymes)
                result = Analysis(rb, Seq(seq)).full()

                st.divider()
                rows = [{
                    "Enzyme": str(enzyme), "Recognition Site": enzyme.site,
                    "Cut Count": len(positions),
                    "Positions": ", ".join(str(p) for p in positions) if positions else "—",
                } for enzyme, positions in result.items()]
                df = pd.DataFrame(rows).sort_values("Cut Count", ascending=False)
                st.dataframe(df, use_container_width=True)

                chart_df = df[df["Cut Count"] > 0].set_index("Enzyme")[["Cut Count"]]
                if not chart_df.empty:
                    st.subheader("Cut Site Counts")
                    st.bar_chart(chart_df)

                st.download_button("Download Results (CSV)", df.to_csv(index=False), file_name="restriction_sites.csv", mime="text/csv")
            except Exception as e:
                st.error(f"Restriction analysis failed: {e}")
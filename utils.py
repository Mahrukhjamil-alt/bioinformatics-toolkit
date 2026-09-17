import re
import streamlit as st

def clean_sequence(seq: str) -> str:
    """Strip whitespace, headers, and non-letter characters from a pasted sequence."""
    if not seq:
        return ""
    lines = seq.strip().splitlines()
    lines = [line for line in lines if not line.startswith(">")]
    joined = "".join(lines)
    return re.sub(r"[^A-Za-z]", "", joined).upper()


def parse_fasta(text: str) -> dict:
    """Parse multi-sequence FASTA text into an ordered {name: sequence} dict."""
    seqs = {}
    name = None
    buf = []
    for raw_line in text.strip().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = clean_sequence("".join(buf))
            name = line[1:].strip() or f"seq{len(seqs) + 1}"
            buf = []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = clean_sequence("".join(buf))
    return seqs
import streamlit as st  # add this import at the very top of the file


def load_uploaded_once(uploaded_file, session_key: str, tracker_key: str, multi: bool = False):
    """
    Load an uploaded file's content into st.session_state[session_key], but only once
    per unique file — so it doesn't overwrite manual edits on every rerun.
    """
    if uploaded_file is None:
        return
    identifier = f"{uploaded_file.name}-{uploaded_file.size}"
    if st.session_state.get(tracker_key) == identifier:
        return
    content = uploaded_file.read().decode("utf-8", errors="ignore")
    if multi:
        st.session_state[session_key] = content
    else:
        if ">" in content:
            parsed = parse_fasta(content)
            st.session_state[session_key] = next(iter(parsed.values())) if parsed else ""
        else:
            st.session_state[session_key] = clean_sequence(content)
    st.session_state[tracker_key] = identifier
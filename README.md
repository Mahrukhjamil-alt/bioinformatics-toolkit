# Bioinformatics Toolkit

🔗 **Live App:** [bioinformatics-toolkit.streamlit.app](https://bioinformatics-toolkit.streamlit.app/)

A Streamlit web app for common bioinformatics tasks...
# Bioinformatics Toolkit

A Streamlit web app for common bioinformatics tasks: sequence analysis, pairwise/multiple sequence alignment, phylogenetic trees, ORF finding, restriction site search, and motif search.

## Features
- Single Sequence Analysis (GC content, molecular weight, transcription/translation, reverse complement)
- Pairwise Sequence Alignment (BLOSUM62/50, PAM250, NUC.4.4)
- Multiple Sequence Alignment
- Phylogenetic Tree (UPGMA / Neighbor-Joining)
- ORF Finder & Codon Usage
- Restriction Enzyme Site Finder
- Motif / Pattern Search

## How to Run

### 1. Clone this repository
```bash
git clone https://github.com/Mahrukhjamil-alt/bioinformatics-toolkit.git
cd bioinformatics-toolkit
```

### 2. Create a virtual environment (recommended)
```bash
python -m venv venv
```
Activate it:
- Windows: `venv\Scripts\activate`
- Mac/Linux: `source venv/bin/activate`

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

### 4. Run the app
```bash
streamlit run app.py
```

This will open the app in your browser at `http://localhost:8501`.
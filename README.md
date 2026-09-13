# CRISPR Guide RNA Designer

A Python tool that designs and ranks candidate CRISPR/Cas9 single guide RNAs (sgRNAs) for any human gene, given just a gene symbol.

## What it does

1. **Fetch** — retrieves the coding mRNA sequence for a given gene from NCBI (via Biopython/Entrez)
2. **Scan** — identifies all candidate guide sites by locating NGG PAM sequences on both DNA strands
3. **Filter** — removes guides with unsuitable GC content (outside 40–60%) or poly-T runs (which terminate Pol III transcription)
4. **Score** — ranks remaining guides using a heuristic inspired by published sequence-efficiency features (Doench et al. 2016, Rule Set 2) — not a reimplementation of the trained model, but based on the same known biological signals (GC content, seed region composition, position-specific nucleotide preferences)

## Usage

### Command line
```bash
pip install -r requirements.txt
python src/fetch_sequence.py   # fetches example genes
python src/find_guides.py      # scans for PAM sites
python src/filter_guides.py    # filters by GC/poly-T
python src/score_guides.py     # scores and ranks
```

### Interactive app
```bash
python -m streamlit run app/streamlit_app.py
```
Enter any human gene symbol (e.g. `TP53`, `BRCA1`) and get a ranked, downloadable list of candidate guides with an interactive position plot.

## Limitations

- On-target scoring is a simplified heuristic, not the full trained Doench Rule Set 2 / Azimuth model
- No off-target scoring in this version (planned as a follow-up phase using genome alignment)
- Currently supports human genes with a fetchable RefSeq (NM_) mRNA record

## Tech stack

Python, Biopython, pandas, Streamlit, matplotlib

## Author

Anjali Bhatt — B.Tech Bioinformatics, IILM University
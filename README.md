# CRISPR Guide RNA Designer

A Python tool that designs and ranks candidate CRISPR/Cas9 single guide RNAs (sgRNAs) for any human gene, given just a gene symbol.

**Live app:** https://a-crispr-guide-designer.streamlit.app/


## Roadmap

- [ ] Off-target scoring: align each candidate guide against a reference genome (planned: Bowtie2 against a single chromosome or small model genome) and score using mismatch-based CFD scoring
- [ ] Support for additional Cas variants (Cas12a) with different PAM requirements
- [ ] Whole-genome off-target search
- [ ] Genome-absolute coordinates for on-target site detection

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


## Off-target scoring (Phase 2)

Off-target risk is precomputed for demo genes (TP53, BRCA1) and loaded by the app when available.

**Pipeline:** guides → Bowtie2 alignment against the chr17 reference (run in Colab) → SAM parsing → simplified CFD-inspired scoring → CSVs in `data/` → merged into the Streamlit table.

**Output columns:** `offtarget_hits`, `max_offtarget_risk`, `offtarget_status` (`scored`, `no off-targets found`, `high repeat content`, `not scored`).

**Why precomputed:** Bowtie2 does not run on Windows or typical deployment environments.


## Limitations

- On-target scoring is a simplified heuristic, not the full trained Doench Rule Set 2 / Azimuth model
- No off-target scoring in this version (planned as a follow-up phase using genome alignment)
- Currently supports human genes with a fetchable RefSeq (NM_) mRNA record
- Off-target search covers **chromosome 17 only**, not the whole genome. "No off-targets found" means none on chr17.
- Only the **top 20** on-target guides per demo gene are scored. Other genes show on-target scores only.
- On-target site detection is a placeholder: a guide with exactly one 0-mismatch hit is treated as its on-target site. A full fix needs converting mRNA-relative positions to genome coordinates.
- Off-target scoring is a simplified CFD-inspired heuristic, not validated for lab use.

## Tech stack

Python, Biopython, pandas, Streamlit, matplotlib

## Author

Anjali Bhatt — B.Tech Bioinformatics, IILM University
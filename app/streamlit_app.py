import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

APP_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(APP_DIR, "..", "src"))
DATA_DIR = os.path.join(APP_DIR, "..", "data")

from fetch_sequence import fetch_gene_sequence
from find_guides import find_guides
from filter_guides import filter_guides
from score_guides import score_all_guides

st.set_page_config(page_title="CRISPR Guide Designer", layout="wide")

OFFTARGET_COLS = ["guide_seq", "num_offtarget_hits", "max_offtarget_risk", "total_offtarget_risk"]
MAX_HITS_DISPLAY = 100  # above this, show ">100" and flag as repeat-like


@st.cache_data
def load_offtarget_scores(gene_symbol):
    """
    Load precomputed off-target scores if available for this gene.
    Returns (dataframe or None, warning message or None).
    """
    path = os.path.join(DATA_DIR, f"{gene_symbol}_guides_with_offtarget.csv")
    if not os.path.exists(path):
        return None, None
    try:
        df = pd.read_csv(path)
    except Exception as e:
        return None, f"Could not read the off-target file ({e}). Showing on-target scores only."
    missing = [c for c in OFFTARGET_COLS if c not in df.columns]
    if missing:
        return None, f"Off-target file is missing columns {missing}. Showing on-target scores only."
    # One row per guide sequence so the merge cannot duplicate rows
    df = df[OFFTARGET_COLS].drop_duplicates(subset="guide_seq")
    return df, None


def offtarget_status(hits):
    if pd.isna(hits):
        return "not scored"
    if hits > MAX_HITS_DISPLAY:
        return "high repeat content"
    if hits == 0:
        return "no off-targets found"
    return "scored"


def hits_display(hits):
    if pd.isna(hits):
        return ""
    if hits > MAX_HITS_DISPLAY:
        return f">{MAX_HITS_DISPLAY}"
    return str(int(hits))


st.title("🧬 CRISPR Guide RNA Designer")
st.write(
    "Enter a human gene symbol to find and rank candidate CRISPR/Cas9 guide RNAs. "
    "Guides are scanned for NGG PAM sites, filtered by GC content and poly-T runs, "
    "and scored using a heuristic inspired by Doench et al. 2016 (Rule Set 2). "
    "This is an educational/portfolio tool, not validated for lab use."
)

gene_symbol = st.text_input("Gene symbol (e.g. TP53, BRCA1)", "")

if st.button("Design Guides"):
    if not gene_symbol.strip():
        st.warning("Please enter a gene symbol.")
        st.stop()

    clean_symbol = gene_symbol.strip().upper()

    with st.spinner(f"Fetching sequence for {clean_symbol}..."):
        try:
            seq_record = fetch_gene_sequence(clean_symbol)
        except ValueError as e:
            st.error(f"Gene not found: {e}. Check the spelling, or try the official HGNC symbol.")
            st.stop()
        except Exception as e:
            st.error(f"Unexpected error while fetching from NCBI: {e}. This may be a temporary network/API issue — try again in a moment.")
            st.stop()

    if len(seq_record.seq) < 25:
        st.error(f"Fetched sequence is too short ({len(seq_record.seq)} bp) to design guides. This gene record may be incomplete.")
        st.stop()

    st.success(f"Fetched {seq_record.id} — {len(seq_record.seq)} bp")

    with st.spinner("Scanning for PAM sites..."):
        guide_df = find_guides(seq_record.seq)

    if guide_df.empty:
        st.warning("No candidate guides found — no NGG PAM sites detected in this sequence. This is unusual for a gene of typical length; double check the gene symbol.")
        st.stop()

    st.write(f"Found **{len(guide_df)}** candidate guides before filtering.")

    with st.spinner("Filtering by GC content and poly-T..."):
        filtered_df = filter_guides(guide_df)

    if filtered_df.empty:
        st.warning("No guides passed GC content / poly-T filtering. This can happen with very short or unusual sequences.")
        st.stop()

    with st.spinner("Scoring guides..."):
        scored_df = score_all_guides(filtered_df)

    offtarget_df, offtarget_warning = load_offtarget_scores(clean_symbol)
    if offtarget_warning:
        st.warning(offtarget_warning)

    has_offtarget = offtarget_df is not None

    if has_offtarget:
        scored_df = scored_df.merge(offtarget_df, on="guide_seq", how="left")
        n_scored = int(scored_df["num_offtarget_hits"].notna().sum())
        st.info(
            f"Off-target scores available for {clean_symbol} "
            f"(precomputed against chromosome 17, top 20 guides by on-target score). "
            f"{n_scored} of {len(scored_df)} listed guides have off-target results."
        )
        scored_df["offtarget_hits"] = scored_df["num_offtarget_hits"].apply(hits_display)
        scored_df["offtarget_status"] = scored_df["num_offtarget_hits"].apply(offtarget_status)
    else:
        st.caption(f"Off-target scoring not precomputed for {clean_symbol} — showing on-target scores only.")

    st.write(f"**{len(scored_df)}** guides passed filtering, ranked by predicted efficiency below:")

    base_cols = ["guide_seq", "strand", "position", "pam_seq", "gc_content", "efficiency_score"]
    if has_offtarget:
        display_df = scored_df[base_cols + ["offtarget_hits", "max_offtarget_risk", "offtarget_status"]]
    else:
        display_df = scored_df[base_cols]

    st.dataframe(display_df, width="stretch")

    csv = display_df.to_csv(index=False)
    st.download_button(
        "Download results as CSV",
        data=csv,
        file_name=f"{clean_symbol}_guides.csv",
        mime="text/csv"
    )

    st.subheader("Guide positions along the gene")
    fig, ax = plt.subplots(figsize=(10, 2))
    plus_guides = scored_df[scored_df["strand"] == "+"]
    minus_guides = scored_df[scored_df["strand"] == "-"]
    ax.scatter(plus_guides["position"], [1] * len(plus_guides), marker="|", color="blue", label="+ strand")
    ax.scatter(minus_guides["position"], [0] * len(minus_guides), marker="|", color="red", label="- strand")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["- strand", "+ strand"])
    ax.set_xlabel("Position (bp)")
    ax.legend(loc="upper right")
    st.pyplot(fig)
    plt.close(fig)
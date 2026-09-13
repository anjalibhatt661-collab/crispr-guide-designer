import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))

from fetch_sequence import fetch_gene_sequence
from find_guides import find_guides
from filter_guides import filter_guides
from score_guides import score_all_guides

st.set_page_config(page_title="CRISPR Guide Designer", layout="wide")

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

    st.write(f"**{len(filtered_df)}** guides passed filtering, ranked by predicted efficiency below:")

    display_df = scored_df[["guide_seq", "strand", "position", "pam_seq", "gc_content", "efficiency_score"]]
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
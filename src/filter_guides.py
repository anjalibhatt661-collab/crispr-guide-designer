import pandas as pd

def calculate_gc_content(seq):
    """Return GC percentage of a sequence."""
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return (gc_count / len(seq)) * 100


def has_poly_t(seq, max_t_run=4):
    """Check if sequence contains a poly-T stretch (terminates Pol III transcription)."""
    return "T" * max_t_run in seq.upper()


def filter_guides(guide_df, gc_min=40, gc_max=60, max_t_run=4):
    """
    Filter candidate guides by GC content range and poly-T presence.
    Adds gc_content column, then filters.
    """
    df = guide_df.copy()
    df["gc_content"] = df["guide_seq"].apply(calculate_gc_content)
    df["has_polyt"] = df["guide_seq"].apply(lambda s: has_poly_t(s, max_t_run))

    filtered = df[
        (df["gc_content"] >= gc_min) &
        (df["gc_content"] <= gc_max) &
        (~df["has_polyt"])
    ].reset_index(drop=True)

    return filtered


if __name__ == "__main__":
    for gene in ["TP53", "BRCA1"]:
        print(f"Filtering guides for {gene}...")
        guide_df = pd.read_csv(f"data/{gene}_guides.csv")
        print(f"  Starting guides: {len(guide_df)}")

        filtered_df = filter_guides(guide_df)
        print(f"  After GC + poly-T filtering: {len(filtered_df)}")
        print(filtered_df.head())

        filtered_df.to_csv(f"data/{gene}_guides_filtered.csv", index=False)
        print()
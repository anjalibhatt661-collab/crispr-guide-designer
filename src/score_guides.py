import pandas as pd

def score_guide(guide_seq, pam_seq):
    """
    Simplified on-target efficiency score, inspired by published
    sequence features from Doench et al. 2016 (Rule Set 2) — NOT
    a reimplementation of the trained model, but a heuristic using
    the same general biological signals, with finer granularity
    across the full 20bp guide rather than just a few checkpoints.
    """
    guide_seq = guide_seq.upper()
    score = 50.0

    # 1. GC content: peak around 50%, penalty scaled by distance from it
    gc = (guide_seq.count("G") + guide_seq.count("C")) / len(guide_seq) * 100
    score -= abs(gc - 50) * 0.4

    # 2. Position-specific nucleotide contribution across ALL 20 positions.
    #    Bases closer to the PAM (3' end, the "seed region") matter more
    #    for Cas9 binding, so weight increases toward position 20.
    for i, base in enumerate(guide_seq):
        position_weight = (i + 1) / 20  # 0.05 at position 1, up to 1.0 at position 20
        if base == "G":
            score += 0.6 * position_weight
        elif base == "C":
            score += 0.3 * position_weight
        elif base == "T":
            score -= 0.5 * position_weight
        elif base == "A":
            score += 0.1 * position_weight

    # 3. Specific known signals from literature
    if guide_seq[19] == "G":       # position 20, right before PAM
        score -= 4
    if guide_seq[0] == "T":        # position 1, 5' end
        score -= 3

    # 4. Penalize any run of 3+ identical bases (secondary structure risk)
    for base in "ACGT":
        if base * 3 in guide_seq:
            score -= 2

    score = max(0, min(100, score))
    return round(score, 2)


def score_all_guides(guide_df):
    df = guide_df.copy()
    df["efficiency_score"] = df.apply(
        lambda row: score_guide(row["guide_seq"], row["pam_seq"]), axis=1
    )
    return df.sort_values("efficiency_score", ascending=False).reset_index(drop=True)


if __name__ == "__main__":
    for gene in ["TP53", "BRCA1"]:
        print(f"Scoring guides for {gene}...")
        filtered_df = pd.read_csv(f"data/{gene}_guides_filtered.csv")

        scored_df = score_all_guides(filtered_df)
        print(f"  Scored {len(scored_df)} guides")
        print(f"  Unique scores: {scored_df['efficiency_score'].nunique()} out of {len(scored_df)}")
        print(scored_df[["guide_seq", "pam_seq", "gc_content", "efficiency_score"]].head(10))

        scored_df.to_csv(f"data/{gene}_guides_scored.csv", index=False)
        print()
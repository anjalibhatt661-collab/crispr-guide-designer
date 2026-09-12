from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re

def find_guides(sequence, pam="NGG", guide_length=20):
    """
    Scan a DNA sequence (both strands) for PAM sites and extract
    candidate guide RNA sequences.

    Returns a DataFrame with columns:
        guide_seq, strand, position, pam_seq
    """
    sequence = str(sequence).upper()
    guides = []

    # --- Sense strand: look for NGG ---
    pam_pattern = pam.replace("N", ".")
    for match in re.finditer(f"(?=({pam_pattern}))", sequence):
        pam_start = match.start()
        guide_start = pam_start - guide_length
        if guide_start >= 0:
            guide_seq = sequence[guide_start:pam_start]
            pam_seq = sequence[pam_start:pam_start + len(pam)]
            guides.append({
                "guide_seq": guide_seq,
                "strand": "+",
                "position": guide_start,
                "pam_seq": pam_seq
            })

    # --- Antisense strand: look for CCN, which is NGG on the reverse complement ---
    rev_comp = str(Seq(sequence).reverse_complement())
    for match in re.finditer(f"(?=({pam_pattern}))", rev_comp):
        pam_start = match.start()
        guide_start = pam_start - guide_length
        if guide_start >= 0:
            guide_seq = rev_comp[guide_start:pam_start]
            pam_seq = rev_comp[pam_start:pam_start + len(pam)]
            original_pos = len(sequence) - pam_start
            guides.append({
                "guide_seq": guide_seq,
                "strand": "-",
                "position": original_pos,
                "pam_seq": pam_seq
            })

    return pd.DataFrame(guides)


if __name__ == "__main__":
    for gene in ["TP53", "BRCA1"]:
        print(f"Scanning {gene} for guides...")
        record = SeqIO.read(f"data/{gene}.fasta", "fasta")
        guide_df = find_guides(record.seq)
        print(f"  Found {len(guide_df)} candidate guides")
        print(guide_df.head())
        guide_df.to_csv(f"data/{gene}_guides.csv", index=False)
        print()
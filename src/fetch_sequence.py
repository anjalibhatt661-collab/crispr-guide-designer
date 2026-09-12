from Bio import Entrez, SeqIO

Entrez.email = "anjalibhatt661@gmail.com"

def fetch_gene_sequence(gene_symbol, organism="Homo sapiens"):
    """
    Fetch the coding (NM_) mRNA sequence for a gene by symbol.
    Returns a Bio.SeqRecord object.
    """
    search_term = f"{gene_symbol}[Gene Name] AND {organism}[Organism]"
    handle = Entrez.esearch(db="gene", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError(f"No gene found for '{gene_symbol}' in {organism}")

    gene_id = record["IdList"][0]

    handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id,
                           linkname="gene_nuccore_refseqrna")
    link_record = Entrez.read(handle)
    handle.close()

    if not link_record[0]["LinkSetDb"]:
        raise ValueError(f"No RefSeq mRNA found for gene ID {gene_id}")

    nuccore_ids = [link["Id"] for link in link_record[0]["LinkSetDb"][0]["Link"]]

    # Fetch summaries first to find which ID is an NM_ (coding) accession
    handle = Entrez.esummary(db="nuccore", id=",".join(nuccore_ids))
    summaries = Entrez.read(handle)
    handle.close()

    nm_id = None
    for s in summaries:
        if s["Caption"].startswith("NM_"):
            nm_id = s["Id"]
            break

    # Fall back to the first result if no NM_ found
    chosen_id = nm_id if nm_id else nuccore_ids[0]

    handle = Entrez.efetch(db="nuccore", id=chosen_id, rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()

    return seq_record


if __name__ == "__main__":
    for gene in ["TP53", "BRCA1"]:
        print(f"Fetching {gene}...")
        try:
            seq = fetch_gene_sequence(gene)
            print(f"  ID: {seq.id}")
            print(f"  Length: {len(seq.seq)} bp")
            print(f"  First 60 bp: {seq.seq[:60]}")
            SeqIO.write(seq, f"data/{gene}.fasta", "fasta")
        except Exception as e:
            print(f"  Failed: {e}")
        print()
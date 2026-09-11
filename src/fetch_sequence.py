from Bio import Entrez, SeqIO

# NCBI requires you to identify yourself
Entrez.email = "anjalibhatt661@gmail.com"

def fetch_gene_sequence(gene_symbol, organism="Homo sapiens"):
    """
    Fetch the coding sequence for a gene by symbol.
    Returns a Bio.SeqRecord object.
    """
    # Step 1: search Gene database to get the Gene ID
    search_term = f"{gene_symbol}[Gene Name] AND {organism}[Organism]"
    handle = Entrez.esearch(db="gene", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError(f"No gene found for '{gene_symbol}' in {organism}")

    gene_id = record["IdList"][0]

    # Step 2: link Gene ID to nucleotide records, get the RefSeq mRNA
    handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id,
                           linkname="gene_nuccore_refseqrna")
    link_record = Entrez.read(handle)
    handle.close()

    if not link_record[0]["LinkSetDb"]:
        raise ValueError(f"No RefSeq mRNA found for gene ID {gene_id}")

    nuccore_id = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]

    # Step 3: fetch the actual sequence
    handle = Entrez.efetch(db="nuccore", id=nuccore_id, rettype="fasta", retmode="text")
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
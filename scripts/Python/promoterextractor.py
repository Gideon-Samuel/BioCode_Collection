from Bio import Entrez, SeqIO

#Identify ourselves to NCBI
Entrez.email = "gideonsamuelj23@gmail.com"

#Fetch a gene sequence (example: TP53, genomic region)
accession = "NG_017013.2"  # TP53 genomic sequence
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq).upper()
gene_start = 1000  # Example start position of the gene in the sequence
promoter_length = 500  # Number of bp upstream to extract

#Extract promoter region (upstream of gene start)
promoter_seq = sequence[max(0, gene_start-promoter_length):gene_start]

#Print results
print(f"Gene: TP53 ({accession})")
print(f"Promoter region (last {promoter_length} bp upstream):")
print(promoter_seq)
print(f"Length: {len(promoter_seq)} bp")

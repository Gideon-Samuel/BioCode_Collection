from Bio import Entrez, SeqIO

# Tell NCBI who you are, by your email
Entrez.email = "gideonsamuelj23@gmail.com"

# Sequence to fetch
accession = "NC_000001.11"

# Fetch and parse FASTA
record = SeqIO.read(Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text"), "fasta")

# Save to file
SeqIO.write(record, accession + ".fa", "fasta")

# Print summary
print("ID:", record.id, "Length:", len(record.seq), "bp")
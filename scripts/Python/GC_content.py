from Bio import Entrez, SeqIO

# Tell NCBI who we are via email!
Entrez.email = "gideonsamuelj23@gmail.com"

# Human mitochondrial DNA
accession = "NC_012920.1"

# Fetch and parse FASTA
record = SeqIO.read(
    Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text"),
    "fasta"
)

# Save to file
SeqIO.write(record, accession + ".fa", "fasta")

# Calculate GC content
gc_count = record.seq.count("G") + record.seq.count("C")
gc_content = (gc_count / len(record.seq)) * 100

# Print result
print("ID:", record.id)
print("Length:", len(record.seq), "bp")
print("GC content:", round(gc_content, 2), "%")

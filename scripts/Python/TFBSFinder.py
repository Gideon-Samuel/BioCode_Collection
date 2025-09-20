from Bio import Entrez, SeqIO

# Step 1. Identify ourselves to NCBI
Entrez.email = "gideonsamuelj23@gmail.com"

# Step 2. Fetch a real DNA sequence (TP53 promoter region, example accession)
accession = "NG_017013.2"  
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq)

# Step 3. Define transcription factor motifs
motifs = {
    "TATA_box": "TATAAA",
    "SP1": "GGGCGG",
    "CAAT_box": "GGCCAATCT"
}

# Step 4. Search motifs in the real sequence
for name, motif in motifs.items():
    if motif in sequence:
        print(f"{name} found at position {sequence.find(motif)} in {accession}")
    else:
        print(f"{name} not found in {accession}")
  
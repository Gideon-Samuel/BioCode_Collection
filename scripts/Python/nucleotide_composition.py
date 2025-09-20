from Bio import Entrez, SeqIO

Entrez.email = "gideonsamuelj23@gmail.com"
accession = "NM_000518.5"  # Human HBB gene
record = SeqIO.read(
    Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text"),
    "fasta"
)
SeqIO.write(record, accession + ".fa", "fasta")
A = record.seq.count("A")
T = record.seq.count("T")
G = record.seq.count("G")
C = record.seq.count("C")
total = len(record.seq)
print("ID:", record.id)
print("Length:", total, "bp")
print("A:", round(A/total*100, 2), "%")
print("T:", round(T/total*100, 2), "%")
print("G:", round(G/total*100, 2), "%")
print("C:", round(C/total*100, 2), "%")

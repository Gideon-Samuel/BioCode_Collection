from Bio import Entrez, SeqIO

#Identify ourselves to NCBI
Entrez.email = "gideonsamuelj23@gmail.com"

#Fetch real sequence (TP53 promoter)
accession = "NG_017013.2"  
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq).upper()  # Convert to uppercase

#Compute basic statistics
length = len(sequence)
gc_count = sequence.count("G") + sequence.count("C")
gc_content = (gc_count / length) * 100

#Nucleotide counts
a_count = sequence.count("A")
t_count = sequence.count("T")
g_count = sequence.count("G")
c_count = sequence.count("C")

#Print results
print(f"Sequence ID: {record.id}")
print(f"Description: {record.description}")
print(f"Length: {length} bp")
print(f"GC Content: {gc_content:.2f}%")
print(f"Nucleotide Counts: A={a_count}, T={t_count}, G={g_count}, C={c_count}")

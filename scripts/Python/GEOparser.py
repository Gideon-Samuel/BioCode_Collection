import GEOparse

# Download a small GEO dataset
gse = GEOparse.get_GEO(geo="GSE70970", destdir="./data", silent=True)

# List all sample IDs
samples = list(gse.gsms.keys())
print("Samples in this GEO dataset:")
for s in samples[:20]:
    print("-", s)

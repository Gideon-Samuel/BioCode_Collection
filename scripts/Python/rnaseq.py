import GEOparse
import pandas as pd
import matplotlib.pyplot as plt

#Download GEO dataset
gse = GEOparse.get_GEO("GSE70970", destdir="./data", silent=True)

#Pick first 4 samples (2 Control, 2 Treated)
samples = list(gse.gsms.keys())[:4]

#Build DataFrame with total counts per sample
counts = {}
for s in samples:
    counts[s] = gse.gsms[s].table["VALUE"].sum()  # sum of all gene counts

counts_df = pd.DataFrame.from_dict(counts, orient="index", columns=["Total_Reads"])

#Plot bar chart
plt.figure(figsize=(6,4))
plt.bar(counts_df.index, counts_df["Total_Reads"], color=["skyblue", "skyblue", "orange", "orange"])
plt.ylabel("Total Read Counts")
plt.title("RNA-seq QC: Total Reads per Sample")
plt.xticks(rotation=45)
plt.show()

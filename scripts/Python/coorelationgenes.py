import GEOparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Download GEO dataset
gse = GEOparse.get_GEO("GSE70970", destdir="./data", silent=True)

#Pick first 4 samples (2 Control, 2 Treated)
samples = list(gse.gsms.keys())[:4]

#Pick first 5 probes (genes) for demonstration
probes = list(gse.gsms[samples[0]].table.index)[:5]

#Build DataFrame: rows = genes, columns = samples
data = pd.DataFrame(index=probes)
for s in samples:
    data[s] = gse.gsms[s].table.loc[probes, "VALUE"]

#Compute correlation between genes
corr = data.corr()  # correlation across samples

#Plot heatmap
sns.heatmap(corr, annot=True, cmap="coolwarm")
plt.title("Gene-Gene Correlation Heatmap (Real GEO Data)")
plt.show()

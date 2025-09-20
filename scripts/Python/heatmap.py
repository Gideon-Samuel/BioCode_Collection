import GEOparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#Downloading GEO Dataset
gse = GEOparse.get_GEO("GSE70970", destdir="./data", silent=True)

# Picking First 4 samples
samples = list(gse.gsms.keys())[:4]
groups = ["Control", "Control", "Treated", "Treated"]

#Picking 5 probes
probes = list(gse.gsms[samples[0]].table.index)[:5]

#Building dataframe
data = pd.DataFrame(index=probes)
for s in samples:
    data[s] = gse.gsms[s].table.loc[probes,"VALUE"]

col_colors = ["skyblue" if g=="Control" else "salmon" for g in groups]

#Plotting Heatmap
sns.clustermap(data, cmap="coolwarm", col_colors=col_colors, figsize=(6,5))
plt.title("Gene Expression Heatmap with Groups", pad=40)
plt.show()
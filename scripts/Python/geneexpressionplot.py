#Importing Modules for expression analysis and Visualization
import GEOparse
import matplotlib.pyplot as plt

# Downloading GEO dataset
gse = GEOparse.get_GEO("GSE70970", destdir="./data", silent=True)

# Picking first 20 samples
samples = list(gse.gsms.keys())[:20]

# Automatically choosing a valid probe
probe = list(gse.gsms[samples[0]].table.index)[0]

# Get expression values
values = [float(gse.gsms[s].table.loc[probe, "VALUE"]) for s in samples]

# Simple bar plot
plt.bar(samples, values)
plt.ylabel("Expression")
plt.title(f"{probe}th Expression")
plt.xticks(rotation=90)
plt.show()

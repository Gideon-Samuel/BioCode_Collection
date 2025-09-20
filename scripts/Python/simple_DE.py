import GEOparse
import pandas as pd
import matplotlib.pyplot as plt

# Downloading GEO dataset
gse = GEOparse.get_GEO("GSE70970", destdir="./data", silent=True)

# Pick first 4 samples: first 2 = control, next 2 = treated
samples = list(gse.gsms.keys())[:4]
control_samples = samples[:2]
treated_samples = samples[2:]

# Pick a probe automatically
probe = list(gse.gsms[samples[0]].table.index)[0]

# Build a simple DataFrame
data = {
    "Sample": samples,
    "Group": ["Control", "Control", "Treated", "Treated"],
    "Expression": [float(gse.gsms[s].table.loc[probe, "VALUE"]) for s in samples]
}
df = pd.DataFrame(data)

# Calculate average expression per group
avg_expr = df.groupby("Group")["Expression"].mean()
fold_change = avg_expr["Treated"] / avg_expr["Control"]

# Print results
print(f"Probe: {probe}")
print(df)
print("\nAverage Expression per Group:")
print(avg_expr)
print(f"\nFold Change (Treated/Control): {round(fold_change,2)}")

#Plotting Bar
colors = ["skyblue" if g=="Control" else "salmon" for g in df["Group"]]

plt.bar(df["Sample"], df["Expression"], color=colors)
plt.ylabel("Expression Value")
plt.title(f"Expression of Probe {probe}")
plt.xticks(rotation=45)
plt.show()
import GEOparse

# Download GEO dataset
gse = GEOparse.get_GEO("GSE70970", destdir="./data", silent=True)

# Pick samples 3–6 for this script
samples = list(gse.gsms.keys())[2:6]

# Split into control and treated (first 2 = control, next 2 = treated)
control_samples = samples[:2]
treated_samples = samples[2:]

# Pick a different probe: second probe in first sample
probe = list(gse.gsms[samples[0]].table.index)[1]

# Get expression values
control_values = [float(gse.gsms[s].table.loc[probe, "VALUE"]) for s in control_samples]
treated_values = [float(gse.gsms[s].table.loc[probe, "VALUE"]) for s in treated_samples]

# Calculate averages
avg_control = sum(control_values) / len(control_values)
avg_treated = sum(treated_values) / len(treated_values)

# Fold change (Treated / Control)
fold_change = avg_treated / avg_control

# Print results
print(f"Probe: {probe}")
print(f"Average Control: {avg_control}")
print(f"Average Treated: {avg_treated}")
print(f"Fold Change (Treated / Control): {round(fold_change, 2)}")

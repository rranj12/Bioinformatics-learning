

import pandas as pd
import matplotlib.pyplot as plt

# Load gene counts (2 columns: gene_id, count)
df = pd.read_csv("results/gene_counts.tsv", sep="\t", header=None, names=["gene_id", "count"])

# Basic sanity
df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)
df = df.sort_values("count", ascending=False)

print("Num genes:", len(df))
print("Total assigned reads (sum counts):", int(df["count"].sum()))
print("\nTop 20 genes by count:")
print(df.head(20).to_string(index=False))

# Plot top 50
top = df.head(50).copy()
plt.figure(figsize=(10, 6))
plt.bar(range(len(top)), top["count"])
plt.xticks(range(len(top)), top["gene_id"], rotation=90, fontsize=6)
plt.ylabel("Read count")
plt.title("Top 50 genes by counts (featureCounts)")
plt.tight_layout()
plt.savefig("results/top50_genes_bar.png", dpi=200)
print("\nWrote: results/top50_genes_bar.png")


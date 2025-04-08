import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Original data
unique_clusters = {
    "Zinc Finger": 161,
    "Homeodomain": 151,
    "Basic Domain": 112,
    "Helix-turn-helix": 161,
    "Unidentified DNA Binding Domain": 95,
    "High Motility Group": 55,
}

# Total value for percentage calculation
total_clusters = 161

# Calculate percent_unique
percent_unique = {
    key: (value / total_clusters) * 100 for key, value in unique_clusters.items()
}

print(percent_unique)
quit()

# Convert to a DataFrame for Seaborn
df = pd.DataFrame(
    list(percent_unique.items()),
    columns=["Transcription Factor Class", "Percent Unique"],
)
df = df.sort_values(by="Percent Unique", ascending=False)

my_colors = ["#0bab2e", "#be27ae", "#cf201f", "#0162ce", "#dbb8fe", "#af66f0"]

plt.rcParams["font.family"] = "Helvetica"
sns.barplot(
    x="Percent Unique",
    y="Transcription Factor Class",
    data=df,
    palette=my_colors,
    edgecolor="black",
    linewidth=0.5,
)

"Percent Unique"
# Add labels and a title
plt.xlabel("Percent Unique")
plt.ylabel("Transcription Factor Class")
plt.title("Unique Clusters")

ax = plt.gca()  # Get the current axis
ax.tick_params(axis="x", width=0.5)  # Set width of x-axis ticks
ax.tick_params(axis="y", width=0.5)
ax.spines["top"].set_visible(False)  # Remove the top spine
ax.spines["right"].set_visible(False)  # Remove the right spine
ax.spines["bottom"].set_linewidth(0.5)  # Set bottom axis line width
ax.spines["left"].set_linewidth(0.5)


plt.tight_layout()
plt.savefig("/Users/derek/Desktop/RNAseq/unique_tf_clusters_plot.pdf")
plt.show()

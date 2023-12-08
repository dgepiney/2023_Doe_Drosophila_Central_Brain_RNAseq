import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json
import itertools


data_path = "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/"
df = pd.read_csv(
  data_path + "zftfs_binary_gene_avg.csv"
)

clusters = list(range(175))
glia = (
  7, 13, 14, 15, 17, 18, 20, 23, 31, 34, 38, 46, 70, 80, 87, 92, 142, 146, 147, 
  174
)
neuron_clusters = [num for num in clusters if num not in glia]


with open(data_path + "hdtf_combos.json") as f:
  cluster_combos_list = json.load(f)

for i, cluster_combos in enumerate(cluster_combos_list): 
  cluster_combos_list[i]["Clusters"]=[
    int(s) 
    for s 
    in cluster_combos["Clusters"]
  ]

clusters_list_of_lists = [
  clusters["Clusters"] 
  for clusters 
  in cluster_combos_list
]


cluster_combos = list(itertools.chain(*clusters_list_of_lists))


for i in neuron_clusters:
  if i not in cluster_combos:
    cluster_combos.append(i)
print(cluster_combos)


custom_column_order = cluster_combos
custom_column_order_str = [str(e) for e in custom_column_order] 
custom_column_order_str.append("Gene")
custom_column_order_str = [f"{t:0>3}" for t in custom_column_order_str]
print(custom_column_order)


def indexer(e): 
  print(e)
  print(custom_column_order_str.index(e))
  return custom_column_order_str.index(e)


sums = np.array(df.iloc[:,1:].sum(axis=1))


df.insert(
  loc=0,
  column="temp_sum",
  value=sums
)

df = df.sort_values(by="temp_sum")

print(df["temp_sum"])


df = df.drop(
  labels="temp_sum",
  axis=1
)

df = df.sort_index(
  axis=1, 
  key=np.vectorize(indexer)
)

# Set "Gene" column as the index of the DataFrame
df.set_index("Gene", inplace=True)

# Replace zeros with NaN to handle missing values
df = df.replace(0, np.nan)

# Create the heatmap in black and white
plt.figure(figsize=(40, 20))  # Set the figure size
plt.xticks(fontsize=10)
plt.yticks(fontsize=10) 
ax = sns.heatmap(
  df, cmap="gray",
  cbar=False,
  linewidths = 0.5,
  linecolor="gray",
  square=True,
  xticklabels=True,
  yticklabels=True
)
ax.tick_params(
  axis="both",
  length=3,
  width=0.5,
  pad=1
)
plt.savefig(
  "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq"
  "/T2_zftfs_heatmap.pdf"
)
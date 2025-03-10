import glob
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# seurat clusters
clusters = (
    0,
    2,
    3,
    4,
    5,
    6,
    7,
    12,
    13,
    14,
    15,
    18,
    19,
    20,
    21,
    24,
    25,
    27,
    29,
    30,
    31,
    32,
    33,
    34,
    36,
    38,
    39,
    40,
    41,
    42,
    43,
    45,
    46,
    47,
    48,
    49,
    50,
    51,
    52,
    53,
    54,
    55,
    56,
    58,
    59,
    60,
    62,
    64,
    65,
    66,
    67,
    68,
    69,
    70,
    71,
    72,
    73,
    74,
    75,
    76,
    77,
    79,
    80,
    81,
    82,
    83,
    84,
    85,
    86,
    87,
    88,
    89,
    90,
    91,
    92,
    93,
    95,
    96,
    97,
    98,
    100,
    101,
    102,
    103,
    105,
    106,
    107,
    108,
    109,
    111,
    112,
    113,
    114,
    115,
    117,
    118,
    119,
    120,
    122,
    123,
    124,
    126,
    127,
    128,
    129,
    130,
    131,
    132,
    133,
    135,
    136,
    138,
    139,
    140,
    141,
    142,
    143,
    146,
    147,
    148,
    149,
    150,
    151,
    152,
    153,
    154,
    155,
    156,
    157,
    158,
    160,
    161,
    162,
    163,
    164,
    165,
    166,
    167,
    168,
    169,
    170,
    172,
    173,
    174,
    175,
    177,
    180,
    181,
    182,
    183,
    185,
    186,
    187,
    188,
    189,
    190,
    191,
    192,
    193,
    195,
    197,
)

# data path
processed_file_path = "/Users/derek/Desktop/RNAseq"

# ensure all csvs are processed
processed_file_pattern = f"{processed_file_path}/*.csv"
processed_files = glob.glob(processed_file_pattern)

# set output folder
output_folder = "/Users/derek/Desktop/RNAseq/heatmaps"


for file in processed_files:
    base_name = os.path.basename(file)
    split_name = base_name.split("_")
    print(f"Processing {file}")
    df = pd.read_csv(file)

    # split and explode the 'clusters' column
    df["clusters"] = df["clusters"].str.split(",")
    expanded_df = df.explode("clusters", ignore_index=True)

    # convert clusters to int
    expanded_df["clusters"] = expanded_df["clusters"].astype(int)

    # create a binary matrix
    binary_matrix = pd.DataFrame(
        0, index=expanded_df["gene"].unique(), columns=clusters
    )

    # fill the binary matrix with 1s for significant genes in each cluster
    for _, row in expanded_df.iterrows():
        if row["clusters"] in clusters:
            binary_matrix.at[row["gene"], row["clusters"]] = 1

    # sort rows by row sums
    row_sums = binary_matrix.sum(axis=1)
    sorted_rows_by_gene_frequency = row_sums.sort_values(ascending=True).index
    binary_matrix_sorted_rows = binary_matrix.loc[sorted_rows_by_gene_frequency]

    # compute the Jaccard similarity between each cluster
    def jaccard_similarity(set1, set2):
        union_length = len(set1.union(set2))
        if union_length == 0:
            return 0  # Skip similarity calculation and return 0 instead of raising an error
        return len(set1.intersection(set2)) / union_length

    # get gene combination for each cluster
    cluster_combinations = {}
    for cluster in clusters:
        genes_in_cluster = binary_matrix.index[binary_matrix[cluster] == 1].tolist()
        cluster_combinations[cluster] = frozenset(genes_in_cluster)

    # create a similarity matrix
    similarity_matrix = pd.DataFrame(index=clusters, columns=clusters, dtype=float)

    # ensure each cluster is counted only once
    for cluster1 in clusters:
        for cluster2 in clusters:
            if cluster1 < cluster2:
                similarity_matrix.at[cluster1, cluster2] = jaccard_similarity(
                    cluster_combinations[cluster1], cluster_combinations[cluster2]
                )
            else:
                similarity_matrix.at[cluster1, cluster2] = np.NaN

    # print clusters with Jaccard index of 1
    for cluster1 in clusters:
        for cluster2 in clusters:
            if cluster1 < cluster2:
                if similarity_matrix.at[cluster1, cluster2] == 1:
                    print(
                        f"Cluster {cluster1} and Cluster {cluster2} have a Jaccard index of 1."
                    )

    # create sorted cluster list
    sorted_clusters_by_similarity = []
    visited_clusters = set()

    # add the pairs with Jaccard index of 1 to the sorted list
    for cluster1 in clusters:
        for cluster2 in clusters:
            if cluster1 < cluster2:
                if similarity_matrix.at[cluster1, cluster2] == 1:
                    if cluster1 not in visited_clusters:
                        sorted_clusters_by_similarity.append(cluster1)
                        visited_clusters.add(cluster1)
                    if cluster2 not in visited_clusters:
                        sorted_clusters_by_similarity.append(cluster2)
                        visited_clusters.add(cluster2)

    # find the most similar cluster
    iteration_counter = 0
    while len(visited_clusters) < len(clusters):
        iteration_counter += 1

        # check if the number of iterations exceeds a threshold
        if iteration_counter > len(clusters) * len(clusters):
            print(
                f"Warning: Exceeded maximum iterations while sorting clusters in {file}"
            )
            break
        max_similarity = -1
        best_pair = None
        for cluster1 in clusters:
            if cluster1 in visited_clusters:
                continue
            for cluster2 in clusters:
                if cluster2 in visited_clusters or cluster1 == cluster2:
                    continue
                similarity = similarity_matrix.at[cluster1, cluster2]
                if similarity > max_similarity:
                    max_similarity = similarity
                    best_pair = (cluster1, cluster2)

        if best_pair:
            cluster1, cluster2 = best_pair
            sorted_clusters_by_similarity.append(cluster1)
            sorted_clusters_by_similarity.append(cluster2)
            visited_clusters.add(cluster1)
            visited_clusters.add(cluster2)

    # sort columns based on similarity
    binary_matrix_sorted_columns = binary_matrix_sorted_rows[
        sorted_clusters_by_similarity
    ]

    # create heatmap
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.family"] = "Helvetica"
    plt.figure(figsize=(40, 20))
    plt.xticks(fontsize=4)
    plt.yticks(fontsize=4)
    ax = sns.heatmap(
        binary_matrix_sorted_columns,
        cmap="gray_r",
        cbar=False,
        linewidths=0.25,
        linecolor="gray",
        square=True,
        xticklabels=True,
        yticklabels=True,
    )
    ax.tick_params(axis="both", length=3, width=0.5, pad=1)

    # save heatmap
    plt.savefig(f"{output_folder}/{split_name[0]}_heatmap.pdf")

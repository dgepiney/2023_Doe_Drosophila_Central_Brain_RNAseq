import subprocess
import glob
import os
import shutil
import pandas as pd
import numpy as np

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

# file directory
processed_file_directory = "/Users/derek/Desktop/RNAseq/processed_files"

# file patterns
file_pattern = "*.csv"
proccessed_files = glob.glob(os.path.join(processed_file_directory, file_pattern))

# write number of unique clusters to text file
with open("unique_clusters_count.txt", "w") as f:
    for file in proccessed_files:
        print(f"Processing file: {file}")

        df = pd.read_csv(file)

        # process clusters column
        df["clusters"] = df["clusters"].str.split(",")
        expanded_df = df.explode("clusters", ignore_index=True)
        expanded_df["clusters"] = expanded_df["clusters"].astype(int)

        # create a dictionary to store gene combinations for each cluster
        cluster_genes = {cluster: set() for cluster in clusters}

        # add genes in each cluster to dictionary
        for _, row in expanded_df.iterrows():
            cluster_genes[row["clusters"]].add(row["gene"])

        # create a dictionary to keep track of gene combinations
        combination_counts = {}

        # loop through clusters to create gene combinations
        for cluster, genes in cluster_genes.items():
            gene_combination = tuple(sorted(genes))
            if gene_combination:
                if gene_combination in combination_counts:
                    combination_counts[gene_combination].append(cluster)
                else:
                    combination_counts[gene_combination] = [cluster]

        # create sets for unique and non-unique clusters
        unique_clusters = set()
        non_unique_clusters = set()

        # classify clusters as unique or non-unique
        for combination, clusters_with_combination in combination_counts.items():
            if len(clusters_with_combination) == 1:
                unique_clusters.update(clusters_with_combination)
            else:
                non_unique_clusters.update(clusters_with_combination)

        # print number of unique and non-unique clusters
        unique_cluster_count = len(unique_clusters)
        print(f"Number of unique clusters: {unique_cluster_count}")

        # write number unique clusters to file
        f.write(f"{os.path.basename(file)}: {unique_cluster_count} unique clusters\n")

import pandas as pd
import itertools
import numpy as np
import math

# takes output of binarization_test.r file 
# removes decimals and x's r added to column names,
# adds 'cluster_' prefix to all column names, 
# averages binary score of each gene per cluster
# and writes the output to a csv file.


#read r binarization output
temp = pd.read_csv(
    "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/hdtf_binary_genes.csv",
)


#imports data as a pandas dataframe
temp = pd.DataFrame(temp)


#looks at dataframe
print(temp.head(5))


#removes decimal points and numbers after it from column names added by R
temp = temp.rename(columns=lambda x: x.split('.')[0])


#looks at dataframe
print(temp.head(5))


#removes X's from column names added by R
temp = temp.rename(columns = lambda x: x.strip('X'))

#looks at dataframe
print(temp.head(5))

#change column names to string
temp.columns = [f"{t:0>3}" for t in temp.columns]

print(temp.columns)


# Get the second column (assuming you want to add it back)
second_column = temp.iloc[:, 1]


#removes first two columns (row numbers and gene names) from dataframe
temp = temp.iloc[:, 2:]


#gets mean binary score of each gene in each cluster
temp_mean = temp.groupby(temp.columns, axis=1).mean()

def round_up(temp_mean):
    """Rounds up based on hardcoded input"""
    round_bools = temp_mean >= 0.3
    temp_mean[round_bools] = np.ceil(temp_mean[round_bools])
    not_round_bools = np.logical_not(round_bools)
    temp_mean[not_round_bools] = np.round(temp_mean[not_round_bools])
    return(temp_mean.astype(int))
    


temp_mean = round_up(temp_mean)


# Add the second column back to the resulting Series
temp_mean = pd.concat([second_column, temp_mean], axis=1)

#rename Gene column
temp_mean.columns.values[0] = "Gene"

#set Gene column top index
temp_mean.set_index("Gene", inplace=True)


temp_mean.to_csv (
    "/Users/derek/Dropbox (University of Oregon)/lab_stuffs/RNAseq/hdtf_binary_gene_avg.csv"
)



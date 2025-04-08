import glob
import os
import pandas as pd

# file directory
preprocessed_file_path = "/Users/derek/Desktop/RNAseq/preprocessed_files"

# file patterns
preprocessed_file_pattern = f"{preprocessed_file_path}/*preprocessed.csv"
preprocessed_files = glob.glob(preprocessed_file_pattern)
significant_file_pattern = f"{preprocessed_file_path}/*significant.csv"

# create significant files
for file in preprocessed_files:
    df = pd.read_csv(file)

    # Filter significant genes based on p-value adjustment
    df_significant = df[df["p_val_adj"] < 0.05]

    # Save the significant genes to a new file before further processing
    significant_output_file = file.replace("preprocessed.csv", "significant.csv")
    df_significant.to_csv(significant_output_file, index=False)

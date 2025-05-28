import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt

# Load data (adjust file path and method as needed)
filename = 'xyzeq_all_imputation_results'
df = pd.read_csv(filename+'.csv')  # or pd.read_excel("your_file.xlsx")
long_df_name = filename + '_reordered_rm2d' + '.csv'
technology_label = 'xyzeq'
rm_output_csv = technology_label+"_runtime_memory_summary.csv"

# Step 2: Identify metrics and methods dynamically
complexity_metrics =["Runtime (s)", "Memory (MB)"]
metrics = complexity_metrics
methods = ["magic", "knn", "softimpute", "simpleimpute", "scVI", "gimVI", "tangram"]
base_cols = {
    "ARI": "Base ARI",
    "NMI": "Base NMI",
    "AMI": "Base AMI",
    "HOMO": "Base HOMO",
    "zero Exp val (%)": "Raw zero Exp val (%)"
}

# Step 3: Reshape the DataFrame
long_df_list = []

for metric in metrics:
    for method in methods:
        col_name = f"{method} {metric}"
        # print(col_name)
        if col_name in df.columns:
            temp_df = df[["Dataset Name", "Size After selecting top genes", "top_genes"]].copy()
            # Split the column into two new columns
            temp_df[['cells', 'genes']] = temp_df['Size After selecting top genes'].str.split('x', expand=True).astype(int)

            temp_df["Metric"] = metric
            temp_df["Method"] = method
            temp_df["Value"] = df[col_name]
            # print(df[col_name])
            long_df_list.append(temp_df)


# Combine all
long_df = pd.concat(long_df_list, ignore_index=True)
# Drop the original column
long_df.drop(columns=['Size After selecting top genes'], inplace=True)


# Convert top_genes to sortable keys (numeric, with 'all' as last)
def top_genes_sort_key(val):
    try:
        return int(val)
    except ValueError:
        return float('inf')  # So 'all' comes last


# Add sorting keys
long_df["top_genes_str"] = long_df["top_genes"].astype(str)
long_df["top_genes_sort_key"] = long_df["top_genes_str"].map(top_genes_sort_key)

# Ensure consistent order for metrics and methods
metric_order = ['Runtime (s)', 'Memory (MB)']
method_order_baseless = ['magic', 'knn', 'softimpute', 'simpleimpute', 'scVI', 'gimVI', 'tangram']

# Step 4: Build ordered DataFrame manually using nested loops
sorted_parts = []
for tg in sorted(long_df["top_genes_str"].unique(), key=top_genes_sort_key):
    # print(tg)
    df_tg = long_df[long_df["top_genes_str"] == tg]
    for metric in metric_order:
        for method in method_order_baseless:
            sub_df = df_tg[(df_tg["Metric"] == metric) & (df_tg["Method"] == method)]
            # print(sub_df.head())
            sorted_parts.append(sub_df)

# Step 5: Concatenate all parts
sorted_df = pd.concat(sorted_parts, ignore_index=True)

# Step 6: Drop helper column
sorted_df = sorted_df.drop(columns=["top_genes_str", "top_genes_sort_key"])

# Clean up Dataset Name by removing "processed", "updated", ".h5ad", etc.
sorted_df["Dataset Name"] = (
    sorted_df["Dataset Name"]
    .str.replace(r"_processed|_updated|\.h5ad", "", regex=True)
)

sorted_df.to_csv(long_df_name, index=False)


def summarize_by_top_genes(df, technology_label, output_csv="tech_summary.csv"):
    """
    Aggregates the dataframe to get mean cells, genes, and metric values for each top_genes and Method.

    Parameters:
    - df (pd.DataFrame): Input DataFrame containing columns ['top_genes', 'cells', 'genes', 'Metric', 'Method', 'Value']
    - technology_label (str): Name of the technology (e.g., 'Visium', 'Stereo-seq') to add to result
    - output_csv (str): File name for saving the summarized results
    """

    # Pivot so metrics become columns: Runtime, Memory, etc.
    pivot_df = (
        df.groupby(["top_genes", "Method", "Metric"])
          .agg(
              mean_cells=("cells", "mean"),
              mean_genes=("genes", "mean"),
              mean_value=("Value", "mean")
          )
          .reset_index()
          .pivot(index=["top_genes", "Method", "mean_cells", "mean_genes"], 
                 columns="Metric", 
                 values="mean_value")
          .reset_index()
    )

    # Add the technology label
    pivot_df["Technology"] = technology_label

    # Reorder columns to make it clean
    cols = ["Technology", "top_genes", "Method", "mean_cells", "mean_genes"] + \
           [col for col in pivot_df.columns if col not in ["Technology", "top_genes", "Method", "mean_cells", "mean_genes"]]
    pivot_df = pivot_df[cols]

    # Save to CSV
    pivot_df.to_csv(output_csv, index=False)
    print(f"Saved summarized data to '{output_csv}'")

    return pivot_df


# Assume `sorted_df` is your full dataframe for a given technology
# For example: run for Visium
visium_summary = summarize_by_top_genes(sorted_df, technology_label=technology_label, 
                                        output_csv=rm_output_csv)


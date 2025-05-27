import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt

# Load data (adjust file path and method as needed)
filename = 'visium_all_imputation_results'
df = pd.read_csv(filename+'.csv')  # or pd.read_excel("your_file.xlsx")
long_df_name = filename + '_reordered' + '.csv'

# Step 2: Identify metrics and methods dynamically
acc_metrics = ["ARI", "NMI", "AMI", "HOMO"]
complexity_metrics =["zero Exp val (%)", "Runtime (s)", "Memory (MB)"]
metrics = acc_metrics + complexity_metrics
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
        if metric in acc_metrics:
            col_name = f"{metric}_{method}"
        elif metric in complexity_metrics:
            col_name = f"{method} {metric}"
        # print(col_name)
        if col_name in df.columns:
            temp_df = df[["Dataset Name", "top_genes"]].copy()
            temp_df["Metric"] = metric
            temp_df["Method"] = method
            temp_df["Value"] = df[col_name]
            # print(df[col_name])
            long_df_list.append(temp_df)

    # Add base metric (no method)
    if metric in base_cols:
        base_col = base_cols[metric]
        if base_col in df.columns:
            temp_df = df[["Dataset Name", "top_genes"]].copy()
            temp_df["Metric"] = metric
            temp_df["Method"] = "base"
            temp_df["Value"] = df[base_col]
            long_df_list.append(temp_df)

# Combine all
long_df = pd.concat(long_df_list, ignore_index=True)


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
metric_order = ['ARI', 'NMI', 'AMI', 'HOMO', 'zero Exp val (%)', 'Runtime (s)', 'Memory (MB)']
method_order = ['base', 'magic', 'knn', 'softimpute', 'simpleimpute', 'scVI', 'gimVI', 'tangram']
method_order_baseless = ['magic', 'knn', 'softimpute', 'simpleimpute', 'scVI', 'gimVI', 'tangram']

# Step 4: Build ordered DataFrame manually using nested loops
sorted_parts = []
for tg in sorted(long_df["top_genes_str"].unique(), key=top_genes_sort_key):
    # print(tg)
    df_tg = long_df[long_df["top_genes_str"] == tg]
    for metric in metric_order:
        for method in method_order:
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


# Function for box plot

def plot_metric_boxplot(df, metric, method_order, title="Metric Comparison", palette="Set2"):
    """
    Plots a grouped boxplot for the given metric with vertical separators and subtitles for top_genes groups.

    Parameters:
    - df: DataFrame with columns ["top_genes", "Metric", "Method", "Value"]
    - metric: str, the metric to plot (e.g., "ARI", "NMI", "AMI")
    - method_order: list, order of methods to plot (e.g., ["base", "magic", "knn", ...])
    - title: str, main plot title
    - palette: str or list, seaborn color palette
    """
    # Subset the desired metric
    metric_df = df[df["Metric"] == metric].copy()

    # Set method order
    metric_df["Method"] = pd.Categorical(metric_df["Method"], categories=method_order, ordered=True)

    # Create combined label (used for sorting, not shown)
    metric_df["top_method"] = metric_df["top_genes"].astype(str) + "\n" + metric_df["Method"].astype(str)

    # Sort for consistent plotting
    metric_df = metric_df.sort_values(["top_genes", "Method"])

    # Prepare the plot
    plt.figure(figsize=(14, 6))
    ax = sns.boxplot(data=metric_df, x="top_method", y="Value", palette=palette)

    # Set up axis ticks to show only method names
    unique_top_genes = sorted(metric_df["top_genes"].unique())
    n_methods = len(method_order)

    tick_positions = []
    tick_labels = []
    for i in range(len(unique_top_genes)):
        start = i * n_methods
        for method in method_order:
            tick_positions.append(start + method_order.index(method))
            tick_labels.append(method)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45)

    # Add vertical dashed lines to separate top_gene groups
    for i in range(1, len(unique_top_genes)):
        plt.axvline(x=i * n_methods - 0.5, color='gray', linestyle='--')

    # Annotate top_genes group titles above each group
    y_max = ax.get_ylim()[1]
    for i, tg in enumerate(unique_top_genes):
        center_x = i * n_methods + (n_methods - 1) / 2
        plt.text(center_x, y_max + 0.02 * y_max, f"{tg} genes", 
                 ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Set plot title with padding
    plt.title(title, pad=30)
    plt.ylabel(f"{metric} Score")
    plt.xlabel("")
    plt.tight_layout()
    plt.show()


# Plot ARI
plot_metric_boxplot(sorted_df, metric="ARI", method_order=method_order, title="ARI Comparison for Stereo-seq")

# Plot NMI
plot_metric_boxplot(sorted_df, metric="NMI", method_order=method_order, title="NMI Comparison for Stereo-seq")

# Plot AMI
plot_metric_boxplot(sorted_df, metric="AMI", method_order=method_order, title="AMI Comparison for Stereo-seq")

# Plot HOMO
plot_metric_boxplot(sorted_df, metric="HOMO", method_order=method_order, title="HOMO Comparison for Stereo-seq")

# Plot Runtime
plot_metric_boxplot(sorted_df, metric="Runtime (s)", method_order=method_order_baseless, title="Runtime Comparison for Stereo-seq")

# Plot Memory
plot_metric_boxplot(sorted_df, metric="Memory (MB)", method_order=method_order_baseless, title="Memory Comparison for Stereo-seq")

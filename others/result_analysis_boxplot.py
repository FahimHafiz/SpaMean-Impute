import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt


# Function for box plot
def plot_metric_boxplot(
    df, metric, method_order, labelname = "Score",
    title="Metric Comparison", palette="Set2", box_width = 0.5, fliersize =3,
    fig_width=12, fig_height=6,
    label_fontsize=14, xtick_fontsize=12, ytick_fontsize=14, title_fontsize=16, subtitle_fontsize=13
):
    """
    Plots a grouped boxplot for the given metric with vertical separators and subtitles for top_genes groups.

    Parameters:
    - df: DataFrame with columns ["top_genes", "Metric", "Method", "Value"]
    - metric: str, metric to plot (e.g., "ARI", "NMI", "AMI")
    - method_order: list, method order for box plot
    - title: str, main plot title
    - palette: seaborn color palette
    - fig_width: float, figure width
    - fig_height: float, figure height
    - label_fontsize: int, font size for y-label
    - tick_fontsize: int, font size for x-tick labels
    - title_fontsize: int, font size for title
    - subtitle_fontsize: int, font size for top_genes group subtitle
    """
    # Subset the desired metric
    metric_df = df[df["Metric"] == metric].copy()

    # Set method order
    metric_df["Method"] = pd.Categorical(metric_df["Method"], categories=method_order, ordered=True)

    # Create combined label (used for sorting, not shown)
    metric_df["top_method"] = metric_df["top_genes"].astype(str) + "\n" + metric_df["Method"].astype(str)

    # Sort for consistent plotting
    metric_df = metric_df.sort_values(["top_genes", "Method"])

    # Plotting
    plt.figure(figsize=(fig_width, fig_height))
    ax = sns.boxplot(
    data=metric_df, hue= "top_method", legend= False,
    x="top_method", y="Value",
    palette=palette,
    width=box_width,       # NEW PARAMETER
    fliersize=fliersize            # Optional, adjust outlier dot size
    )


    # X-axis ticks: show method names only
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
    # Set x-tick labels and rotate
    ax.set_xticklabels(tick_labels, rotation=45)

    # Apply separate font sizes to tick labels
    ax.tick_params(axis='x', labelsize=xtick_fontsize)
    ax.tick_params(axis='y', labelsize=ytick_fontsize)


    # Vertical lines to separate top_genes groups
    for i in range(1, len(unique_top_genes)):
        plt.axvline(x=i * n_methods - 0.5, color='gray', linestyle='--')

    # Group subtitles above each gene group
    y_max = ax.get_ylim()[1]
    for i, tg in enumerate(unique_top_genes):
        center_x = i * n_methods + (n_methods - 1) / 2
        plt.text(center_x, y_max + 0.02 * y_max, f"{tg} genes", 
                 ha='center', va='bottom', fontsize=subtitle_fontsize, fontweight='bold')

    # Labels and title
    ax.set_ylabel(f"{metric} {labelname}", fontsize=label_fontsize)
    ax.set_xlabel("")
    ax.set_title(title, fontsize=title_fontsize, pad=30)

    plt.tight_layout()
    plt.show()

### Reordering the Dataframe..................

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


### Plotting different results...................

# Plot ARI
plot_metric_boxplot(sorted_df, metric="ARI", method_order=method_order, 
                            title="ARI Comparison for Visium",  box_width = 0.5, 
                            fliersize =3, fig_width=12, fig_height=6,
                            label_fontsize=14, xtick_fontsize=11, ytick_fontsize=16, 
                            title_fontsize=16, subtitle_fontsize=13)

# Plot NMI
plot_metric_boxplot(sorted_df, metric="NMI", method_order=method_order, 
                            title="NMI Comparison for Visium",  box_width = 0.5, 
                            fliersize =3, fig_width=12, fig_height=6,
                            label_fontsize=14, xtick_fontsize=11, ytick_fontsize=16, 
                            title_fontsize=16, subtitle_fontsize=13)

# Plot AMI
plot_metric_boxplot(sorted_df, metric="AMI", method_order=method_order, 
                            title="AMI Comparison for Visium",  box_width = 0.5, 
                            fliersize =3, fig_width=12, fig_height=6,
                            label_fontsize=14, xtick_fontsize=11, ytick_fontsize=16, 
                            title_fontsize=16, subtitle_fontsize=13)


# Plot HOMO
plot_metric_boxplot(sorted_df, metric="HOMO", method_order=method_order, 
                            title="HOMO Comparison for Visium",  box_width = 0.5, 
                            fliersize =3, fig_width=12, fig_height=6,
                            label_fontsize=14, xtick_fontsize=11, ytick_fontsize=16, 
                            title_fontsize=16, subtitle_fontsize=13)

# Plot Runtime
plot_metric_boxplot(sorted_df, metric="Runtime (s)", method_order=method_order_baseless, labelname=" ",
                            title="Runtime (s) Comparison for Visium",  box_width = 0.3, 
                            fliersize =3, fig_width=12, fig_height=6,
                            label_fontsize=14, xtick_fontsize=11, ytick_fontsize=16, 
                            title_fontsize=16, subtitle_fontsize=13)

# Plot Memory
plot_metric_boxplot(sorted_df, metric="Memory (MB)", method_order=method_order_baseless, labelname=" ",
                            title="Memory (MB) Comparison for Visium",  box_width = 0.3, 
                            fliersize =3, fig_width=12, fig_height=6,
                            label_fontsize=14, xtick_fontsize=11, ytick_fontsize=16, 
                            title_fontsize=16, subtitle_fontsize=13)


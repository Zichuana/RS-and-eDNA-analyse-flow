import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import Patch

# 1. Set font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False

# 2. Load data
data = pd.read_csv('level-2-59_dem.csv', index_col=0)
# data = data.T

# 3. Extract grouping information (assuming the last column is the grouping information)
groups = data.iloc[:, -1]
print(groups)
data = data.iloc[:, :-1].T  # Remove the grouping column

# 4. Calculate relative abundance
rel_abundance = data.div(data.sum(axis=0), axis=1) * 100
print(rel_abundance)
# 5. Calculate total (or mean) abundance by group
grouped_abundance = rel_abundance.T.groupby(groups).mean().T  # Use sum() or mean()
print(grouped_abundance)
# 6. Calculate average abundance for all samples (for legend)
total_abundance = rel_abundance.mean(axis=1).sort_values(ascending=False)

# 7. Set colors
colors = plt.cm.tab20.colors
sorted_phyla = total_abundance.index
color_mapping = {phylum: colors[i % len(colors)] for i, phylum in enumerate(sorted_phyla)}

# try:
#     plt.style.use('seaborn-v0_8')
# except:
#     plt.style.use('seaborn')

# 8. Plot stacked bar chart (by group)
plt.figure(figsize=(12, 10))
bottom = np.zeros(len(grouped_abundance.columns))

for phylum in reversed(sorted_phyla):
    plt.bar(grouped_abundance.columns, grouped_abundance.loc[phylum], 
            bottom=bottom, color=color_mapping[phylum])
    bottom += grouped_abundance.loc[phylum]

# 9. Create legend (showing average abundance for all samples)
legend_handles = [Patch(facecolor=color_mapping[phylum], 
                        label=f"{phylum} ({total_abundance[phylum]:.1f}%)") 
                 for phylum in sorted_phyla]

# 10. Beautify the chart
plt.title('Relative Abundance by Altitude', fontsize=18, pad=20, fontname='Times New Roman')
plt.ylabel('Relative Abundance (%)', fontsize=18, fontname='Times New Roman')
plt.xlabel('Samples', fontsize=18, fontname='Times New Roman')
# plt.xticks(fontsize=16)
plt.xticks(rotation=45, ha='right', fontsize=16, fontname='Times New Roman')
plt.yticks(np.arange(0, 101, 10), [f"{int(y)}%" for y in np.arange(0, 101, 10)], fontsize=16, fontname='Times New Roman')

# 11. Add legend
plt.legend(handles=legend_handles,
           bbox_to_anchor=(1, 1), loc='upper left',
           fontsize=14, ncol=1, framealpha=0.5,
           handlelength=1.0, handleheight=1.0, prop={'family': 'Times New Roman', 'size': 16}) 

plt.tight_layout()
plt.savefig('community_by_group_66.png', dpi=300, bbox_inches='tight')
plt.show()
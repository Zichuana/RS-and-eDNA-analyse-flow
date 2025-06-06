import rasterio
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

plt.rcParams['font.size'] = 14
plt.rcParams['font.sans-serif'] = ['Times New Roman']

# File list
tif_files = [
    "/Users/sumyee/Downloads/last/NDVI_2017_summer.tif",
    "/Users/sumyee/Downloads/last/WET_2017_summer.tif",
    "/Users/sumyee/Downloads/last/NDBSI_2017_summer.tif",
    "/Users/sumyee/Downloads/last/LST_2017_summer.tif",
    "/Users/sumyee/Downloads/last/RSEI_2017_summer.tif",
    "/Users/sumyee/Downloads/last/SI_2017_summer.tif"
]

# Indicator names
labels = ['NDVI', 'WET', 'NDBSI', 'LST', 'RSEI', 'SI']

# X-axis ranges
xlim_ranges = [
    [-1, 1],    # NDVI
    [-1, 1],    # WET
    [-1, 1],    # NDBSI
    [0, 60],    # LST
    [0, 1],     # RSEI
    [0, 1]      # SI
]

# DataFrame to store threshold results
thresholds_df = pd.DataFrame(columns=[
    'Indicator', 
    'Min (1%)', 
    'Max (99%)', 
    'Low/Med Threshold', 
    'Med/High Threshold'
])

plt.figure(figsize=(20, 10), dpi=100)

for i, (tif_path, xlim, label) in enumerate(zip(tif_files, xlim_ranges, labels)):
    # Read data
    with rasterio.open(tif_path) as src:
        data = src.read(1).astype(np.float32)        
        data[data == src.nodata] = np.nan

    valid_data = data[~np.isnan(data)]
    
    # Calculate the 98% data range (excluding 1% lowest and 1% highest)
    lower_bound = np.percentile(valid_data, 2)
    upper_bound = np.percentile(valid_data, 98)
    filtered_data = valid_data[(valid_data >= lower_bound) & (valid_data <= upper_bound)]
    
    # Uniformly divide the 98% range into low/medium/high (by value range)
    value_range = upper_bound - lower_bound
    low_threshold = lower_bound + value_range * 1/3  # Low/Medium boundary
    high_threshold = lower_bound + value_range * 2/3  # Medium/High boundary
    
    # Store threshold results
    thresholds_df.loc[i] = [
        label,
        lower_bound,
        upper_bound,
        low_threshold,
        high_threshold
    ]
    
    # Plot histogram
    ax = plt.subplot(2, 3, i+1)
    counts, bins, _ = plt.hist(valid_data, bins=np.linspace(*xlim, 50), color='blue', alpha=0.7)
    
    # Set axis limits
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(0, np.max(counts)*1.1)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    ax.yaxis.offsetText.set_fontsize(12)
    
    # Add threshold lines
    # plt.axvline(x=low_threshold, color='orange', linestyle='--', 
    #            linewidth=1.5, label=f'Low/Med: {low_threshold:.2f}')
    # plt.axvline(x=high_threshold, color='purple', linestyle='--', 
    #            linewidth=1.5, label=f'Med/High: {high_threshold:.2f}')
    
    # Mark the 98% range boundaries (dashed gray lines)
    plt.axvline(x=lower_bound, color='orange', linestyle=':', 
               linewidth=1.5, label=f'1%: {lower_bound:.2f}')
    plt.axvline(x=upper_bound, color='purple', linestyle=':', 
               linewidth=1.5, label=f'99%: {upper_bound:.2f}')
    
    plt.axvline(x=np.min(valid_data), color='green', linestyle='--', 
               linewidth=1.5, label=f'Min: {np.min(valid_data):.2f}',)
    plt.axvline(x=np.max(valid_data), color='red', linestyle='--', 
               linewidth=1.5, label=f'Max: {np.max(valid_data):.2f}')
    
    # Title and legend
    plt.title(label, fontsize=20)
    plt.xlabel("Value", fontsize=18)
    plt.ylabel("Frequency", fontsize=18)
    plt.legend(fontsize=16)
    plt.grid(True, linestyle=':', alpha=0.5)

plt.tight_layout(pad=2.0)
plt.savefig('multi_histograms_2017_with_98percent_value_ranges.png', dpi=300, bbox_inches='tight')
plt.show()

# Save threshold results to CSV and Excel
thresholds_df.to_csv('thresholds_results.csv', index=False)
thresholds_df.to_excel('thresholds_results.xlsx', index=False)

print("Threshold results have been saved to thresholds_results.csv and thresholds_results.xlsx")
print(thresholds_df)
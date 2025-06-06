import rasterio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# plt.rcParams['font.size'] = 16
# Set the global font to Arial
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Times New Roman']  # Use Arial
plt.rcParams['axes.unicode_minus'] = False  # Fix the display of negative signs
# 1. Define file paths
tif_files = {
    'NDVI': "/Users/sumyee/Downloads/last/NDVI_2017_summer.tif",
    'WET': "/Users/sumyee/Downloads/last/WET_2017_summer.tif",
    'NDBSI': "/Users/sumyee/Downloads/last/NDBSI_2017_summer.tif",
    'LST': "/Users/sumyee/Downloads/last/LST_2017_summer.tif",
    'RSEI': "/Users/sumyee/Downloads/last/RSEI_2017_summer.tif",
    'SI': "/Users/sumyee/Downloads/last/SI_2017_summer.tif"
}

# Assume you have an elevation file
elevation_file = "Qinghai_DEM_Export_500.tif"

# 2. Read all data
data = {}
with rasterio.open(elevation_file) as src:
    elevation = src.read(1)
    mask = elevation == src.nodata  # Create a mask
    
for name, path in tif_files.items():
    with rasterio.open(path) as src:
        data[name] = src.read(1)
        # Ensure data alignment
        data[name][mask] = np.nan

# 3. Define elevation ranges
# elevation_bins = [2500, 3000, 3500, 4000, 4500, 4500, 600, 700, 800, 900, 1000, np.inf]
elevation_labels = elevation_bins = [2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000]
elevation_labels = [
    '2500-3000', '3000-3500', '3500-4000',
    '4000-4500', '4500-5000', '5000-5500',
    '5500-6000', '6000-6500', '6500-7000'
]

# 4. Classify data into elevation ranges
elevation_classes = pd.cut(elevation.flatten(), bins=elevation_bins, labels=elevation_labels)

# 5. Create DataFrame for analysis
df = pd.DataFrame()
for name in tif_files.keys():
    df[name] = data[name].flatten()
    
df['Altitude range'] = elevation_classes
df = df.dropna()  # Remove invalid values

# 6. Analyze statistics for each elevation range
grouped_stats = df.groupby('Altitude range').describe()

# Print statistics
print(grouped_stats)

# 7. Visualize the distribution for each elevation range
plt.figure(figsize=(15, 10))
for i, name in enumerate(tif_files.keys(), 1):
    plt.subplot(2, 3, i)
    sns.boxplot(x='Altitude range', y=name, data=df, showfliers=False)
    plt.title(name, fontsize=18)
    plt.xlabel('Altitude range', fontsize=16)  # Font size for x-axis label
    plt.ylabel(f'{name} index', fontsize=16)     # Font size for y-axis label
    plt.xticks(rotation=45)
plt.tight_layout()


plt.savefig('elevation.png', dpi=300)
plt.show()

# 8. Calculate the correlation matrix for each elevation range
# plt.figure(figsize=(15, 10))
# for i, elev_class in enumerate(elevation_labels[:3]):  # Only show the first 3 elevation ranges as an example
#     subset = df[df['elevation_class'] == elev_class]
#     corr = subset[tif_files.keys()].corr()
    
#     plt.subplot(1, 3, i+1)
#     sns.heatmap(corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
#     plt.title(f'Correlation - Elevation {elev_class}')
# plt.tight_layout()
# plt.savefig('xgxjz.png', dpi=300)
# plt.show()
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import LinearSegmentedColormap

# Set global parameters
plt.rcParams['font.size'] = 16
plt.rcParams['font.sans-serif'] = ['Times New Roman']

tif_files = {
    'NDVI': "/Users/sumyee/Downloads/last/NDVI_2017_summer.tif",
    'WET': "/Users/sumyee/Downloads/last/WET_2017_summer.tif",
    'NDBSI': "/Users/sumyee/Downloads/last/NDBSI_2017_summer.tif",
    'LST': "/Users/sumyee/Downloads/last/LST_2017_summer.tif",
    'RSEI': "/Users/sumyee/Downloads/last/RSEI_2017_summer.tif",
    'SI': "/Users/sumyee/Downloads/last/SI_2017_summer.tif"
}

# Define parameters for each variable
params = {
    'NDVI': {'vmin': -0.2, 'vmax': 0.5, 'cmap': LinearSegmentedColormap.from_list('arid_wet', 
            ['white', 'brown', 'yellow', 'green'])},
    'WET': {'vmin': -0.5, 'vmax': 0.3, 'cmap': LinearSegmentedColormap.from_list('arid_wet', 
            ['#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6'])},
    'NDBSI': {'vmin': -0.3, 'vmax': 0.3, 'cmap': 'bwr'},
    'LST': {'vmin': 5, 'vmax': 50, 'cmap': LinearSegmentedColormap.from_list('lst_cmap', ['blue', 'yellow', 'red'])},
    'RSEI': {'vmin': 0, 'vmax': 1, 'cmap': LinearSegmentedColormap.from_list('rsei_cmap', ['blue', 'cyan', 'green', 'yellow', 'red'])},
    'SI': {'vmin': 0, 'vmax': 0.4, 'cmap': LinearSegmentedColormap.from_list('si_cmap', ['white', 'yellow', 'red'])}
    # 'SAR_SI': {'vmin': 2, 'vmax': 10, 'cmap': LinearSegmentedColormap.from_list('sar_si_cmap', ['darkblue', 'blue', 'lightblue',
    # 'darkorange', 'red', 'darkred'])}
}

# Create a large canvas for concatenation
fig, axes = plt.subplots(2, 3, figsize=(18, 9))  # 2 rows, 3 columns, the last position is empty
axes = axes.flatten()  # Flatten the 2D array to 1D

# Iterate through each TIFF file and plot
for i, (name, file_path) in enumerate(tif_files.items()):
    print(name)
    
    # Read the TIFF file
    with rasterio.open(file_path) as src:
        data = src.read(1)
        # Handle invalid values
        data = np.ma.masked_where(data == src.nodata, data)
    
    # Get current parameters
    param = params[name]
    
    # Plot the image
    im = axes[i].imshow(data, **param)
    axes[i].set_title(name)
    axes[i].axis('off')
    
    # Add color bar
    cbar = fig.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)
    # cbar.set_label(name)

# Remove the axis of the last blank position
# fig.delaxes(axes[-1])

# Adjust layout
plt.tight_layout()
plt.savefig('combined_image_RSEI.png', dpi=300)  # Save the concatenated image

plt.show()
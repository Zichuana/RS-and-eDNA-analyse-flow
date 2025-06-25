import rasterio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd
from matplotlib.patches import Patch

plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False

def load_dem_tif(file_path):
    with rasterio.open(file_path) as src:
        dem_data = src.read(1)  # Read the first band
        profile = src.profile    # Get metadata
        bounds = src.bounds      # Get boundary coordinates
        nodata = src.nodata      # Get nodata value
        transform = src.transform # Get geotransform information
        
        # Handle nodata values
        if nodata is not None:
            dem_data = np.where(dem_data == nodata, np.nan, dem_data)
        
    return dem_data, profile, bounds, transform

def load_sample_points(csv_file):
    """Load sample point coordinates and names from a CSV file"""
    df = pd.read_csv(csv_file)
    coordinates = []
    sample_names = []
    
    for _, row in df.iterrows():
        if 'longitude' in row and 'latitude' in row:
            coordinates.append((row['longitude'], row['latitude']))
            sample_names.append(row.get('name', ''))
        elif 'lat_lon' in row:  # Compatible with combined latitude and longitude format
            parts = str(row['lat_lon']).split()
            if len(parts) == 4:
                lat = float(parts[0]) * (-1 if parts[1] == 'S' else 1)
                lon = float(parts[2]) * (-1 if parts[3] == 'W' else 1)
                coordinates.append((lon, lat))
                sample_names.append(row.get('name', ''))
    
    return coordinates, sample_names

# Elevation range colors
ELEVATION_COLORS = [
        '#F8F1E5', '#F5DEB3', '#E6C88C', '#D2B48C',
        '#BC9A6C', '#A0522D', '#8B4513', '#654321',
        '#3A2F1D'
    ]

def get_elevation_color(elevation):
    """Get the corresponding color based on the elevation value"""
    if np.isnan(elevation):
        return '#000000'  # Invalid data in black
    
    intervals = [
        (2500, 3000), (3000, 3500), (3500, 4000),
        (4000, 4500), (4500, 5000), (5000, 5500),
        (5500, 6000), (6000, 6500), (6500, 7000)
    ]
    
    for i, (low, high) in enumerate(intervals):
        if low <= elevation < high:
            return ELEVATION_COLORS[i]
    
    return '#000000'  # Out of range in black

def get_qinghai_cmap():
    # Colors extracted from GEE code
    colors = [
        '#F8F1E5', '#F5DEB3', '#E6C88C', '#D2B48C',
        '#BC9A6C', '#A0522D', '#8B4513', '#654321',
        '#3A2F1D'
    ]
    
    # Create colormap
    cmap = ListedColormap(colors)
    
    # Define color boundaries (2500-6000m divided into 10 intervals)
    bounds = list(np.linspace(2500, 7000, len(colors)+1))
    bounds = np.array(bounds)
    cmap.set_under('white') 
    return cmap, bounds

def plot_dem(dem_data, cmap, bounds, transform, coordinates, sample_names, title="Elevation map"):
    plt.figure(figsize=(12, 10))
    
    # Plot DEM
    img = plt.imshow(dem_data, cmap=cmap, vmin=2500, vmax=7000)
    
    # Remove axes
    plt.axis('off')
    
    # Add color bar - adjust position and size
    cbar = plt.colorbar(img, boundaries=bounds, ticks=bounds[:-1], 
                        spacing='proportional', shrink=0.5,
                        pad=0.02, aspect=30)
    cbar.set_label('Elevation (m)', fontsize=16)
    
    # Add coordinate points (colored by elevation)
    for (lon, lat), name in zip(coordinates, sample_names):
        try:
            # Get elevation value at the point
            row, col = rasterio.transform.rowcol(transform, lon, lat)
            elevation = dem_data[row, col]
            color = get_elevation_color(elevation)
            
            # Plot point
            plt.scatter(col, row, color=color, s=120, 
                       edgecolor='black', linewidth=1.5, zorder=5)
        except (IndexError, ValueError):
            continue
    
    plt.title(title, fontsize=18, pad=20)
    
    # Adjust layout to ensure all elements are displayed
    plt.savefig('dem_with_points66.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.tight_layout()
    plt.show()

def main():
    # Load DEM data
    tif_file = "Qinghai_DEM_Export_500.tif"
    dem_data, profile, bounds, transform = load_dem_tif(tif_file)
    
    # Load sample points from CSV file
    csv_file = "/Users/sumyee/Documents/Data/SraRunTable.csv"  # Replace with your file path
    coordinates, sample_names = load_sample_points(csv_file)
    
    # Get colormap
    cmap, color_bounds = get_qinghai_cmap()
    
    # Plot image
    plot_dem(dem_data, cmap, color_bounds, transform, 
             coordinates, sample_names)
    
    # Print information
    print(f"Data range: {bounds}")
    print(f"Resolution: {profile['transform'][0]} meters")
    print(f"Minimum elevation: {np.nanmin(dem_data):.1f} meters")
    print(f"Maximum elevation: {np.nanmax(dem_data):.1f} meters")

if __name__ == "__main__":
    main()
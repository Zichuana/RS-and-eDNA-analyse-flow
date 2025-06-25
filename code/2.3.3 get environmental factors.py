import rasterio
import numpy as np
import matplotlib.pyplot as plt
from rasterio.mask import mask
from shapely.geometry import Point, mapping
import pandas as pd
import os

plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False

def parse_lat_lon(lat_lon_str):
    """Parse a string in the format '36.26 N 94.22 E' to a (longitude, latitude) tuple"""
    if pd.isna(lat_lon_str) or not isinstance(lat_lon_str, str):
        return None
    
    try:
        parts = lat_lon_str.split()
        if len(parts) != 4:
            return None
            
        lat = float(parts[0])
        lat_dir = parts[1]
        lon = float(parts[2])
        lon_dir = parts[3]
        
        # Adjust signs based on direction
        if lat_dir == 'S':
            lat = -lat
        if lon_dir == 'W':
            lon = -lon
            
        return (lon, lat)
    except (ValueError, IndexError):
        return None

from pyproj import Transformer, CRS
from shapely.geometry import Point, mapping
from shapely.ops import transform
from math import floor

def create_precise_circle(lon, lat, radius_meters):
    """
    Create a precise circular area (considering the curvature of the Earth)
    
    Parameters:
        lon: Longitude (WGS84)
        lat: Latitude (WGS84)
        radius_meters: Radius (meters)
    
    Returns:
        Circular polygon in WGS84 coordinate system (Shapely Polygon)
    """
    # 1. Determine UTM projection
    utm_zone = floor((lon + 180) / 6) + 1
    hemisphere = 'north' if lat >= 0 else 'south'
    utm_epsg = f"EPSG:326{utm_zone}" if hemisphere == 'north' else f"EPSG:327{utm_zone}"
    utm_crs = CRS(utm_epsg)
    
    # 2. Define coordinate transformer
    wgs84_crs = CRS("EPSG:4326")
    to_utm = Transformer.from_crs(wgs84_crs, utm_crs, always_xy=True)
    to_wgs84 = Transformer.from_crs(utm_crs, wgs84_crs, always_xy=True)
    
    # 3. Generate circle in UTM
    point_utm = transform(to_utm.transform, Point(lon, lat))
    circle_utm = point_utm.buffer(radius_meters)
    
    # 4. Transform back to WGS84
    circle_wgs84 = transform(to_wgs84.transform, circle_utm)
    
    return circle_wgs84


# Used in the extract_circle_stats function:
def extract_circle_stats(tif_path, center_coords, radius=500):
    """
    Extract statistical values within a circular area (precise projection version)
    """
    with rasterio.open(tif_path) as src:
        lon, lat = center_coords
        
        # Create precise circle
        circle = create_precise_circle(lon, lat, radius)
        
        # Extract data
        try:
            out_image, out_transform = mask(src, [mapping(circle)], crop=True, nodata=np.nan)
            data = out_image[0]
            
            if np.all(np.isnan(data)):
                return None
                
            return {
                'mean': np.nanmean(data),
                'std': np.nanstd(data),
                'min': np.nanmin(data),
                'max': np.nanmax(data)
            }
        except ValueError as e:
            print(f"Error masking data: {e}")
            return None

def process_tif_files(csv_path, tif_files, output_csv, radius=500):
    """
    Main processing function
    :param csv_path: Path to input CSV file
    :param tif_files: List of TIFF files
    :param output_csv: Path to output CSV file
    :param radius: Radius of the circle (meters)
    """
    # Read CSV file
    df = pd.read_csv(csv_path)
    
    # Check if necessary columns exist
    if 'feature-id' not in df.columns or 'lat_lon' not in df.columns:
        raise ValueError("CSV file must contain 'feature-id' and 'lat_lon' columns")
    
    # Parse coordinates
    coords = []
    valid_indices = []
    for idx, row in df.iterrows():
        coord = parse_lat_lon(row['lat_lon'])
        if coord is not None:
            coords.append(coord)
            valid_indices.append(idx)
    
    # Create result DataFrame
    result_df = df.loc[valid_indices, ['feature-id']].copy()
    
    # Process each TIFF file
    for tif_path in tif_files:
        var_name = os.path.basename(tif_path).split('_')[0]
        print(f"Processing {var_name}...")
        
        # Initialize statistical columns
        for stat in ['mean', 'std']:
            result_df[f"{var_name}_{stat}"] = np.nan
        
        # Extract data for each valid coordinate
        for i, coord in enumerate(coords):
            stats = extract_circle_stats(tif_path, coord, radius)
            
            if stats is not None:
                for stat in ['mean', 'std']:
                    result_df.loc[valid_indices[i], f"{var_name}_{stat}"] = stats[stat]
    
    # Save results
    result_df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")
    
    return result_df

if __name__ == "__main__":
    # Input file path
    csv_path = "/Users/sumyee/Documents/Data/SraRunTable.csv"
    
    # List of TIFF files
    tif_files = [
        "/Users/sumyee/Downloads/NDVI_2017_summer.tif",
        "/Users/sumyee/Downloads/WET_2017_summer.tif",
        "/Users/sumyee/Downloads/NDBSI_2017_summer.tif",
        "/Users/sumyee/Downloads/LST_2017_summer.tif",
        "/Users/sumyee/Downloads/RSEI_2017_summer.tif",
        "/Users/sumyee/Downloads/SI_2017_summer.tif"
    ]
    
    # Output file path
    output_csv = "/Users/sumyee/Documents/remote/environmental_data_results_513_500m.csv"
    
    # Process data
    result_df = process_tif_files(csv_path, tif_files, output_csv, radius=500)
    
    # Display results
    print("\nProcessing complete, result preview:")
    print(result_df.head())
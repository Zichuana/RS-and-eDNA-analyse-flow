// Import shapefile and use it in GEE

// 1. Define the region of interest (replace with your SHP file)
var roi = table.geometry(); // or var roi = shp.geometry();
Map.centerObject(roi, 6);
Map.addLayer(roi, {color: 'red', fillColor: '00000000'}, 'Qinghai Province Boundary');
var startDate = '2018-06-21';
var endDate = '2018-09-22';

// 2. Load DEM data (a combination of NASADEM and SRTM)
var dem = ee.ImageCollection([
  ee.Image('NASA/NASADEM_HGT/001').select('elevation').filterDate(startDate, endDate),
  ee.Image('CGIAR/SRTM90_V4').select('elevation').filterDate(startDate, endDate)
]).mosaic().clip(roi);

// 3. Optimize the color palette (matching the characteristics of the Qinghai Plateau)
var optimizedPalette = [
  '#F5DEB3', // Light yellow-brown (elevation 2500-3000m)
  '#D2B48C', // Standard brown (3000-4000m)
  '#A0522D', // Dark brown (4000m+)
  '#8DB6CD'  // Light blue for water bodies
];

// Optimized color palette (matching the elevation gradient of the reference map)
var optimizedPalette = [
  '#F8F1E5', // 2500-2800m (very light yellow-white)
  '#F5DEB3', // 2800-3100m (light yellow-brown)
  '#E6C88C', // 3100-3400m (transition yellow-brown)
  '#D2B48C', // 3400-3700m (standard yellow soil)
  '#BC9A6C', // 3700-4000m (dark yellow soil)
  '#A0522D', // 4000-4300m (brown)
  '#8B4513', // 4300-4600m (dark brown)
  '#654321', // 4600-5000m (very dark brown)
  '#3A2F1D', // 5000m+ (near black)
  '#8DB6CD'  // Water blue (unchanged)
];

// 5. Terrain visualization parameters
var elevationVis = {
  min: 2500,  // The lowest point in Qinghai is around 1650m, start grading from 2500m
  max: 6000,  // The highest peak in the Kunlun Mountains is 6168m
  palette: optimizedPalette // Use the first 3 colors (for land areas)
};

// 7. Layer overlay display
Map.addLayer(dem, elevationVis, 'Terrain Elevation');

// Add at the end of your code:
var thumbnailURL = dem.visualize(elevationVis)
  .getThumbURL({
    dimensions: 2048,  // Maximum size
    format: 'jpg',     // Supported formats: jpg/png
    region: roi
  });

print('Click this link to download JPG:', thumbnailURL);  // Click the link in the Console to download

// 6.2 Complete data export to Google Drive (optional)
Export.image.toDrive({
  image: dem,
  description: 'Qinghai_DEM_Export',
  scale: 30,          // Export resolution (meters)
  region: roi,
  maxPixels: 1e13,
  crs: 'EPSG:4326',   // WGS84 coordinate system
  fileFormat: 'GeoTIFF'
});

print('Complete DEM data export task has been created, please run the export in the Tasks panel');
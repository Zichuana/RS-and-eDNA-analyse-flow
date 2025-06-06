// Import shapefile and use it in GEE

// Define the region of interest
var roi = table.geometry();
Map.centerObject(roi, 5);
Map.addLayer(roi, {color: 'red'}, 'Region Boundary');

// 3. Define a cloud masking function (optimized version)
function cloudMask(img) {
  var qa = img.select('StateQA');
  var cloudState = qa.bitwiseAnd(0x3);       // Check bits 0-1 (cloud state)
  var cloudShadow = qa.bitwiseAnd(1 << 2);   // Check bit 2 (cloud shadow)
  var cirrus = qa.bitwiseAnd(0x300);         // Check bits 8-9 (cirrus clouds)
  
  return cloudState.eq(0)         // Clear sky
      .and(cloudShadow.eq(0))     // No cloud shadow
      .and(cirrus.eq(0));         // No cirrus clouds
}

// 2. Set the time range (example uses 2017)
// Modify the time range section
var startDate = '2017-05-30';
var endDate = '2017-08-01';

// Load MODIS reflectance data
var modis = ee.ImageCollection('MODIS/061/MOD09A1')
  .filterDate(startDate, endDate)
  .filterBounds(roi)
  // Red, Green, NIR

// Calculate the Salt Index
var addIndices = function(image) {
  
  var scaled = image.select(['sur_refl_b01','sur_refl_b03','sur_refl_b04'])
               .multiply(0.0001)
               .updateMask(cloudMask(image));
  
  var si = scaled.expression(
    'sqrt(Green * Red)', {
      'Green': scaled.select('sur_refl_b03'),
      'Red': scaled.select('sur_refl_b01')
    }).rename('SI');
  
  // var ndsi = image.normalizedDifference(['sur_refl_b01', 'sur_refl_b04'])
  //   .rename('NDSI');
    
  return scaled.addBands([si]);
};

var withIndices = modis.map(addIndices);

var composite = withIndices.mean().clip(roi);
var composite_median = withIndices.median().clip(roi);

// Visualization
Map.centerObject(roi, 7);
Map.addLayer(composite.select('SI'), {min: -1, max: 1, palette: ['green', 'yellow', 'red']}, 'Salt Index');
Map.addLayer(composite_median.select('SI'), {min: 500, max: 3000, palette: ['green', 'yellow', 'red']}, 'Salt Index median');
// Map.addLayer(composite.select('NDSI'), {min: -0.5, max: 0.5, palette: ['blue', 'white', 'brown']}, 'NDSI');
Export.image.toDrive({
  image: composite.select('SI'),
  description: 'SI_2017_summer',
  fileNamePrefix: 'SI_2017_summer',
  scale: 500,
  region: roi,
  maxPixels: 1e13,
  folder: 'MODIS_Indices',
  crs: 'EPSG:4326'
});
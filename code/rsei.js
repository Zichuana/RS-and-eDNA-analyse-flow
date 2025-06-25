// Import shapefile and use it in GEE

// Define the research area
var roi = table.geometry();
Map.centerObject(roi, 5);
Map.addLayer(roi, {color: 'red'}, '研究区域边界');

// 2. Set the time range (the example uses 2017)

// Modify the time range section
var startDate = '2017-5-5';
var endDate = '2017-8-7';

var mod09a1 = ee.ImageCollection('MODIS/061/MOD09A1')
    .filterBounds(roi)
    .filterDate(startDate, endDate);

// Print the number of original images
print('Number of original MOD09A1 images', mod09a1.size()); 

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

// 4. Calculate ecological indices
var calculateIndices = function(img) {
  
  // 2. Band scaling and masking
  var scaled = img.select(['sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'sur_refl_b04', 
                          'sur_refl_b05', 'sur_refl_b06', 'sur_refl_b07'])
                 .multiply(0.0001)
                 .updateMask(cloudMask(img));
  
  // 3. Improved WET index
  // Using a validated MODIS WET index formula
  var wet = scaled.expression(
    '0.1147*b1 + 0.2489*b2 + 0.2408*b3 + 0.3132*b4 -0.3122*b5- 0.6416*b6 - 0.5087*b7', {
      'b1': scaled.select('sur_refl_b01'), // Red (620-670nm)
      'b2': scaled.select('sur_refl_b02'), // Near Infrared (841-876nm)
      'b3': scaled.select('sur_refl_b03'), // Blue (459-479nm)
      'b4': scaled.select('sur_refl_b04'), // Green (545-565nm)
      'b5': scaled.select('sur_refl_b05'),
      'b6': scaled.select('sur_refl_b06'), // Mid-Infrared 1 (1628-1652nm)
      'b7': scaled.select('sur_refl_b07')  // Mid-Infrared 2 (2105-2155nm)
  }).rename('WET');
  

  var evi = img.expression(
  '(NIR - RED) / (NIR + RED)', {
    'NIR': img.select('sur_refl_b02'),
    'RED': img.select('sur_refl_b01')  // L value can be adjusted (0.5 is a common default value)
  }).rename('EVI').clamp(-1, 1);

  
  var ndbsi = scaled.expression(
  '(IBI + SI) / 2', {
    'IBI': scaled.expression(
      '(2*SWIR1/(SWIR1+NIR) - (NIR/(NIR+RED) + GREEN/(GREEN+SWIR1))) / (2*SWIR1/(SWIR1+NIR) + (NIR/(NIR+RED) + GREEN/(GREEN+SWIR1)))', {
        'SWIR1': scaled.select('sur_refl_b06'),
        'NIR': scaled.select('sur_refl_b02'),
        'RED': scaled.select('sur_refl_b01'),
        'GREEN': scaled.select('sur_refl_b04')
      }),
    'SI': scaled.expression(
      '((SWIR1+RED)-(NIR+BLUE))/((SWIR1+RED)+(NIR+BLUE))', {
        'SWIR1': scaled.select('sur_refl_b06'),
        'RED': scaled.select('sur_refl_b01'),
        'NIR': scaled.select('sur_refl_b02'),
        'BLUE': scaled.select('sur_refl_b03')
      })
  }).rename('NDBSI'); 

  var lst = ee.ImageCollection('MODIS/061/MOD11A2')
    .filterDate(img.date().advance(-4, 'day'), img.date().advance(4, 'day'))
    .filterBounds(roi)
    .map(function(img) {
      // Apply quality control - using decimal or hexadecimal notation
      var qc = img.select('QC_Day');
      var quality = qc.bitwiseAnd(3).eq(0); // Use only highest quality data (3 = 0b11)
      return img.select('LST_Day_1km').updateMask(quality);
    })
    .median() // Use median instead of first()
    .multiply(0.02)
    .subtract(273.15) // Convert to Celsius
    .rename('LST');
    

  return img.addBands([evi, wet, ndbsi])
           .addBands(lst)
           .set('system:time_start', img.get('system:time_start'));
};

// Print the number of images after cloud masking
var mod09a1_cloudMasked = mod09a1.map(cloudMask);
print('Number of available images after cloud masking', mod09a1_cloudMasked.size());  // Check the number of images after cloud masking

var mod09a1_processed = mod09a1.map(calculateIndices);
print('Number of final processed images', mod09a1_processed.size());

// 5. Calculate the mean of each index (summer average)
var meanNDVI = mod09a1.map(calculateIndices).select('EVI').median().clip(roi);
var meanWET = mod09a1.map(calculateIndices).select('WET').median().clip(roi);
var meanNDBSI = mod09a1.map(calculateIndices).select('NDBSI').median().clip(roi);
var meanLST = mod09a1.map(calculateIndices).select('LST').median().clip(roi);

// 6. Set individual visualization parameters for each index
var ndviParams = {
  min: -0.2,  // NDVI theoretical range [-1,1], vegetation areas generally > 0.2
  max: 0.8,
  palette: ['white', 'brown', 'yellow', 'green']
};

// Based on the provided statistical results, set parameters
var wetParams = {
  min: -0.5,  // Min: -0.47 from your statistics
  max: 0.7,   // Max: 0.7 from your statistics
  palette: ['#d7191c', '#fdae61', '#ffffbf', '#abd9e9', '#2c7bb6'],
  opacity: 0.8
};

var ndbsiParams = {
  min: -1,    // NDBSI range [-1,1]
  max: 1,
  palette: ['blue', 'white', 'red']  // Red indicates high building/bare soil
};

// Adjust parameters based on actual statistics
var lstParams = {
  min: 10,    // Adjust based on actual statistics
  max: 45,    // Adjust based on actual statistics
  palette: ['blue', 'yellow', 'red']
};

// 7. Add to map layers (control display separately)
Map.addLayer(meanNDVI, ndviParams, 'EVI (Vegetation Index)');
Map.addLayer(meanNDBSI, ndbsiParams, 'NDBSI (Built-up/Bare Soil Index)');
Map.addLayer(meanLST, lstParams, 'LST (Land Surface Temperature)');

// Add a legend to the map
var legend = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px'
  }
});

var legendTitle = ui.Label('WET Index Legend', {fontWeight: 'bold'});
var legendItems = [
  ui.Label('> 0.3: Very Wet', {color: '0066ff'}),
  ui.Label('0-0.3: Wet', {color: '99ccff'}),
  ui.Label('-0.3-0: Dry', {color: 'ff9999'}),
  ui.Label('< -0.3: Very Dry', {color: 'ff0000'})
];

legend.add(legendTitle);
legendItems.forEach(function(item) {
  legend.add(item);
});

Map.add(legend);

// 1. Calculate the mean of each index for the summer (single image)
var summerMean = mod09a1.map(calculateIndices).select(['EVI', 'WET', 'NDBSI', 'LST']).median().clip(roi);

// Calculate the min/max values for each band
var stats = summerMean.reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: roi,
  scale: 500,  // MODIS resolution
  bestEffort: true,
  maxPixels: 1e13
});

// Print the results
print('Summer index statistics', stats);

// 2. Improved normalization function (directly process a single image)
var normalize = function(img) {
  var bandNames = ['EVI', 'WET', 'NDBSI', 'LST'];
  
  // Get the min/max values for the entire area
  var minMax = summerMean.select(bandNames).reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: roi,
    scale: 500,
    maxPixels: 1e13
  });

  // Normalize each band
  var normalizedBands = bandNames.map(function(name) {
    return img.select(name)
      .unitScale(
        ee.Number(minMax.get(name + '_min')),
        ee.Number(minMax.get(name + '_max'))
      )
      .rename(name);
  });

  return ee.Image.cat(normalizedBands);
};
// 3. Normalize the summer mean values
var normalizedMean = normalize(summerMean);
print('Summer Mean Image', summerMean);
print('Normalized Mean Image', normalizedMean);

// Calculate the min/max values for each band
var stats = normalizedMean.reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: roi,
  scale: 500,  // MODIS resolution
  bestEffort: true,
  maxPixels: 1e13
});

// Print the results
print('Statistics of normalizedMean indices', stats);

// 4. Improved PCA function (designed for a single image)
// Define a function to perform PCA and calculate contribution rates
function doPCAandCR(image) {
    // Get the band names of the image
    var bandNames = image.bandNames();
    var region = image.geometry();
    var meanDict = image.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: roi,
        scale: 500,
        tileScale: 16,
        maxPixels: 1e9
    });
    var means = ee.Image.constant(meanDict.values(bandNames));
    var centered = image.subtract(means);
    var arrays = centered.toArray();
    var covar = arrays.reduceRegion({
        reducer: ee.Reducer.centeredCovariance(),
        geometry: roi,
        scale: 500,
        tileScale: 16,
        maxPixels: 1e9
    });
    var covarArray = ee.Array(covar.get('array'));
    var eigens = covarArray.eigen();
    var eigenValues = eigens.slice(1, 0, 1);
    var eigenVectors = eigens.slice(1, 1);
    var arrayImage = arrays.toArray(1);
    var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);
    var sdImage = ee.Image(eigenValues.sqrt())
        .arrayProject([0])
        .arrayFlatten([getNewBandNames('sd', bandNames)]);
    var pcImage = principalComponents
        .arrayProject([0])
        .arrayFlatten([getNewBandNames('pc', bandNames)])
        .divide(sdImage);
    // Calculate the contribution rate of each band by dividing each eigenvalue by the sum of all eigenvalues
    var sumEigenValues = eigenValues.reduce(ee.Reducer.sum(), [0]).get([0, 0]);
    var contributionRate = eigenValues.divide(sumEigenValues).project([0]);
    
    var pc1Loadings = eigenVectors.slice(1, 0, 1).project([0]);
    var squaredLoadings = pc1Loadings.pow(2);
    var sumSquaredLoadings = squaredLoadings.reduce(ee.Reducer.sum(), [0]).get([0]);

    var bandContributions = squaredLoadings.divide(sumSquaredLoadings).project([0]);
    
    var eigenvalues = eigens.slice(1, 0, 1).project([0]); // Get the eigenvalues array
    var totalVariance = eigenvalues.reduce(ee.Reducer.sum(), [0]).get([0]); // Calculate total variance
    var contributionrate = eigenvalues.divide(totalVariance).project([0]); // Calculate the contribution rate of each principal component
    
    var contributionDict = ee.Dictionary.fromLists(
        bandNames,
        bandContributions.multiply(100).toList()
    );

    print('Eigenvalues:', eigenValues);
    print('PC1 loadings (eigenvector):', pc1Loadings);
    print('Band contributions to PC1 (%):', contributionDict);
    print('Principal components contribution rate:', contributionrate);

    return { pcImage: pcImage, contributionRate: contributionRate };
}

function getNewBandNames(prefix, bandNames) {
    var seq = ee.List.sequence(1, bandNames.length());
    return seq.map(function (b) {
        return ee.String(prefix).cat(ee.Number(b).int());
    });
}

var result = doPCAandCR(normalizedMean);
print("PCA", result.pcImage);

var visParams = {
  min: 0,
  max: 1,
  palette: ['blue', 'cyan', 'green', 'yellow', 'red']
};

// Normalize the PC1 image
var pc1Stats = result.pcImage.select(0).reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: roi,
  scale: 500,
  maxPixels: 1e13
});
print(pc1Stats)
var minPC1 = ee.Number(pc1Stats.get('pc1_min'));
var maxPC1 = ee.Number(pc1Stats.get('pc1_max'));
var normalizedPC1Image = result.pcImage.select(0).subtract(minPC1).divide(maxPC1.subtract(minPC1));

// caculate RSEI
// var rsei = ee.Image(1).subtract(normalizedPC1Image).rename('RSEI');
// var rsei = ee.Image(1).subtract(normalizedPC1Image).rename('RSEI');

Map.addLayer(normalizedPC1Image, visParams, 'resultPCA1');
print("The contribution rate of each band", result.contributionRate);


// 6. Visualization
// Map.addLayer(rsei.clip(roi), visParams, 'RSEI');

// 7. Export the result
Export.image.toDrive({
  image: normalizedPC1Image,
  description: 'RSEI_2017_summer',
  fileNamePrefix: 'RSEI_2017_summer',
  scale: 500,
  region: roi,
  maxPixels: 1e13,
  folder: 'MODIS_Indices',
  crs: 'EPSG:4326'
});

var indicesToExport = ['EVI', 'WET', 'NDBSI', 'LST'];

indicesToExport.forEach(function(index) {
  var yearlyIndex = mod09a1.map(calculateIndices).select(index).median().clip(roi);
  
  Export.image.toDrive({
    image: yearlyIndex,
    description: index + '_2017_summer',
    fileNamePrefix: index + '_2017_summer',
    scale: 500,
    region: roi,
    maxPixels: 1e13,
    folder: 'MODIS_Indices',
    crs: 'EPSG:4326'
  });
});
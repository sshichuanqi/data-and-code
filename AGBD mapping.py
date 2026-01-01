// The entire workflow was split into separate parts to avoid user memory limit exceeded or computation timed-out errors.

//Part1. data preprocessing

var geometry = ee.FeatureCollection("projects/ee-shichuanqi2002/assets/Beijing/beijing");
Map.centerObject(geometry);

//Sentinel-2 bands and DEM were used as independent variables, and GEDI L4A data was the dependent variable.
var gedi = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY');
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var glo30 = ee.ImageCollection('COPERNICUS/DEM/GLO30');

var startDate = ee.Date.fromYMD(2022, 5, 1);
//var endDate = startDate.advance(1, 'year');
var endDate = ee.Date.fromYMD(2022, 8, 31);

//Mask clouds by Cloudscore; Calculate spectral indices

var filteredS2 = s2
  .filter(ee.Filter.date(startDate, endDate))
  .filter(ee.Filter.bounds(geometry))
var s2Projection = ee.Image(filteredS2.first()).select('B4').projection();
var scaleBands = function(image) {
  return image.multiply(0.0001)
    .copyProperties(image, ['system:time_start']);
};

var csPlus = ee.ImageCollection(
    'GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var csPlusBands = csPlus.first().bandNames();

function maskLowQA(image) {
  var qaBand = 'cs';
  var clearThreshold = 0.5;
  var mask = image.select(qaBand).gte(clearThreshold);
  return image.updateMask(mask);
}

var addIndices = function(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4'])
    .rename('ndvi');

  var mndwi = image.normalizedDifference(['B3', 'B11'])
    .rename('mndwi'); //Modified Normalized Difference Water Index

  var ndbi = image.normalizedDifference(['B11', 'B8'])
    .rename('ndbi'); //Normalized Difference Built-up Index

  var evi = image.expression(
    '2.5 * ((NIR - RED)/(NIR + 6*RED - 7.5*BLUE + 1))', {
      'NIR': image.select('B8'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
    }).rename('evi'); //Enhanced Vegetation Index

  var bsi = image.expression(
      '(( X + Y ) - (A + B)) /(( X + Y ) + (A + B)) ', {
        'X': image.select('B11'),
        'Y': image.select('B4'),
        'A': image.select('B8'),
        'B': image.select('B2'),
    }).rename('bsi'); //Bare Soil Index
  
  return image
    .addBands(ndvi)
    .addBands(mndwi)
    .addBands(ndbi)
    .addBands(evi)
    .addBands(bsi);
};

var filteredS2WithCs = filteredS2.linkCollection(
    csPlus, csPlusBands);
 
var s2Processed = filteredS2WithCs
  .map(maskLowQA)
  .select('B.*')
  .map(scaleBands);
  //.map(addIndices);
  
var s2Composite = s2Processed.median()
  .setDefaultProjection(s2Projection);
Map.addLayer(s2Composite)

//Calculate slope and elevation

var glo30Filtered = glo30.filter(ee.Filter.bounds(geometry))
.select('DEM');

var demProj = glo30Filtered.first().select(0).projection();
var elevation = glo30Filtered.mosaic().rename('dem').setDefaultProjection(demProj);

var slope = ee.Terrain.slope(elevation);
var demBands = elevation.addBands(slope);

//Remove invalid or unreliable GEDI measurements

var qualityMask = function(image) {
  return image.updateMask(image.select('l4_quality_flag').eq(1))
      .updateMask(image.select('degrade_flag').eq(0));
};

// with a relative standard error > 50% 
// agbd_se / agbd > 0.5
var errorMask = function(image) {
  var relative_se = image.select('agbd_se')
    .divide(image.select('agbd'));
  return image.updateMask(relative_se.lte(0.5));
};

// Function to mask GEDI measurements on slopes > 30%
var slopeMask = function(image) {
  return image.updateMask(slope.lt(30));
};

var gediFiltered = gedi
  .filter(ee.Filter.date(startDate, endDate))
  .filter(ee.Filter.bounds(geometry));

var gediProjection = ee.Image(gediFiltered.first())
  .select('agbd').projection();

//gedi footprints before quality check ('Filtered' means shot date filtering)
var gediPre = gediFiltered.mosaic()
  .select('agbd').setDefaultProjection(gediProjection);
  
var gediProcessed = gediFiltered.map(qualityMask).map(errorMask).map(slopeMask);

var gediMosaic = gediProcessed.mosaic()
  .select('agbd').setDefaultProjection(gediProjection);

var rgbVis = {min: 0.0, max: 0.3, gamma: 1.2, bands: ['B4', 'B3', 'B2'],};
Map.addLayer(s2Composite.clip(geometry), rgbVis, 'Sentinel-2 Composite');
  
Map.addLayer(elevation.clip(geometry),{min:0, max: 1000}, 'Elevation', false);
Map.addLayer(slope.clip(geometry),{min: 0, max: 45}, 'Slope', false);
  
var gediVis = {min: 0,max: 200,palette: ['#edf8fb','#b2e2e2','#66c2a4','#2ca25f','#006d2c'],bands: ['agbd']};

Map.addLayer(gediFiltered.mosaic().clip(geometry),gediVis, 'GEDI L4A (Raw)', false);
  
Map.addLayer(gediMosaic.clip(geometry), gediVis, 'GEDI L4A (Filtered)');

//Export data to asset or to drive

var exportPath = 'projects/ee-shichuanqi2002/assets/Beijing/';

Export.image.toDrive({
  image: s2Composite.clip(geometry),
  description: 'S2_Composite_Export',
  //assetId: exportPath + 's2_composite_2022_5-8',
  region: geometry,
  scale: 10,
  maxPixels:1e13
});

Export.image.toDrive({
  image: demBands.clip(geometry),
  description: 'DEM_Bands_Export',
  //assetId: exportPath + 'dem_bands',
  region: geometry,
  scale: 30,
  maxPixels:1e13
});

Export.image.toDrive({
  image: gediMosaic.clip(geometry),
  description: 'GEDI_Mosaic_Export',
  //assetId: exportPath + 'gedi_mosaic_2022_5-8',
  region: geometry,
  scale: 25,
  maxPixels:1e13
});



//Part2. model construction and evaluation


//Import data generated in part1

var exportPath='projects/ee-shichuanqi2002/assets/Beijing/';
var s2Composite=ee.Image(exportPath+'s2_composite_2022_5-8');
var demBands=ee.Image(exportPath+'dem_bands');
var gediMosaic=ee.Image(exportPath +'gedi_mosaic_2022_5-8');

var geometry=s2Composite.geometry();
Map.centerObject(geometry);

//Resample data to common grids

var gridScale =30;
var gridProjection=ee.Projection('EPSG:3857').atScale(gridScale);
var stacked=s2Composite.addBands(demBands).addBands(gediMosaic);
var stacked=stacked.resample('bilinear')
// Aggregate pixels with 'mean' statistics
var stackedResampled=stacked
.reduceResolution({
  reducer:ee.Reducer.mean(),
  maxPixels:1024
})
.reproject({
crs:gridProjection
});
var stackedResampled=stackedResampled
.updateMask(stackedResampled.mask().gt(0));

//Use mask and statifiedsample to extract training samples because there were many null values

var predictors=s2Composite.bandNames().cat(demBands.bandNames());
var predicted=gediMosaic.bandNames().get(0);
print('predictors',predictors);
print('predicted',predicted);

var predictorImage=stackedResampled.select(predictors);
var predictedImage=stackedResampled.select([predicted]);

var classMask= predictedImage.mask().toInt().rename('class');

var numSamples =1000;

var training=stackedResampled.addBands(classMask)
.stratifiedSample({
numPoints:numSamples,
classBand:'class',
region:geometry,
scale:gridScale,
classValues:[0,1],
classPoints:[0,numSamples],
dropNulls:true,
tileScale:16,
});
print('Number of Features Extracted',training.size());
print('Sample Training Feature',training.first());


//Train a Random Forest Regression Model

var model=ee.Classifier.smileRandomForest(500)
.setOutputMode('REGRESSION')
.train({
features:training,
classProperty:predicted,
inputProperties:predictors
});
// Get model's predictions for training samples
var predicted_collection=training.classify({
  classifier:model,
  outputName:'agbd_predicted'
});

//Calculate RMSE

var calculateRmse=function(input){
  var observed=ee.Array(
    input.aggregate_array('agbd'));
  var predicted=ee.Array(
    input.aggregate_array('agbd_predicted'));
  var rmse =observed.subtract(predicted).pow(2)
    .reduce('mean',[0]).sqrt().get([0]); 
  return rmse;
};
var rmse=calculateRmse(predicted_collection);
print('RMSE',rmse)

// Create a plot of observed vs. predicted values

var chart = ui.Chart.feature.byFeature({
    features:predicted_collection.select(['agbd','agbd_predicted']),
    xProperty:'agbd',
    yProperties:['agbd_predicted'],
    }).setChartType('ScatterChart')
      .setOptions({});
print(chart);

var predictedImage = stackedResampled.classify({classifier: model,outputName: 'agbd'});

//Export data

//var exportPath='projects/ee-shichuanqi2002/assets/';
Export.image.toDrive({
image:predictedImage.clip(geometry),
description:'Predicted_Image_Export',
//assetId:exportPath+'predicted_agbd_xuyu1_2021',
region:geometry,
scale:gridScale
});

Export.image.toDrive({
image:stackedResampled.clip(geometry),
description:'stackedResampled_Export',
region:geometry,
scale:gridScale
});

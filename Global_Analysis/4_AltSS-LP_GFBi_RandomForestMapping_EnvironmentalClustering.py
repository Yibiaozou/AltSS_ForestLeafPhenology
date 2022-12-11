# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 18:31:53 2022

@author: Yibiao.Z
"""

## This code is adapted from the pipeline of Van den Hoogen et al., BioRxiv (2021). A geospatial mapping pipeline for ecologists, doi: 10.1101/2021.07.07.451145

# Import the modules of interest
import pandas as pd
import numpy as np
import subprocess
import time
import datetime
import ee
import os
from functools import partial
from pathlib import Path
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
from itertools import combinations
from itertools import repeat
from pathlib import Path

ee.Initialize()

####################################################################################################################################################################
# Configuration and project-specific settings
####################################################################################################################################################################
# Input the name of the username that serves as the home folder for asset storage
usernameFolderString = 'ybzou'

# Input the Cloud Storage Bucket that will hold the bootstrap collections when uploading them to Earth Engine
# !! This bucket should be pre-created before running this script
bucketOfInterest = 'altss_lp_gmp'

# Specify file name of raw point collection (without extension); must be a csv file; don't include '.csv'
titleOfRawPointCollection = 'GMP_Df_PCA_Kmeans_FSD'

# Input the name of the classification property: Bimodality index (BI)
classProperty = 'BI'

# Input the name of the project folder inside which all of the assets will be stored
# This folder will be generated automatically in GEE
projectFolder = 'AltssLP_GMP_PCA'

# Specify the column names where the latitude and longitude information is stored: these columns must be present in the csv containing the observations
latString = 'latitude'
longString = 'longitude'

# Name of a local folder holding input data
holdingFolder = 'D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Code/GEE_RF/Input_PCA'

# Name of a local folder for output data
outputFolder = 'D:/Zeus/ETH_zurich_MSc/ETHz_S4/MatserThesis_Crowther_Lab/Code/GEE_RF/Output_PCA'

# Create directory to hold training data
Path(outputFolder).mkdir(parents=True, exist_ok=True)

# Path to location of ee and gsutil python dependencies
bashFunction_EarthEngine = 'C:/Users/Yibiao.Z/anaconda3/Scripts/earthengine.exe'
bashFunctionGSUtil = 'C:/Users/Yibiao.Z/anaconda3/Scripts/gsutil.exe'

# Perform modeling in log space? (True or False)
log_transform_classProperty = False

# Ensemble of top 10 models from grid search? (True or False)
ensemble = False

# Proportion of variance to be covered by the PCA for interpolation/extrapolation
propOfVariance = 90

####################################################################################################################################################################
# Export settings
####################################################################################################################################################################
# Set pyramidingPolicy for exporting purposes
pyramidingPolicy = 'mean'

# Specify CRS to use (of both raw csv and final maps)
CRStoUse = 'EPSG:4326'

# Geometry to use for export
exportingGeometry = ee.Geometry.Polygon([[[-180, 88], [180, 88], [180, -88], [-180, -88]]], None, False);

# Set resolution of final image in arc seconds (30 arc seconds equals to ± 927m at the equator)
export_res = 30

# Convert resolution to degrees
res_deg = export_res/3600


####################################################################################################################################################################
# General settings
####################################################################################################################################################################

# Input the normal wait time (in seconds) for "wait and break" cells
normalWaitTime = 5

# Input a longer wait time (in seconds) for "wait and break" cells
longWaitTime = 10


####################################################################################################################################################################
# RF and Cross validation settings
####################################################################################################################################################################
# Grid search parameters; specify range
# variables per split
varsPerSplit_list = [1, 2, 4, 5, 8, 10]
# minium leaf population
leafPop_list = [1, 2, 5, 10, 20]
# Set number of trees in RF models
nTrees = 250

####################################################################################################################################################################
# Bash and Google Cloud Bucket settings
####################################################################################################################################################################
# Specify the necessary arguments to upload the files to a Cloud Storage bucket
# I.e., create bash variables in order to create/check/delete Earth Engine Assets

# Specify the arguments to these functions
arglist_preEEUploadTable = ['upload','table']
arglist_initEEUploadTable = ['--x_column', longString, '--y_column', latString, '--crs', CRStoUse]
arglist_postEEUploadTable = ['--x_column', 'Pixel_Long', '--y_column', 'Pixel_Lat']
arglist_preGSUtilUploadFile = ['cp']
formattedBucketOI = 'gs://'+bucketOfInterest
assetIDStringPrefix = '--asset_id='
arglist_CreateCollection = ['create','collection']
arglist_CreateFolder = ['create','folder']
arglist_Detect = ['asset','info']
arglist_ls = ['ls']
arglist_Delete = ['rm','-r']
stringsOfInterest = ['Asset does not exist or is not accessible']

# Compose the arguments into lists that can be run via the subprocess module
bashCommandList_Detect = [bashFunction_EarthEngine]+arglist_Detect
bashCommandList_ls = [bashFunction_EarthEngine]+arglist_ls
bashCommandList_Delete = [bashFunction_EarthEngine]+arglist_Delete
bashCommandList_CreateCollection = [bashFunction_EarthEngine]+arglist_CreateCollection
bashCommandList_CreateFolder = [bashFunction_EarthEngine]+arglist_CreateFolder

####################################################################################################################################################################
# Covariate data settings
####################################################################################################################################################################

full_composite = ee.Image("users/ybzou/AltSSLP_LabComposite/Lab_Composite_Selected")
pathOfPointCollection = holdingFolder+'/'+titleOfRawPointCollection+'.csv'
rawPointCollection = pd.read_csv(pathOfPointCollection)

	# Load rawPointCollection
    
covariateList = list(rawPointCollection.columns)


Hansen_dataset = ee.Image('UMD/hansen/global_forest_change_2020_v1_8');
GlobalTreeCover = Hansen_dataset.select('treecover2000');
forestMask = GlobalTreeCover.gte(10);

compositeToClassify = full_composite.select(covariateList).mask(forestMask)

####################################################################################################################################################################
	# Upload data to GEE
	####################################################################################################################################################################
	# Specify path of raw point collection


	# Print basic information on the csv
print('Number of observations in original Collection', rawPointCollection.shape[0])

	# Drop NAs
preppedCollection = rawPointCollection.dropna(how='any')
print('Number of observations after dropping NAs', preppedCollection.shape[0])


	# Format the bash call to upload the file to the Google Cloud Storage bucket
gsutilBashUploadList = [bashFunctionGSUtil]+arglist_preGSUtilUploadFile+[pathOfPointCollection]+[formattedBucketOI]
subprocess.run(gsutilBashUploadList)
print('GMP_Df_PCA'+' uploaded to a GCSB!')

	# Wait for a short period to ensure the command has been received by the server
time.sleep(normalWaitTime/2)

	# Wait for the GSUTIL uploading process to finish before moving on
while not all(x in subprocess.run([bashFunctionGSUtil,'ls',formattedBucketOI],stdout=subprocess.PIPE).stdout.decode('utf-8') for x in [titleOfRawPointCollection]):
	print('Not everything is uploaded...')
	time.sleep(normalWaitTime)
print('Everything is uploaded moving on...')

	# Upload the file into Earth Engine as a table asset
assetIDForSampling = 'users/'+usernameFolderString+'/'+projectFolder+'/'+titleOfRawPointCollection
earthEngineUploadTableCommands = [bashFunction_EarthEngine]+arglist_preEEUploadTable+[assetIDStringPrefix+assetIDForSampling]+[formattedBucketOI+'/'+titleOfRawPointCollection+'.csv']+arglist_initEEUploadTable
subprocess.run(earthEngineUploadTableCommands)
print('Upload to EE queued!')

	# Wait for a short period to ensure the command has been received by the server
time.sleep(normalWaitTime/2)

	# !! Break and wait
count = 1
while count >= 1:
	taskList = [str(i) for i in ee.batch.Task.list()]
	subsetList = [s for s in taskList if titleOfRawPointCollection in s]
	subsubList = [s for s in subsetList if any(xs in s for xs in ['RUNNING', 'READY'])]
	count = len(subsubList)
	print(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), 'Waiting for upload to complete...', end = '\r')
	time.sleep(normalWaitTime)
print('Upload to GEE complete! Moving on...')
    




##################################################################################################################################################################
# Hyperparameter tuning: which is done in a separate R file 'AltSS-LP_GFBi_Data_Preprocess.R'
##################################################################################################################################################################
classifierList = []

for vps in varsPerSplit_list:
	for lp in leafPop_list:

		model_name = classProperty + '_rf_VPS' + str(vps) + '_LP' + str(lp)

		rf = ee.Feature(ee.Geometry.Point([0,0])).set('cName',model_name,'c',ee.Classifier.smileRandomForest(
		numberOfTrees = nTrees,
		variablesPerSplit = vps,
		minLeafPopulation = lp,
		bagFraction = 0.632,
		seed = 42
		).setOutputMode('REGRESSION'))

		classifierList.append(rf)
        
# bestModelName = 'BI_rf_VPS5_LP1' for method 1 partition
bestModelName = 'BI_rf_VPS10_LP1'

##################################################################################################################################################################
# Classify image
##################################################################################################################################################################
fcOI = ee.FeatureCollection("users/ybzou/AltssLP_GMP_PCA/GMP_Df_PCA_Kmeans_FSD_3PC")

#covariateList_PCA = list(rawPointCollection.columns)[3:13]
#covariateList_PCA = list(rawPointCollection.columns)[2:12]
if ensemble == False:
	# Load the best model from the classifier list
	classifier = ee.Classifier(ee.Feature(ee.FeatureCollection(classifierList).filterMetadata('cName', 'equals', bestModelName).first()).get('c'))

	# Train the classifier with the collection
	trainedClassifer = classifier.train(fcOI, classProperty, covariateList)

	# Classify the image
	classifiedImage = compositeToClassify.classify(trainedClassifer,classProperty+'_Predicted')
    

if ensemble == True:
	def classifyImage(classifierName):
		# Load the best model from the classifier list
		classifier = ee.Classifier(ee.Feature(ee.FeatureCollection(classifierList).filterMetadata('cName', 'equals', classifierName).first()).get('c'))

		# Train the classifier with the collection
		trainedClassifer = classifier.train(fcOI, classProperty, covariateList)

		# Classify the image
		classifiedImage = compositeToClassify.classify(trainedClassifer,classProperty+'_Predicted')

		return classifiedImage

# export the predicted map by the best tuned model
exportTask = ee.batch.Export.image.toAsset(
	image = classifiedImage.toFloat(),
	description = classProperty+'_PredictedEnsembleMean_FSD',
	assetId = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_PredictedEnsembleMean_PCA_Kmeans_FSD_3PC' ,
	crs = CRStoUse,
	crsTransform = '['+str(res_deg)+',0,-180,0,'+str(-res_deg)+',90]',
	region = exportingGeometry,
	maxPixels = int(1e13),
	pyramidingPolicy = {".default": pyramidingPolicy}
);
exportTask.start()
print('Image export task started, moving on')


##################################################################################################################################################################
# Multivariate (PCA) int-ext analysis
##################################################################################################################################################################
# PCA interpolation/extrapolation helper function
def assessExtrapolation(fcOfInterest, propOfVariance):
	# Compute the mean and standard deviation of each band, then standardize the point data
# =============================================================================
#     fcOfInterest = preppedCollection[covariateList]
#     propOfVariance = propOfVariance
# =============================================================================
    
	meanVector = fcOfInterest.mean()
	stdVector = fcOfInterest.std()
	standardizedData = (fcOfInterest-meanVector)/stdVector

	# Then standardize the composite from which the points were sampled
	meanList = meanVector.tolist()
	stdList = stdVector.tolist()
	bandNames = list(meanVector.index)
	meanImage = ee.Image(meanList).rename(bandNames)
	stdImage = ee.Image(stdList).rename(bandNames)
	standardizedImage = compositeToClassify.subtract(meanImage).divide(stdImage)

	# Run a PCA on the point samples
	pcaOutput = PCA()
	pcaOutput.fit(standardizedData)

	# Save the cumulative variance represented by each PC
	cumulativeVariance = np.cumsum(np.round(pcaOutput.explained_variance_ratio_, decimals=4)*100)

	# Make a list of PC names for future organizational purposes
	pcNames = ['PC'+str(x) for x in range(1,fcOfInterest.shape[1]+1)]

	# Get the PC loadings as a data frame
	loadingsDF = pd.DataFrame(pcaOutput.components_,columns=[str(x)+'_Loads' for x in bandNames],index=pcNames)

	# Get the original data transformed into PC space
	transformedData = pd.DataFrame(pcaOutput.fit_transform(standardizedData,standardizedData),columns=pcNames)

	# Make principal components images, multiplying the standardized image by each of the eigenvectors
	# Collect each one of the images in a single image collection

	# First step: make an image collection wherein each image is a PC loadings image
	listOfLoadings = ee.List(loadingsDF.values.tolist())
	eePCNames = ee.List(pcNames)
	zippedList = eePCNames.zip(listOfLoadings)
	def makeLoadingsImage(zippedValue):
		return ee.Image.constant(ee.List(zippedValue).get(1)).rename(bandNames).set('PC',ee.List(zippedValue).get(0))
	loadingsImageCollection = ee.ImageCollection(zippedList.map(makeLoadingsImage))

	# Second step: multiply each of the loadings image by the standardized image and reduce it using a "sum"
	# to finalize the matrix multiplication
	def finalizePCImages(loadingsImage):
		PCName = ee.String(ee.Image(loadingsImage).get('PC'))
		return ee.Image(loadingsImage).multiply(standardizedImage).reduce('sum').rename([PCName]).set('PC',PCName)
	principalComponentsImages = loadingsImageCollection.map(finalizePCImages)

	# Choose how many principal components are of interest in this analysis based on amount of
	# variance explained
	numberOfComponents = sum(i < propOfVariance for i in cumulativeVariance)+1
	print('Number of Principal Components being used:',numberOfComponents)

	# Compute the combinations of the principal components being used to compute the 2-D convex hulls
	tupleCombinations = list(combinations(list(pcNames[0:numberOfComponents]),2))
	print('Number of Combinations being used:',len(tupleCombinations))

	# Generate convex hulls for an example of the principal components of interest
	cHullCoordsList = list()
	for c in tupleCombinations:
		firstPC = c[0]
		secondPC = c[1]
		outputCHull = ConvexHull(transformedData[[firstPC,secondPC]])
		listOfCoordinates = transformedData.loc[outputCHull.vertices][[firstPC,secondPC]].values.tolist()
		flattenedList = [val for sublist in listOfCoordinates for val in sublist]
		cHullCoordsList.append(flattenedList)

	# Reformat the image collection to an image with band names that can be selected programmatically
	pcImage = principalComponentsImages.toBands().rename(pcNames)

	# Generate an image collection with each PC selected with it's matching PC
	listOfPCs = ee.List(tupleCombinations)
	listOfCHullCoords = ee.List(cHullCoordsList)
	zippedListPCsAndCHulls = listOfPCs.zip(listOfCHullCoords)

	def makeToClassifyImages(zippedListPCsAndCHulls):
		imageToClassify = pcImage.select(ee.List(zippedListPCsAndCHulls).get(0)).set('CHullCoords',ee.List(zippedListPCsAndCHulls).get(1))
		classifiedImage = imageToClassify.rename('u','v').classify(ee.Classifier.spectralRegion([imageToClassify.get('CHullCoords')]))
		return classifiedImage

	classifedImages = ee.ImageCollection(zippedListPCsAndCHulls.map(makeToClassifyImages))
	finalImageToExport = classifedImages.sum().divide(ee.Image.constant(len(tupleCombinations)))

	return finalImageToExport

# PCA interpolation-extrapolation image
PCA_int_ext = assessExtrapolation(rawPointCollection[covariateList], propOfVariance).rename('PCA_pct_int_ext')
exportTask = ee.batch.Export.image.toAsset(
	image = PCA_int_ext.toFloat(),
	description = classProperty+'_PCA_pct_int_ext_FSD',
	assetId = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_PCA_pct_int_ext_Kmeans_FSD',
	crs = CRStoUse,
	crsTransform = '['+str(res_deg)+',0,-180,0,'+str(-res_deg)+',90]',
	region = exportingGeometry,
	maxPixels = int(1e13),
	pyramidingPolicy = {".default": pyramidingPolicy}
);
exportTask.start()
print('Image export task started, moving on')
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 18:09:03 2022

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

# ee.Authenticate(auth_mode='paste')
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
titleOfRawPointCollection = 'GMP_Df_FSD'

# Input the name of the classification property: Bimodality index (BI)
classProperty = 'BI'

# Input the name of the project folder inside which all of the assets will be stored
# This folder will be generated automatically in GEE
projectFolder = 'AltssLP_GMP_FSD'

# Specify the column names where the latitude and longitude information is stored: these columns must be present in the csv containing the observations
latString = 'latitude'
longString = 'longitude'

# Name of a local folder holding input data
holdingFolder = '**/Global_Analysis/Data/GEE_RF/Input_FSD'

# Name of a local folder for output data
outputFolder = '**/Global_Analysis/Data/GEE_RF/Output_FSD'

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
# Covariate data settings
####################################################################################################################################################################


full_composite = ee.Image("users/ybzou/AltSSLP_LabComposite/Lab_Composite_Selected")


pathOfPointCollection = holdingFolder+'/'+titleOfRawPointCollection+'.csv'
rawPointCollection = pd.read_csv(pathOfPointCollection)

	# Load rawPointCollection
    
covariateList = list(rawPointCollection.columns)


Hansen_dataset = ee.Image('UMD/hansen/global_forest_change_2020_v1_8');
GlobalTreeCover = Hansen_dataset.select('treecover2000');
forestMask = GlobalTreeCover.gte(10); # get global forest regions with tree cover >= 10%

compositeToClassify = full_composite.select(covariateList).mask(forestMask) # mask all covariates layers with the global forest regions


# Scale of composite
scale = full_composite.projection().nominalScale().getInfo()

####################################################################################################################################################################
# Additional settings
####################################################################################################################################################################

####################################################################################################################################################################
# RF and Cross validation settings
####################################################################################################################################################################
# Grid search parameters; specify range
# variables per split
varsPerSplit_list = [1, 2, 4, 5, 8, 10, 15, 20]

# minium leaf population
leafPop_list = [1, 2, 5, 10, 20]

## model number
num_model = len(varsPerSplit_list)*len(leafPop_list)
if (num_model > 10): num_model=10

# Set k for k-fold CV and make a list of the k-fold CV assignments to use
k = 10
kList = list(range(1,k+1))

# Metric to use for sorting k-fold CV hyperparameter tuning
sort_acc_prop = 'Mean_R2' # (either one of 'Mean_R2', 'Mean_MAE', 'Mean_RMSE')

# Set number of trees in RF models
nTrees = 250

# Input the name of the property that holds the CV fold assignment
cvFoldString = 'CV_Fold'

# Input the title of the CSV that will hold all of the data that has been given a CV fold assignment
titleOfCSVWithCVAssignments = classProperty+"_training_data"

assetIDForCVAssignedColl = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_training_data'

####################################################################################################################################################################
# Bootstrap settings
####################################################################################################################################################################
# Number of bootstrap iterations
bootstrapIterations = 100

# Generate the seeds for bootstrapping
seedsToUseForBootstrapping = list(range(1, bootstrapIterations+1))

# Input the header text that will name the bootstrapped dataset
bootstrapSamples = classProperty+'_bootstrapSamples'

# Stratification inputs
# Write the name of the variable used for stratification
# !! This variable should be included in the input dataset
stratificationVariableString = "WWF_Biome"

# Input the dictionary of values for each of the stratification category levels
# !! This area breakdown determines the proportion of each biome to include in every bootstrap
strataDict = {
	1: 14.900835665820974,
	2: 2.941697660221864,
	3: 0.526059731441294,
	4: 9.56387696566245,
	5: 2.865354077500338,
	6: 11.519674266872787,
	7: 16.26999434439293,
	8: 8.047078485979089,
	9: 0.861212221078014,
	10: 3.623974712557433,
	11: 6.063922959332467,
	12: 2.5132866428302836,
	13: 20.037841544639985,
	98: 0.0026519072167008,
}

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
# Helper functions
####################################################################################################################################################################
# Function to convert FeatureCollection to Image
def fcToImg(f):
	# Reduce to image, take mean per pixel
   img = sampledFC.reduceToImage(
		properties = [f],
		reducer = ee.Reducer.mean()
	)
   return img

# Function to convert GEE FC to pd.DataFrame. Not ideal as it's calling .getInfo(), but does the job
def GEE_FC_to_pd(fc):
	result = []
	# Fetch data as a list
	values = fc.toList(100000).getInfo()
	# Fetch column names
	BANDS = fc.first().propertyNames().getInfo()
	# Remove system:index if present
	if 'system:index' in BANDS: BANDS.remove('system:index')

	# Convert to data frame
	for item in values:
		values = item['properties']
		row = [str(values[key]) for key in BANDS]
		row = ",".join(row)
		result.append(row)

	df = pd.DataFrame([item.split(",") for item in result], columns = BANDS)
	df.replace('None', np.nan, inplace = True)

	return df

# R^2 function
def coefficientOfDetermination(fcOI,propertyOfInterest,propertyOfInterest_Predicted):
	# Compute the mean of the property of interest
	propertyOfInterestMean = ee.Number(ee.Dictionary(ee.FeatureCollection(fcOI).select([propertyOfInterest]).reduceColumns(ee.Reducer.mean(),[propertyOfInterest])).get('mean'));

	# Compute the total sum of squares
	def totalSoSFunction(f):
		return f.set('Difference_Squared',ee.Number(ee.Feature(f).get(propertyOfInterest)).subtract(propertyOfInterestMean).pow(ee.Number(2)))
	totalSumOfSquares = ee.Number(ee.Dictionary(ee.FeatureCollection(fcOI).map(totalSoSFunction).select(['Difference_Squared']).reduceColumns(ee.Reducer.sum(),['Difference_Squared'])).get('sum'))

	# Compute the residual sum of squares
	def residualSoSFunction(f):
		return f.set('Residual_Squared',ee.Number(ee.Feature(f).get(propertyOfInterest)).subtract(ee.Number(ee.Feature(f).get(propertyOfInterest_Predicted))).pow(ee.Number(2)))
	residualSumOfSquares = ee.Number(ee.Dictionary(ee.FeatureCollection(fcOI).map(residualSoSFunction).select(['Residual_Squared']).reduceColumns(ee.Reducer.sum(),['Residual_Squared'])).get('sum'))

	# Finalize the calculation
	r2 = ee.Number(1).subtract(residualSumOfSquares.divide(totalSumOfSquares))

	return ee.Number(r2)

# RMSE function
def RMSE(fcOI,propertyOfInterest,propertyOfInterest_Predicted):
	# Compute the squared difference between observed and predicted
	def propDiff(f):
		diff = ee.Number(f.get(propertyOfInterest)).subtract(ee.Number(f.get(propertyOfInterest_Predicted)))

		return f.set('diff', diff.pow(2))

	# calculate RMSE from squared difference
	rmse = ee.Number(fcOI.map(propDiff).reduceColumns(ee.Reducer.mean(), ['diff']).get('mean')).sqrt()

	return rmse

# MAE function
def MAE(fcOI,propertyOfInterest,propertyOfInterest_Predicted):
	# Compute the absolute difference between observed and predicted
	def propDiff(f):
		diff = ee.Number(f.get(propertyOfInterest)).subtract(ee.Number(f.get(propertyOfInterest_Predicted)))

		return f.set('diff', diff.abs())

	# calculate RMSE from squared difference
	mae = ee.Number(fcOI.map(propDiff).reduceColumns(ee.Reducer.mean(), ['diff']).get('mean'))

	return mae

# Function to add folds stratified per biome
def assignFolds(biome):
	fc_filtered = fc_agg.filter(ee.Filter.eq(stratificationVariableString, biome))

	cvFoldsToAssign = ee.List.sequence(0, fc_filtered.size()).map(lambda i: ee.Number(i).mod(k).add(1))

	fc_sorted = fc_filtered.randomColumn(seed = biome).sort('random')

	fc_wCVfolds = ee.FeatureCollection(cvFoldsToAssign.zip(fc_sorted.toList(fc_filtered.size())).map(lambda f: ee.Feature(ee.List(f).get(1)).set(cvFoldString, ee.List(f).get(0))))

	return fc_wCVfolds


# Define a function to take a feature with a classifier of interest
def computeCVAccuracyAndRMSE(featureWithClassifier):
	# Pull the classifier from the feature
	cOI = ee.Classifier(featureWithClassifier.get('c'))

	# Create a function to map through the fold assignments and compute the overall accuracy
	# for all validation folds
	def computeAccuracyForFold(foldFeature):
		# Organize the training and validation data
		foldNumber = ee.Number(ee.Feature(foldFeature).get('Fold'))
		trainingData = fcOI.filterMetadata(cvFoldString,'not_equals',foldNumber)
		validationData = fcOI.filterMetadata(cvFoldString,'equals',foldNumber)

		# Train the classifier and classify the validation dataset
		trainedClassifier = cOI.train(trainingData,classProperty,covariateList)
		outputtedPropName = classProperty+'_Predicted'
		classifiedValidationData = validationData.classify(trainedClassifier,outputtedPropName)

		# Compute accuracy metrics
		r2ToSet = coefficientOfDetermination(classifiedValidationData,classProperty,outputtedPropName)
		rmseToSet = RMSE(classifiedValidationData,classProperty,outputtedPropName)
		maeToSet = MAE(classifiedValidationData,classProperty,outputtedPropName)
		return foldFeature.set('R2',r2ToSet).set('RMSE', rmseToSet).set('MAE', maeToSet)

	# Compute the mean and std dev of the accuracy values of the classifier across all folds
	accuracyFC = kFoldAssignmentFC.map(computeAccuracyForFold)
	meanAccuracy = accuracyFC.aggregate_mean('R2')
	sdAccuracy = accuracyFC.aggregate_total_sd('R2')

	# Calculate mean and std dev of RMSE
	RMSEvals = accuracyFC.aggregate_array('RMSE')
	RMSEvalsSquared = RMSEvals.map(lambda f: ee.Number(f).multiply(f))
	sumOfRMSEvalsSquared = RMSEvalsSquared.reduce(ee.Reducer.sum())
	meanRMSE = ee.Number.sqrt(ee.Number(sumOfRMSEvalsSquared).divide(k))

	sdRMSE = accuracyFC.aggregate_total_sd('RMSE')

	# Calculate mean and std dev of MAE
	meanMAE = accuracyFC.aggregate_mean('MAE')
	sdMAE= accuracyFC.aggregate_total_sd('MAE')

	# Compute the feature to return
	featureToReturn = featureWithClassifier.select(['cName']).set('Mean_R2',meanAccuracy,'StDev_R2',sdAccuracy, 'Mean_RMSE',meanRMSE,'StDev_RMSE',sdRMSE, 'Mean_MAE',meanMAE,'StDev_MAE',sdMAE)
	return featureToReturn

####################################################################################################################################################################
# Start of pipeline
####################################################################################################################################################################
# Check if FC is already present, if not perform data aggregation and upload

files_present = subprocess.run(bashCommandList_ls+['users/'+usernameFolderString+'/'+projectFolder], stdout=subprocess.PIPE).stdout.decode('utf-8').split('\n')
files_present = [fileName.replace('projects/earthengine-legacy/assets/', '') for fileName in files_present]

if (assetIDForCVAssignedColl in files_present):
	# Load the collection with the assigned K-Fold assignments
	fcOI = ee.FeatureCollection(assetIDForCVAssignedColl)

	print('Size of training data:', fcOI.size().getInfo())
	print('Data prep already performed. Moving on with grid search')

else:
	####################################################################################################################################################################
	# Initialization
	####################################################################################################################################################################

	# Turn the folder string into an assetID and perform the folder creation
	assetIDToCreate_Folder = 'users/'+usernameFolderString+'/'+projectFolder
	print(assetIDToCreate_Folder,'being created...')


	# Create the folder within Earth Engine
	subprocess.run(bashCommandList_CreateFolder+[assetIDToCreate_Folder])
	while any(x in subprocess.run(bashCommandList_Detect+[assetIDToCreate_Folder],stdout=subprocess.PIPE).stdout.decode('utf-8') for x in stringsOfInterest):
		print('Waiting for asset to be created...')
		time.sleep(normalWaitTime)
	print('Asset created!')

	# Sleep to allow the server time to receive incoming requests
	time.sleep(normalWaitTime/2)


	####################################################################################################################################################################
	# Upload data to GEE
	####################################################################################################################################################################
	# Specify path of raw point collection

	pathOfPointCollection = holdingFolder+'/'+titleOfRawPointCollection+'.csv'

	# Load rawPointCollection
	rawPointCollection = pd.read_csv(pathOfPointCollection)[[latString, longString, classProperty]]

	# Print basic information on the csv
	print('Number of observations in original Collection', rawPointCollection.shape[0])

	# Drop NAs
	preppedCollection = pd.read_csv(pathOfPointCollection).dropna(how='any')
	print('Number of observations after dropping NAs', preppedCollection.shape[0])

	# Write the CSV to disk and upload it to Earth Engine as a Feature Collection
	localPathToCVAssignedData = holdingFolder+'/'+titleOfCSVWithCVAssignments+'.csv'
	preppedCollection.to_csv(localPathToCVAssignedData,index=False)

	# Format the bash call to upload the file to the Google Cloud Storage bucket
	gsutilBashUploadList = [bashFunctionGSUtil]+arglist_preGSUtilUploadFile+[pathOfPointCollection]+[formattedBucketOI]
	subprocess.run(gsutilBashUploadList)
	print(titleOfCSVWithCVAssignments+' uploaded to a GCSB!')

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
    

	####################################################################################################################################################################
	# Sample dataset, pixel aggregation, assign cross validation folds
	####################################################################################################################################################################
	# Raw dataset as FeatureCollection
	pointCollection = ee.FeatureCollection('users/'+usernameFolderString+'/'+projectFolder+'/'+titleOfRawPointCollection)

    # To my understanding, sampledFC can be replace as: ee.FeatureCollection("self-processed composite")
    
	# Sample composite
	sampledFC = ee.FeatureCollection('users/ybzou/AltssLP_GMP_FSD/GMP_Df_FSD')
    
    full_composite.reduceRegions(
		reducer = ee.Reducer.first(),
		collection = pointCollection,
		scale = scale,
		tileScale = 16)
    
	def stringToNumeric(feature):
        num = ee.Number.parse(sampledFC.get(feature))
        return sampledFC.set(feature, num)
    
	# Instantiate list of properties to select
    fullPropList = covariateList + [stratificationVariableString] + [classProperty]


	# Perform pixel aggregation by converting FeatureCollection to Image
	fc_asImg = ee.ImageCollection(list(map(fcToImg, fullPropList))).toBands().rename(fullPropList)

	# And sampling that image to get pixel values
	fc_agg = fc_asImg.sample(
		region = exportingGeometry,
		scale = scale,
		projection = 'EPSG:4326',
		geometries = True,
	)

	# Log transform classProperty, if specified
	if log_transform_classProperty == True:
		fc_agg = fc_agg.filter(ee.Filter.gt(classProperty, 0)).map(lambda f: f.set(classProperty, ee.Number(f.get(classProperty)).log()))

	# Convert biome column to int
	fc_agg = fc_agg.map(lambda f: f.set(stratificationVariableString, ee.Number(f.get(stratificationVariableString)).toInt()))

	# Retrieve biome classes present in dataset
	biome_list = fc_agg.aggregate_array(stratificationVariableString).distinct()

	# Assign folds to each feature, stratified by biome
	fcToExport = ee.FeatureCollection(biome_list.map(assignFolds)).flatten()

	# Export to assets
	fcOI_exportTask = ee.batch.Export.table.toAsset(
		collection = fcToExport,
		description = classProperty+'_training_data_export',
		assetId = assetIDForCVAssignedColl
	)
	fcOI_exportTask.start()

	# Wait for a short period to ensure the command has been received by the server
	time.sleep(normalWaitTime/2)

	# !! Break and wait
	count = 1
	while count >= 1:
		taskList = [str(i) for i in ee.batch.Task.list()]
		subsetList = [s for s in taskList if classProperty in s]
		subsubList = [s for s in subsetList if any(xs in s for xs in ['RUNNING', 'READY'])]
		count = len(subsubList)
		print(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), 'Number of running jobs:', count)
		time.sleep(normalWaitTime)
	print('Moving on...')

	# Load the collection with the pre-assigned K-Fold assignments
    fcOI = ee.FeatureCollection(assetIDForCVAssignedColl)

	# Write to file
	GEE_FC_to_pd(fcOI).to_csv('**/Global_Analysis/Data/GEE_RF/Output_FSD/'+classProperty+'_training_data.csv')

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

	# Get top model name from hyperparameter tuning
bestModelName = 'BI_rf_VPS20_LP2'
grid_search_results.limit(1, sort_acc_prop, False).first().get('cName')

	# Get top 10 models from hyperparameter tuning
top_10Models = ee.List(['BI_rf_VPS20_LP5', 'BI_rf_VPS20_LP2', 'BI_rf_VPS15_LP2', 'BI_rf_VPS15_LP5', 'BI_rf_VPS20_LP1',
                'BI_rf_VPS20_LP10', 'BI_rf_VPS15_LP1', 'BI_rf_VPS15_LP10', 'BI_rf_VPS10_LP2', 'BI_rf_VPS10_LP1'])
grid_search_results.limit(10, sort_acc_prop, False).aggregate_array('cName')


##################################################################################################################################################################
# Classify image
##################################################################################################################################################################
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

	# Classify the images, return mean
	classifiedImage = ee.ImageCollection(top_10Models.map(classifyImage)).mean().rename(classProperty+'_PredictedEnsembleMean')

# export the predicted map by the best tuned model
exportTask = ee.batch.Export.image.toAsset(
	image = classifiedImage.toFloat(),
	description = classProperty+'_PredictedEnsembleMean',
	assetId = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_PredictedEnsembleMean' ,
	crs = CRStoUse,
	crsTransform = '['+str(res_deg)+',0,-180,0,'+str(-res_deg)+',90]',
	region = exportingGeometry,
	maxPixels = int(1e13),
	pyramidingPolicy = {".default": pyramidingPolicy}
);
exportTask.start()
print('Image export task started, moving on')


##################################################################################################################################################################
# Bootstrapping
##################################################################################################################################################################

# Input the number of points to use for each bootstrap model: equal to number of observations in training dataset
preppedCollection = pd.read_csv(outputFolder+'/'+classProperty+'_training_data.csv')
bootstrapModelSize = preppedCollection.shape[0]

# Create an empty dataframe
stratSample = preppedCollection.head(0)

# Bootstrap sampling, add to df
for n in seedsToUseForBootstrapping:
	# Perform the subsetting
	# stratSample = preppedCollection.groupby(stratificationVariableString, group_keys=False).apply(lambda x: x.sample(n=int(round((strataDict.get(x.name)/100)*bootstrapModelSize)), replace=True, random_state=n))
	sampleToConcat = preppedCollection.groupby(stratificationVariableString, group_keys=False).apply(lambda x: x.sample(n=int(round((strataDict.get(x.name)/100)*bootstrapModelSize)), replace=True, random_state=n))
	sampleToConcat['bootstrapIteration'] = n
	stratSample = pd.concat([stratSample, sampleToConcat])

# Format the title of the CSV and export it to a holding location
fullLocalPath = holdingFolder+'/'+bootstrapSamples+'.csv'
stratSample.to_csv(holdingFolder+'/'+bootstrapSamples+'.csv',index=False)

# Format the bash call to upload the files to the Google Cloud Storage bucket
gsutilBashUploadList = [bashFunctionGSUtil]+arglist_preGSUtilUploadFile+[fullLocalPath]+[formattedBucketOI]
subprocess.run(gsutilBashUploadList)
print(bootstrapSamples+' uploaded to a GCSB!')

# Wait for the GSUTIL uploading process to finish before moving on
while not all(x in subprocess.run([bashFunctionGSUtil,'ls',formattedBucketOI],stdout=subprocess.PIPE).stdout.decode('utf-8') for x in [bootstrapSamples]):
	print('Not everything is uploaded...')
	time.sleep(5)
print('Everything is uploaded moving on...')

# Upload the file into Earth Engine as a table asset
assetIDForCVAssignedColl = 'users/'+usernameFolderString+'/'+projectFolder+'/'+bootstrapSamples 
earthEngineUploadTableCommands = [bashFunction_EarthEngine]+arglist_preEEUploadTable+[assetIDStringPrefix+assetIDForCVAssignedColl]+[formattedBucketOI+'/'+bootstrapSamples+'.csv']+arglist_postEEUploadTable
subprocess.run(earthEngineUploadTableCommands)
print('Upload to EE queued!')

# Wait for a short period to ensure the command has been received by the server
time.sleep(normalWaitTime/2)

# !! Break and wait
count = 1
while count >= 1:
	taskList = [str(i) for i in ee.batch.Task.list()]
	subsetList = [s for s in taskList if classProperty in s]
	subsubList = [s for s in subsetList if any(xs in s for xs in ['RUNNING', 'READY'])]
	count = len(subsubList)
	print(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), 'Number of running jobs:', count)
	time.sleep(normalWaitTime)
print('Moving on...')

# Load the best model from the classifier list
classifierToBootstrap = ee.Classifier(ee.Feature(ee.FeatureCollection(classifierList).filterMetadata('cName','equals',bestModelName).first()).get('c'))

# Create empty list to store all fcs
fcList = []
# Run a for loop to create multiple bootstrap iterations
for n in seedsToUseForBootstrapping:

	# Format the title of the CSV and export it to a holding location
	# titleOfColl = bootstrapSamples+str(n).zfill(3)
	collectionPath = 'users/'+usernameFolderString+'/'+projectFolder+'/'+bootstrapSamples

	# Load the collection from the path
	fcToTrain = ee.FeatureCollection(collectionPath).filter(ee.Filter.eq('bootstrapIteration', n))

	# Append fc to list
	fcList.append(fcToTrain)

# Helper fucntion to train a RF classifier and classify the composite image
def bootstrapFunc(fc):
	# Train the classifier with the collection
	trainedClassifer = classifierToBootstrap.train(fc,classProperty,covariateList)

	# Classify the image
	classifiedImage = compositeToClassify.classify(trainedClassifer,classProperty+'_Predicted')

	return classifiedImage

# Reduce bootstrap images to mean
meanImage = ee.ImageCollection.fromImages(list(map(bootstrapFunc, fcList))).reduce(
	reducer = ee.Reducer.mean()
)

# export bootstrapping mean
exportTask = ee.batch.Export.image.toAsset(
	image = meanImage.toFloat(),
	description = classProperty+'_PredictedBootstrapMean',
	assetId = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_PredictedBootstrapMean' ,
	crs = CRStoUse,
	crsTransform = '['+str(res_deg)+',0,-180,0,'+str(-res_deg)+',90]',
	region = exportingGeometry,
	maxPixels = int(1e13),
	pyramidingPolicy = {".default": pyramidingPolicy}
);
exportTask.start()
print('Image export task started, moving on')

# Reduce bootstrap images to lower and upper CIs
upperLowerCIImage = ee.ImageCollection.fromImages(list(map(bootstrapFunc, fcList))).reduce(
	reducer = ee.Reducer.percentile([2.5,97.5],['lower','upper'])
)

# export bootstrapping upper and lower 95% CI bands
exportTask = ee.batch.Export.image.toAsset(
	image = upperLowerCIImage.toFloat(),
	description = classProperty+'_PredictedBootstrapUpperLowerCI',
	assetId = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_PredictedBootstrapUpperLowerCI' ,
	crs = CRStoUse,
	crsTransform = '['+str(res_deg)+',0,-180,0,'+str(-res_deg)+',90]',
	region = exportingGeometry,
	maxPixels = int(1e13),
	pyramidingPolicy = {".default": pyramidingPolicy}
);
exportTask.start()
print('Image export task started, moving on')

# Reduce bootstrap images to standard deviation
stdDevImage = ee.ImageCollection.fromImages(list(map(bootstrapFunc, fcList))).reduce(
	reducer = ee.Reducer.stdDev()
)

# export bootstrapping standard deviation
exportTask = ee.batch.Export.image.toAsset(
	image = stdDevImage.toFloat(),
	description = classProperty+'_PredictedBootstrapStdDev',
	assetId = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_PredictedBootstrapStdDev' ,
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
    fcOfInterest = preppedCollection[covariateList]
    propOfVariance = propOfVariance
    
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
PCA_int_ext = assessExtrapolation(preppedCollection[covariateList], propOfVariance).rename('PCA_pct_int_ext')
exportTask = ee.batch.Export.image.toAsset(
	image = PCA_int_ext.toFloat(),
	description = classProperty+'_PCA_pct_int_ext',
	assetId = 'users/'+usernameFolderString+'/'+projectFolder+'/'+classProperty+'_PCA_pct_int_ext',
	crs = CRStoUse,
	crsTransform = '['+str(res_deg)+',0,-180,0,'+str(-res_deg)+',90]',
	region = exportingGeometry,
	maxPixels = int(1e13),
	pyramidingPolicy = {".default": pyramidingPolicy}
);
exportTask.start()
print('Image export task started, moving on')

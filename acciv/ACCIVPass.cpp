#include "ACCIVPass.h"

#include "MakeDirectory.h"

#include "ParameterFileReader.h"

#include "TiePointSet.h"

#include "Correlator.h"

#include "ScatteredVectorData2D.h"
#include "ScatteredScalarData2D.h"
#include "VectorField2D.h"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <math.h>
#include <cstdlib>
#include <algorithm>

bool ACCIVPass::doPass(const UString & inFolder)
{
	folder = inFolder;
printf("beginning pass\n");
	if(!getParameters())
		return false;
printf("parameters read\n");

	if(!setUpWorkDirectory())
		return false;
printf("created work directory\n");


	if(!getGridInfo())
		return false;
printf("grid read\n");

	setUpTiePointAndVelocityTemplates();

	if(!getImageTimes())
		return false;
printf("image times read\n");

	if(!getImagePairs())
		return false;

	if(advectImages)
	{
		computeAdvectionTime();
	}

	if(advectImages && doImageAdvectionStep)
	{
		if(!doImageAdvection())
			return false;
printf("images advected\n");
	}

	if(doImageCorrelationStep)
	{
		if(!correlateImages())
			return false;
printf("images correlated\n");
		if(!combineTiePoints())
			return false;
printf("tie points combined\n");
	}
	else
	{
		combinedTiePointsFileName = combinedCorrelationTiePointsFileName;
	}

	if(advectImages)
	{
		if(doUndoTiePointAdvectionStep)
		{
			if(!undoTiePointAdvection())
				return false;
printf("tie points unadvected\n");
		}
		else
		{
			combinedTiePointsFileName = deadvectedTiePointFileName;
		}
	}


	if(doVelocityConstructionStep)
	{
		if(!constructVelocities())
			return false;
printf("velocities constructed\n");
	}


	return true;
}

bool ACCIVPass::getParameters()
{
	UString fileName = UString::makeFormat("%s/parameters.ascii",(const char *)folder);
	ParameterFileReader reader(fileName);

	UString defaultParameterFile;
	UString tag;
	tag = "defaultParameterFile";
	if(reader.getString(tag, defaultParameterFile))
	{
		defaultParameterFile = UString::makeFormat("%s/%s",(const char *)folder,(const char *)defaultParameterFile);
		reader.addFile(defaultParameterFile);
	}

	// make a list of image file names based on the template
	//   and indices
	UString imageTemplate;
	tag = "earlierImageIndices";
	if(!reader.getIntegerList(tag, earlierImageIndices))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	tag = "earlierImageFileTemplate";
	if(!reader.getString(tag, imageTemplate))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}

	for(SInt32 index=0; index < earlierImageIndices.getSize(); index++)
	{
		UString fileName = UString::makeFormat(imageTemplate, earlierImageIndices[index]);
		fileName = UString::makeFormat("%s/%s", (const char *)folder, (const char *)fileName);
		earlierImageFileNames.add(fileName);
	}

	tag = "laterImageIndices";
	if(!reader.getIntegerList(tag, laterImageIndices))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	tag = "laterImageFileTemplate";
	if(!reader.getString(tag, imageTemplate))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}

	for(SInt32 index=0; index < laterImageIndices.getSize(); index++)
	{
		UString fileName = UString::makeFormat(imageTemplate, laterImageIndices[index]);
		fileName = UString::makeFormat("%s/%s", (const char *)folder, (const char *)fileName);
		laterImageFileNames.add(fileName);
	}

	tag = "gridGeometryFactorsFileName";
	if(!reader.getString(tag, gridGeometryFactorsFileName))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	gridGeometryFactorsFileName = UString::makeFormat("%s/%s", (const char *)folder, (const char *)gridGeometryFactorsFileName);
	tag = "outGridVelocityFileName";
	if(!reader.getString(tag, outGridVelocityFileName))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	outGridVelocityFileName = UString::makeFormat("%s/%s", (const char *)folder, (const char *)outGridVelocityFileName);
	tag = "outScatteredVelocityFileName";
	if(!reader.getString(tag, outScatteredVelocityFileName))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	outScatteredVelocityFileName = UString::makeFormat("%s/%s", (const char *)folder, (const char *)outScatteredVelocityFileName);

	tag = "advectImages";
	if(!reader.getBool(tag, advectImages))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	if(advectImages)
	{
		tag = "advectionTimeFlag";
		if(!reader.getString(tag, advectionTimeFlag))
		{
			fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
			return false;
		}

		tag = "advectionTimeFlag";
		if(!reader.getString(tag, advectionTimeFlag))
		{
			fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
			return false;
		}

		UString fileName;
		tag = "inEarlierVelocityFileName";
		if(!reader.getString(tag, fileName))
		{
			fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
			return false;
		}
		fileName = UString::makeFormat("%s/%s", (const char *)folder, (const char *)fileName);
		inGridVelocityFileNames.add(fileName);
		tag = "inLaterVelocityFileName";
		if(!reader.getString(tag, fileName))
		{
			fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
			return false;
		}
		fileName = UString::makeFormat("%s/%s", (const char *)folder, (const char *)fileName);
		inGridVelocityFileNames.add(fileName);
	}

	tag = "correlationBoxSize";
	if(!reader.getIntegerList(tag, correlationBoxSize) || (correlationBoxSize.getSize() != 2))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Problem with parsing tag %s\n", (const char *)tag);
		return false;
	}
	tag = "searchRange";
	if(!reader.getIntegerList(tag, searchRange) || (searchRange.getSize() != 4))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Problem with parsing tag %s\n", (const char *)tag);
		return false;
	}
	tag = "stride";
	if(!reader.getIntegerList(tag, stride) || (stride.getSize() != 2))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Problem with parsing tag %s\n", (const char *)tag);
		return false;
	}
	tag = "correlationTolerance";
	if(!reader.getDouble(tag, correlationTolerance))
	{
		correlationTolerance = 1e-2;
	}

	tag = "minimumCorrelationCoefficient";
	if(!reader.getDouble(tag, minimumCorrelationCoefficient))
	{
		minimumCorrelationCoefficient = 0.5;
	}

	tag = "minimumTimeSeparation";
	if(!reader.getDouble(tag, minimumTimeSeparation))
	{
		minimumTimeSeparation = 0.0;
	}

	tag = "maximumTimeSeparation";
	if(!reader.getDouble(tag, maximumTimeSeparation))
	{
		maximumTimeSeparation = 1e8;
	}

	tag = "smoothFitOutlierRemovalIterationCount";
	if(!reader.getInteger(tag, smoothFitOutlierRemovalIterationCount))
	{
		smoothFitOutlierRemovalIterationCount = 1;
	}
	tag = "streamlineFollowingSmoothFitIterationCount";
	if(!reader.getInteger(tag, streamlineFollowingSmoothFitIterationCount))
	{
		streamlineFollowingSmoothFitIterationCount = 2;
	}
	tag = "smoothFitMinControlPointScatteredNeighbors";
	if(!reader.getInteger(tag, smoothFitMinControlPointScatteredNeighbors))
	{
		smoothFitMinControlPointScatteredNeighbors = 16;
	}
	tag = "pathVectorCount";
	if(!reader.getInteger(tag, pathVectorCount))
	{
		pathVectorCount = 13;
	}
	tag = "advectionErrorTolerance";
	if(!reader.getDouble(tag, advectionErrorTolerance))
	{
		advectionErrorTolerance = 1e-4;
	}
	tag = "advectionMaxTimeStepCount";
	if(!reader.getInteger(tag, advectionMaxTimeStepCount))
	{
		advectionMaxTimeStepCount = 10000;
	}
	tag = "outlierThresholdConstant";
	if(!reader.getDouble(tag, outlierThresholdConstant))
	{
		outlierThresholdConstant = 6.0;
	}


	tag = "maskFinalVelocity";
	if(!reader.getBool(tag, maskFinalVelocity))
	{
		maskFinalVelocity = false;
	}

	tag = "workDirectory";
	if(!reader.getString(tag, workDirectory))
	{
		workDirectory = "_work";
	}
	workDirectory = UString::makeFormat("%s/%s", (const char *)folder, (const char *)workDirectory);

	tag = "doImageAdvectionStep";
	if(!reader.getBool(tag, doImageAdvectionStep))
	{
		doImageAdvectionStep = true;
	}
	tag = "doImageCorrelationStep";
	if(!reader.getBool(tag, doImageCorrelationStep))
	{
		doImageCorrelationStep = true;
	}
	tag = "doUndoTiePointAdvectionStep";
	if(!reader.getBool(tag, doUndoTiePointAdvectionStep))
	{
		doUndoTiePointAdvectionStep = true;
	}
	tag = "doVelocityConstructionStep";
	if(!reader.getBool(tag, doVelocityConstructionStep))
	{
		doVelocityConstructionStep = true;
	}
	tag = "overwriteWorkFiles";
	if(!reader.getBool(tag, overwriteWorkFiles))
	{
		fprintf(stderr, "ACCIVPass::getParameters: Could not read tag %s\n", (const char *)tag);
		fprintf(stderr, "  Using default value of overwriteWorkFiles = false\n");
		overwriteWorkFiles = false;
	}
	return true;
}

bool ACCIVPass::setUpWorkDirectory()
{
	return makeDirectory(workDirectory);
}

bool ACCIVPass::getGridInfo()
{

	if(!gridVelocityScaling.read(gridGeometryFactorsFileName))
	{
		fprintf(stderr, "ACCIVPass::getGridInfo: Could not read grid velocity scaling from file %s\n", (const char *)gridGeometryFactorsFileName);
		return false;
	}

	imageSize.setSize(2);
	imageSize[0] = gridVelocityScaling.getXSize();
	imageSize[1] = gridVelocityScaling.getYSize();

	imageBounds.setSize(4);
	imageBounds[0] = gridVelocityScaling.getXLower();
	imageBounds[1] = gridVelocityScaling.getXUpper();
	imageBounds[2] = gridVelocityScaling.getYLower();
	imageBounds[3] = gridVelocityScaling.getYUpper();
	
	return true;
}

bool ACCIVPass::getImageTimes()
{
	for(SInt32 index=0; index < earlierImageFileNames.getSize(); index++)
	{
		double time;
		if(!MaskedImage::readTimeHDF5File(earlierImageFileNames[index], time))
			return false;
		earlierImageTimes.add(time);
	}
	for(SInt32 index=0; index < laterImageFileNames.getSize(); index++)
	{
		double time;
		if(!MaskedImage::readTimeHDF5File(laterImageFileNames[index], time))
			return false;
		laterImageTimes.add(time);
	}
	return true;
}

bool ACCIVPass::getImagePairs()
{
	// find the pairs of images that need to be correlated
	//   the first image should always have time < that of the second
	//   image and the images cannot be the same
	maxSeparationTime = 0;

	for(SInt32 earlierIndex = 0; earlierIndex < earlierImageIndices.getSize(); earlierIndex++)
	{
		for(SInt32 laterIndex = 0; laterIndex < laterImageIndices.getSize(); laterIndex++)
		{
			double deltaT = laterImageTimes[laterIndex] - earlierImageTimes[earlierIndex];
			if((deltaT > minimumTimeSeparation) && (deltaT < maximumTimeSeparation))
			{
				imagePairs.add(earlierIndex);
				imagePairs.add(laterIndex);
				maxSeparationTime = std::max(maxSeparationTime,deltaT);
			}
		}
	}
	printf("maximum time separating valid image pairs: %g\n", maxSeparationTime);

	return true;
}

void ACCIVPass::computeAdvectionTime()
{
	// compute the time to advect the images to
	double meanEarlierTime = 0.0;
	double meanLaterTime = 0.0;
	for(SInt32 earlierIndex = 0; earlierIndex < earlierImageTimes.getSize(); earlierIndex++)
	{
		meanEarlierTime += earlierImageTimes[earlierIndex];
	}
	meanEarlierTime /= earlierImageTimes.getSize();

	for(SInt32 laterIndex = 0; laterIndex < laterImageTimes.getSize(); laterIndex++)
	{
		meanLaterTime += laterImageTimes[laterIndex];
	}
	meanLaterTime /= laterImageTimes.getSize();

	advectionTime = 0.5*(meanEarlierTime + meanLaterTime);
	if(advectionTimeFlag == "early")
	{
		advectionTime = meanEarlierTime;
	}
	else if(advectionTimeFlag == "middle")
	{
		advectionTime = 0.5*(meanEarlierTime + meanLaterTime);
	}
	else if(advectionTimeFlag == "late")
	{
		advectionTime = meanLaterTime;
	}
	else
	{
		fprintf(stderr, "ACCIVPass::doImageAdvection: bad value for advectionTimeFlag\n");
		U_ASSERT(false);
	}
}

bool ACCIVPass::doImageAdvection()
{
	UArray<UString> earlierAdvectedImageFileNames, laterAdvectedImageFileNames;
	for(SInt32 earlierIndex = 0; earlierIndex < earlierImageIndices.getSize(); earlierIndex++)
	{
		earlierAdvectedImageFileNames.add(UString::makeFormat(earlierAdvectedImageTemplate, earlierImageIndices[earlierIndex]));
	}
	for(SInt32 laterIndex = 0; laterIndex < laterImageIndices.getSize(); laterIndex++)
	{
		laterAdvectedImageFileNames.add(UString::makeFormat(laterAdvectedImageTemplate, laterImageIndices[laterIndex]));
	}

	// load the advection velocity
	VectorField2D velocity;
	if(!velocity.read(inGridVelocityFileNames[0]))
		return false;

printf("advecting images\n");
	// for each early image
	for(SInt32 earlierIndex = 0; earlierIndex < earlierImageIndices.getSize(); earlierIndex++)
	{
printf(" %s\n", (const char *)earlierAdvectedImageFileNames[earlierIndex]);
		if(!overwriteWorkFiles)
		{
			if(hdf5FileExists(earlierAdvectedImageFileNames[earlierIndex]))
				continue; // the advected image alread exists
		}
		MaskedImage image;
		// load the image
		if(!image.read(earlierImageFileNames[earlierIndex]))
			return false;


		double deltaT = advectionTime - earlierImageTimes[earlierIndex];

		// call the image advection function
		if(!image.advect(velocity, deltaT, advectionErrorTolerance, advectionMaxTimeStepCount))
			return false;

		// save the image
		if(!image.write(earlierAdvectedImageFileNames[earlierIndex]))
			return false;
	}

	// load the advection velocity
	if(!velocity.read(inGridVelocityFileNames[1]))
		return false;

	// for each late image
	for(SInt32 laterIndex = 0; laterIndex < laterImageIndices.getSize(); laterIndex++)
	{
printf(" %s\n", (const char *)laterAdvectedImageFileNames[laterIndex]);
		if(!overwriteWorkFiles)
		{
			if(hdf5FileExists(laterAdvectedImageFileNames[laterIndex]))
				continue; // the advected image alread exists
		}
		// make sure this isn't redundant to advection that has already occurred
		bool alreadyAdvected = false;
		SInt32 earlierIndex = 0;
		while(earlierIndex < earlierImageIndices.getSize())
		{
			if(laterImageFileNames[laterIndex] == earlierImageFileNames[earlierIndex])
			{
				alreadyAdvected = true;
				break;
			}
			earlierIndex++;
		}
		if(alreadyAdvected)
		{
			MaskedImage image;
			if(!image.read(earlierAdvectedImageFileNames[earlierIndex]))
				return false;
			if(!image.write(laterAdvectedImageFileNames[laterIndex]))
				return false;
			continue;
		}
		// load the images
		MaskedImage image;
		if(!image.read(laterImageFileNames[laterIndex]))
			return false;


		double deltaT = advectionTime - laterImageTimes[laterIndex];

		// call the image advection function
		image.advect(velocity, deltaT, advectionErrorTolerance, advectionMaxTimeStepCount);

		// save the image
		if(!image.write(laterAdvectedImageFileNames[laterIndex]))
			return false;
	}

	earlierImageFileNames = earlierAdvectedImageFileNames;
	laterImageFileNames = laterAdvectedImageFileNames;
	return true;
}

void ACCIVPass::setUpTiePointAndVelocityTemplates()
{
	earlierAdvectedImageTemplate = UString::makeFormat("%s/earlierAdvectedImage_%%i.h5", (const char *)workDirectory);
	laterAdvectedImageTemplate = UString::makeFormat("%s/laterAdvectedImage_%%i.h5", (const char *)workDirectory);

	correlationTiePointTemplate = UString::makeFormat("%s/correlationTiePoints_%%i_%%i.h5", (const char *)workDirectory);
	combinedCorrelationTiePointsFileName = UString::makeFormat("%s/combinedCorrelationTiePoints.h5", (const char *)workDirectory);
	deadvectedTiePointFileName = UString::makeFormat("%s/deadvectedTiePoints.h5", (const char *)workDirectory);

	noCurvedPathScatteredVelocityFileName = UString::makeFormat("%s/noCurvedPathScatteredVelocity.h5", (const char *)workDirectory);
	noCurvedPathGridVelocityFileName = UString::makeFormat("%s/noCurvedPathGridVelocity.h5", (const char *)workDirectory);

	curvedPathScatteredVelocityTemplate = UString::makeFormat("%s/curvedPathScatteredVelocity_%%i.h5", (const char *)workDirectory);
	curvedPathGridVelocityTemplate = UString::makeFormat("%s/curvedPathGridVelocity_%%i.h5", (const char *)workDirectory);
	noOutlierScatteredVelocityTemplate = UString::makeFormat("%s/noOutlierScatteredVelocity_%%i.h5", (const char *)workDirectory);
	noOutlierGridVelocityTemplate = UString::makeFormat("%s/noOutlierGridVelocity_%%i.h5", (const char *)workDirectory);
	noOutlierTiePointsTemplate = UString::makeFormat("%s/noOutlierTiePoints_%%i.h5", (const char *)workDirectory);
}

bool ACCIVPass::correlateImages()
{
	Correlator correlator(correlationTolerance);

printf("correlating images\n");
	// for each image pair:
	for(SInt32 pairIndex = 0; pairIndex < imagePairs.getSize(); pairIndex += 2)
	{
printf(" image pair %i %i\n", earlierImageIndices[imagePairs[pairIndex]], laterImageIndices[imagePairs[pairIndex+1]]);
		UString fileName = UString::makeFormat(correlationTiePointTemplate,
			earlierImageIndices[imagePairs[pairIndex]],
			laterImageIndices[imagePairs[pairIndex+1]]);

		if(!overwriteWorkFiles)
		{
			if(hdf5FileExists(fileName))
				continue;
		}
		// load the images
		MaskedImage image1, image2;
		if(!image1.read(earlierImageFileNames[imagePairs[pairIndex]]))
			return false;
		if(!image2.read(laterImageFileNames[imagePairs[pairIndex+1]]))
			return false;

		double deltaT = laterImageTimes[imagePairs[pairIndex+1]]
		    - earlierImageTimes[imagePairs[pairIndex]];
		UArray<SInt32> localSearchRange(4);
		for(SInt32 index = 0; index < 4; index++)
		{
			SInt32 sign = (searchRange[index] > 0 ? 1 : -1);
			SInt32 magnitude = std::max(1,(SInt32)(fabs(searchRange[index])*deltaT/maxSeparationTime + 0.5));
			localSearchRange[index] = sign*magnitude;
		}
		printf("  separation time: %f\n", deltaT);
		printf("  range scaling: %f\n", deltaT/maxSeparationTime);
		printf("  search range [%i %i %i %i]\n", localSearchRange[0], localSearchRange[1], localSearchRange[2],
				localSearchRange[3]);

		// call the image correlation function
		TiePointSet tiePoints;
		if(!correlator.correlateImages(image1, image2, tiePoints,
			correlationBoxSize[0], correlationBoxSize[1], localSearchRange[0], localSearchRange[1], localSearchRange[2],
			localSearchRange[3], stride[0], stride[1]))
		{
			// write a dummy file
			tiePoints.getCorrelationCoefficients().add(-1);
			tiePoints.getX1().add(-1);
			tiePoints.getX2().add(-1);
			tiePoints.getY1().add(-1);
			tiePoints.getY2().add(-1);
			tiePoints.getDeltaTs().add(-1);
			tiePoints.getLowerIndexDeltaT().add(-1);
			tiePoints.getUpperIndexDeltaT().add(-1);
			if(!tiePoints.write(fileName))
				return false;
			continue;
		}

		tiePoints.getDeltaTs().add(deltaT);
		tiePoints.getLowerIndexDeltaT().add(0);
		tiePoints.getUpperIndexDeltaT().add(tiePoints.getX1().getSize()-1);

		// store the tie points and correlaiton coefficients in a file in the work directory

		if(!tiePoints.write(fileName))
			return false;
	}

	return true;
}

bool ACCIVPass::undoTiePointAdvection()
{
	if(!overwriteWorkFiles)
	{
		if(hdf5FileExists(deadvectedTiePointFileName))
		{
			combinedTiePointsFileName = deadvectedTiePointFileName;
			return true;
		}
	}

	// load the advection velocity
	VectorField2D velocity1, velocity2;
	if(!velocity1.read(inGridVelocityFileNames[0]))
		return false;
	if(!velocity2.read(inGridVelocityFileNames[1]))
		return false;

	UArray<double> deltaT1s, deltaT2s;

	// for each image pair:
	for(SInt32 pairIndex = 0; pairIndex < imagePairs.getSize(); pairIndex += 2)
	{
		deltaT1s.add(earlierImageTimes[imagePairs[pairIndex]] - advectionTime);
		deltaT2s.add(laterImageTimes[imagePairs[pairIndex+1]] - advectionTime);
	}
	TiePointSet tiePoints;
	if(!tiePoints.read(combinedTiePointsFileName))
		return false;

	if(!tiePoints.unadvect(velocity1, velocity2, deltaT1s, deltaT2s, advectionErrorTolerance, advectionMaxTimeStepCount))
		return false;

	if(!tiePoints.write(deadvectedTiePointFileName))
		return false;

	combinedTiePointsFileName = deadvectedTiePointFileName;

	return true;
}

bool ACCIVPass::constructVelocities()
{
	printf("constructing velocities\n");

	// call the function to convert tie points to preliminary
	// velocities (without curved paths)
	if(!constructNoCurvatureVelocity())
		return false;

	UString currentGridVelocityFileName = noCurvedPathGridVelocityFileName;
	UString currentScatteredVelocityFileName = noCurvedPathScatteredVelocityFileName;
	UString currentTiePointsFileName = combinedTiePointsFileName;

	printf("constructing grid from scattered velocity (no curvature)\n");
	UInt32 minControlPointScatteredNeighbors = smoothFitMinControlPointScatteredNeighbors; //max(2*imagePairs.getSize(),hardMinPointsPerBin); // max of 4 points per image pair or 16 points
	double meanRes = 0.0;
	UArray<vector2d> residuals;
	if(!constructGridVelocityFromScatteredVelocity(currentScatteredVelocityFileName,
		currentGridVelocityFileName, minControlPointScatteredNeighbors, 1, residuals, meanRes))
		return false;
	printf("  smooth fit mean residual: %g\n", meanRes);

	UInt32 index = 0;
	UInt32 outerCount = (smoothFitOutlierRemovalIterationCount == 0) ? 1 : 3;
	for(UInt32 outerIndex = 0; outerIndex < outerCount; outerIndex++)
	{
		UInt32 count = streamlineFollowingSmoothFitIterationCount;
		if(outerIndex == 1)
		{
			count = smoothFitOutlierRemovalIterationCount;
		}

		double previousRes = 0;
		const double errorThreshold = 0.01;
		for(UInt32 innerIndex = 0; innerIndex < count; innerIndex++)
		{
			printf(" curved-path/smooth-fit iteration %i\n", index);
			if(outerIndex == 1)
			{
				printf(" removing outliers\n");
				UString newTiePointsFileName = UString::makeFormat(noOutlierTiePointsTemplate, index);
				// call the outlier removal function
				if(!removeOutliers(residuals,
					currentTiePointsFileName,
					newTiePointsFileName))
					return false;
				currentTiePointsFileName = newTiePointsFileName;
			}
			else
			{
				currentTiePointsFileName = combinedTiePointsFileName;
			}
			currentScatteredVelocityFileName = UString::makeFormat(curvedPathScatteredVelocityTemplate, index);
			// call the function that computes the new scattered velocity along curved paths
			if(!constructVelocityAlongCurvedPaths(currentTiePointsFileName, currentGridVelocityFileName,
				currentScatteredVelocityFileName))
				return false;

			currentGridVelocityFileName = UString::makeFormat(curvedPathGridVelocityTemplate, index);
				// set the input scattered velocity file name to be the result
			// call the function that computes the smooth velocity
			if(!constructGridVelocityFromScatteredVelocity(
				currentScatteredVelocityFileName,
				currentGridVelocityFileName,
				minControlPointScatteredNeighbors*pathVectorCount,
				pathVectorCount,
				residuals,
				meanRes))
				return false;

			index++;
			double error = abs(previousRes-meanRes)/std::max(previousRes,meanRes);
			printf("  smooth fit mean residual: %g\n", meanRes);
			printf("  residual change: %.4f\n", error);
			if(error < errorThreshold)
				break;
			previousRes = meanRes;
		}
	}

	// store the final smooth velocity and the final scattered velocity without outlier removal

	bool writeGridVelocity = true;
	bool writeScatteredVelocity = true;
	if(!overwriteWorkFiles)
	{
		if(hdf5FileExists(outGridVelocityFileName))
		{
			writeGridVelocity = false;
			if(hdf5FileExists(outScatteredVelocityFileName))
			{
				writeScatteredVelocity = false;
			}
		}
	}
	if(writeGridVelocity)
	{
		// write the final velocity with zeros along the edges
		writeFinalGridVelocity(currentGridVelocityFileName,currentScatteredVelocityFileName,
				outGridVelocityFileName);
	}

	if(writeScatteredVelocity)
	{
		// use the velocity *without* the zeroed-out edges to compute the final scattered velocities
		if(!constructVelocityAlongCurvedPaths(combinedTiePointsFileName, currentGridVelocityFileName,
			outScatteredVelocityFileName))
			return false;
		if(!appendOneSigmaUncertainty(residuals,
				meanRes,
				outScatteredVelocityFileName))
			return false;
	}

	return true;
}

bool ACCIVPass::combineTiePoints()
{
	if(!overwriteWorkFiles)
	{
		if(hdf5FileExists(combinedCorrelationTiePointsFileName))
		{
			combinedTiePointsFileName = combinedCorrelationTiePointsFileName;
			return true;
		}
	}

	TiePointSet allTiePoints;
	UArray<double> & allX1 = allTiePoints.getX1();
	UArray<double> & allY1 = allTiePoints.getY1();
	UArray<double> & allX2 = allTiePoints.getX2();
	UArray<double> & allY2 = allTiePoints.getY2();
	UArray<double> & allCoeffs = allTiePoints.getCorrelationCoefficients();
	UArray<double> & allDeltaTs = allTiePoints.getDeltaTs();
	UArray<SInt32> & allLowers = allTiePoints.getLowerIndexDeltaT();
	UArray<SInt32> & allUppers = allTiePoints.getUpperIndexDeltaT();


	for(SInt32 pairIndex = 0; pairIndex < imagePairs.getSize(); pairIndex += 2)
	{
		UString fileName = UString::makeFormat(correlationTiePointTemplate, 
			earlierImageIndices[imagePairs[pairIndex]],
			laterImageIndices[imagePairs[pairIndex+1]]);
		TiePointSet tiePoints;
		if(!tiePoints.read(fileName))
			continue;

		UArray<double> & x1 = tiePoints.getX1();
		UArray<double> & y1 = tiePoints.getY1();
		UArray<double> & x2 = tiePoints.getX2();
		UArray<double> & y2 = tiePoints.getY2();
		UArray<double> & coeffs = tiePoints.getCorrelationCoefficients();
		UArray<double> & deltaTs = tiePoints.getDeltaTs();
		UArray<SInt32> & lowers = tiePoints.getLowerIndexDeltaT();
		UArray<SInt32> & uppers = tiePoints.getUpperIndexDeltaT();

		if((coeffs.getSize() == 1) && (coeffs[0] == -1))
			// this is just a dummy empty file
			continue;

		for(SInt32 timeIndex = 0; timeIndex < deltaTs.getSize(); timeIndex++)
		{
			SInt32 lowerIndex = allX1.getSize();
			for(SInt32 tiePointIndex = lowers[timeIndex]; tiePointIndex <= uppers[timeIndex]; tiePointIndex++)
			{
				if(coeffs[tiePointIndex] < minimumCorrelationCoefficient)
					continue;
				allX1.add(x1[tiePointIndex]);
				allY1.add(y1[tiePointIndex]);
				allX2.add(x2[tiePointIndex]);
				allY2.add(y2[tiePointIndex]);
				allCoeffs.add(coeffs[tiePointIndex]);
			}
			SInt32 upperIndex = allX1.getSize()-1;
			if(lowerIndex <= upperIndex)
			{
				allDeltaTs.add(deltaTs[timeIndex]);
				allLowers.add(lowerIndex);
				allUppers.add(upperIndex);
			}
		}
	}

	if((allX1.getSize() == 0) || !allTiePoints.write(combinedCorrelationTiePointsFileName))
		return false;

	combinedTiePointsFileName = combinedCorrelationTiePointsFileName;
	if(advectImages && !appendCorrelationUncertainties())
		return false;

	return true;
}

bool ACCIVPass::constructNoCurvatureVelocity()
{
	if(!overwriteWorkFiles)
	{
		if(hdf5FileExists(noCurvedPathScatteredVelocityFileName))
		{
			return true;
		}
	}

	ScatteredVectorData2D scatteredVelocity;

	UArray<double> & X = scatteredVelocity.getX();
	UArray<double> & Y = scatteredVelocity.getY();
	UArray<vector2d> & V = scatteredVelocity.getData();
	UArray<double> & weights = scatteredVelocity.getWeights();

	UString tiePointsFileName = combinedTiePointsFileName;

	TiePointSet tiePoints;
	if(!tiePoints.read(tiePointsFileName))
		return false;
	UArray<double> & x1 = tiePoints.getX1();
	UArray<double> & y1 = tiePoints.getY1();
	UArray<double> & x2 = tiePoints.getX2();
	UArray<double> & y2 = tiePoints.getY2();
	UArray<double> & coeffs = tiePoints.getCorrelationCoefficients();
	UArray<double> & deltaTs = tiePoints.getDeltaTs();
	UArray<SInt32> & lowers = tiePoints.getLowerIndexDeltaT();
	UArray<SInt32> & uppers = tiePoints.getUpperIndexDeltaT();
	for(SInt32 timeIndex = 0; timeIndex < deltaTs.getSize(); timeIndex++)
	{
		double deltaT = deltaTs[timeIndex];
		for(SInt32 tiePointIndex = lowers[timeIndex]; tiePointIndex <= uppers[timeIndex]; tiePointIndex++)
		{
			double x = 0.5*(x1[tiePointIndex] + x2[tiePointIndex]);
			double y = 0.5*(y1[tiePointIndex] + y2[tiePointIndex]);
			V.add(vector2d((x2[tiePointIndex] - x1[tiePointIndex])/deltaT,
					(y2[tiePointIndex] - y1[tiePointIndex])/deltaT));
			X.add(x);
			Y.add(y);
			weights.add(coeffs[tiePointIndex]*deltaT/maxSeparationTime);
		}
		scatteredVelocity.setSize(V.getSize());
	}

	if(!scatteredVelocity.write(noCurvedPathScatteredVelocityFileName))
		return false;

	if(!appendScaledScatteredVelocity(noCurvedPathScatteredVelocityFileName))
		return false;

	return true;
}


bool ACCIVPass::constructGridVelocityFromScatteredVelocity(const UString & inOutScatteredVelocityFileName,
															const UString & outGridVelocityFileName,
															UInt32 minControlPointScatteredNeighbors,
															UInt32 pathVectorCount,
															UArray<vector2d> & outResiduals,
															double & meanRes)
{
	if(!overwriteWorkFiles)
	{
		if(hdf5FileExists(outGridVelocityFileName))
		{
			hid_t fileID = H5Fopen(outGridVelocityFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
			if(fileID < 0)
			{
				fprintf(stderr, "ACCIVPass::constructGridVelocityFromScatteredVelocity: Could not open file %s\n", (const char *)outGridVelocityFileName);
				return false;
			}

			UString dataSetName = "/meanResidual";
			hsize_t dimensions[1];
			H5T_class_t typeClass;
			size_t typeSize;
			herr_t status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
			if(status < 0)
			{
				dataSetName.makeUpper();
				status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
			}
			if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || (dimensions[0] != 1))
			{
				fprintf(stderr, "ACCIVPass::constructGridVelocityFromScatteredVelocity: Could not get dataset info for %s in file %s\n",
					(const char *)dataSetName, (const char *)outGridVelocityFileName);
				return false;
			}

			status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, &meanRes);
			if(status < 0)
			{
				fprintf(stderr, "ACCIVPass::constructGridVelocityFromScatteredVelocity: Could not read dataset %s in file %s\n",
					(const char *)dataSetName, (const char *)outGridVelocityFileName);
				return false;
			}
			return true;
		}
	}

	ScatteredVectorData2D scatteredVelocity;
	if(!scatteredVelocity.read(inOutScatteredVelocityFileName))
		return false;

	VectorField2D gridVelocity(imageSize[0], imageSize[1], imageBounds[0], imageBounds[1], imageBounds[2], imageBounds[3]);
	
	meanRes = 0.0;

	scatteredVelocity.smoothFit(gridVelocity, minControlPointScatteredNeighbors, outResiduals, meanRes);

	if(!gridVelocity.write(outGridVelocityFileName))
		return false;

	// add the mean residual
	hid_t fileID = H5Fopen(outGridVelocityFileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ACCIVPass::constructGridVelocityFromScatteredVelocity: Could not open file %s\n",
			(const char *)outGridVelocityFileName);
		return false;
	}

	UString dataSetName = "/meanResidual";
	hsize_t dimensions[2];
	dimensions[0] = 1;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, &meanRes);
	if((status < 0))
	{
		fprintf(stderr, "ACCIVPass::constructGridVelocityFromScatteredVelocity: Could not make dataset %s in file %s\n",
			(const char *)dataSetName, (const char *)outGridVelocityFileName);
		return false;
	}

	status = H5Fclose(fileID);

	if(!appendOneSigmaUncertainty(outResiduals,
			meanRes,
			inOutScatteredVelocityFileName))
		return false;

	if(!appendScaledGridVelocity(outGridVelocityFileName))
		return false;

	return true;
}

bool ACCIVPass::writeFinalGridVelocity(const UString & inGridVelocityFileName,
		const UString & inScatteredVelocityFileName,
		const UString & inFinalGridVelocityFileName)
{
	VectorField2D gridVelocity;
	ScatteredVectorData2D scatteredVelocity;
	if(!gridVelocity.read(inGridVelocityFileName))
		return false;

	if(!scatteredVelocity.read(inScatteredVelocityFileName))
		return false;

	if(maskFinalVelocity)
	{
		UInt32 nx = gridVelocity.getXSize();
		UInt32 ny = gridVelocity.getYSize();
		double xl = gridVelocity.getXLower();
		double xu = gridVelocity.getXUpper();
		double yl = gridVelocity.getYLower();
		double yu = gridVelocity.getYUpper();
		double deltaX = (xu - xl)/double(nx-1);
		double deltaY = (yu - yl)/double(ny-1);



		ScatteredScalarData2D scatteredMask(scatteredVelocity.getSize()+2*nx+2*ny-4);
		for(UInt32 index = 0; index < scatteredVelocity.getSize(); index++)
		{
			scatteredMask.getX(index) = scatteredVelocity.getX(index);
			scatteredMask.getY(index) = scatteredVelocity.getY(index);
			scatteredMask.getData(index) = 1.0;
		}
		for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
		{
			UInt32 index = scatteredVelocity.getSize() + xIndex;
			scatteredMask.getX(index) = gridVelocity.getX()[xIndex];
			scatteredMask.getY(index) = yl;
			scatteredMask.getData(index) = 0.0;
			index = scatteredVelocity.getSize() + xIndex + nx;
			scatteredMask.getX(index) = gridVelocity.getX()[xIndex];
			scatteredMask.getY(index) = yu;
			scatteredMask.getData(index) = 0.0;
		}
		for(UInt32 yIndex = 1; yIndex < ny-1; yIndex++)
		{
			UInt32 index = scatteredVelocity.getSize() + 2*nx + yIndex-1;
			scatteredMask.getX(index) = xl;
			scatteredMask.getY(index) = gridVelocity.getY()[yIndex];
			scatteredMask.getData(index) = 0.0;
			index = scatteredVelocity.getSize() + 2*nx + (ny-2) + yIndex-1;
			scatteredMask.getX(index) = xu;
			scatteredMask.getY(index) = gridVelocity.getY()[yIndex];
			scatteredMask.getData(index) = 0.0;
		}

		ScalarField2D gridMask(nx,ny,xl,xu,yl,yu);
		UArray<double> residuals;
		double meanRes;
		scatteredMask.smoothFit(gridMask,1,residuals,meanRes);

		for(UInt32 index = 0; index < gridVelocity.getDataSize(); index++)
		{
			// make sure to avoid overshoot
			double mask = std::min(1.0,std::max(0.0,gridMask[index]));
			gridVelocity[index] *= mask;
		}
	}



/*	// compute the maximum and minimum scattered x and y
	double maxX = scatteredVelocity.getX(0);
	double minX = scatteredVelocity.getX(0);
	double maxY = scatteredVelocity.getY(0);
	double minY = scatteredVelocity.getY(0);
	for(UInt32 index = 1; index < scatteredVelocity.getSize(); index++)
	{
		maxX = std::max(maxX,scatteredVelocity.getX(index));
		minX = std::min(minX,scatteredVelocity.getX(index));
		maxY = std::max(maxY,scatteredVelocity.getY(index));
		minY = std::min(minY,scatteredVelocity.getY(index));
	}

	// use a cubic function to zero out the edges of the velocity field where there is no data
	// f(x) = a*x^3 + b*x^2 + c*x + d
	// f(0) = 0 = d
	// f'(0) = 0 = c
	// f(1) = 1 = a+b+c+d
	// f'(1) = 0 = 3*a+2*b+c
	// a = -2; b = 3; c = 0; d = 0
	// f(x) = x^2*(3-2*x)

	// make the region for zeroing out the edges at least 12 pixels wide
	double xlPixels = std::max((minX-xl)/deltaX,12.0);
	double xuPixels = std::max((maxX-xu)/deltaX,12.0);
	double ylPixels = std::max((minY-yl)/deltaY,12.0);
	double yuPixels = std::max((maxY-yu)/deltaY,12.0);

	for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
	{
		double x = xIndex/xlPixels;
		if(x > 1.0)
			break;
		double scale = x*x*(3-2*x);

		for(UInt32 yIndex = 0; yIndex < ny; yIndex++)
		{
			gridVelocity(xIndex,yIndex) *= scale;
		}
	}

	for(UInt32 xIndex = nx-1; xIndex >= 0; xIndex--)
	{
		double x = (nx-1-xIndex)/xuPixels;
		if(x > 1.0)
			break;
		double scale = x*x*(3-2*x);

		for(UInt32 yIndex = 0; yIndex < ny; yIndex++)
			gridVelocity(xIndex,yIndex) *= scale;
	}

	for(UInt32 yIndex = 0; yIndex < ny; yIndex++)
	{
		double y = yIndex/ylPixels;
		if(y > 1.0)
			break;
		double scale = y*y*(3-2*y);

		for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
			gridVelocity(xIndex,yIndex) *= scale;
	}

	for(UInt32 yIndex = ny-1; yIndex >= 0; yIndex--)
	{
		double y = (ny-1-yIndex)/yuPixels;
		if(y > 1.0)
			break;
		double scale = y*y*(3-2*y);

		for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
			gridVelocity(xIndex,yIndex) *= scale;
	}*/

	if(!gridVelocity.write(outGridVelocityFileName))
		return false;
	if(!appendScaledGridVelocity(outGridVelocityFileName))
		return false;

	return true;

}


bool ACCIVPass::appendOneSigmaUncertainty(const UArray<vector2d> & residuals,
									double meanResidual,
									const UString & inOutScatteredVelocityFileName)
{
	ScatteredVectorData2D scatteredVelocity;
	if(!scatteredVelocity.read(inOutScatteredVelocityFileName))
		return false;

	UArray<vector2d> & vScatter = scatteredVelocity.getData();

	SInt32 velocityPointCount = vScatter.getSize();

	UArray<vector2d> vScale(velocityPointCount);
	// interpolate to get a velocity scaling at each scattered point
	gridVelocityScaling.interpolate(scatteredVelocity.getX(),scatteredVelocity.getY(), vScale);

	printf("  smooth fit RMS one-sigma uncertainty: %g\n", meanResidual);
	
	hid_t fileID = H5Fopen(inOutScatteredVelocityFileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ACCIVPass::appendOneSigmaUncertainty: Could not open file %s\n",
                  (const char *)inOutScatteredVelocityFileName);
		return false;
	}

	UString dataSetName = "/rmsOneSigmaUncertainty";
	hsize_t dimensions[1];
	dimensions[0] = 1;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, &meanResidual);
	if((status < 0))
	{
		fprintf(stderr, "ACCIVPass::appendOneSigmaUncertainty: Could not make dataset %s in file %s\n",
                  (const char *)dataSetName, (const char *)inOutScatteredVelocityFileName);
		return false;
	}

	UArray<double> data(velocityPointCount);
	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		data[index] = residuals[index].x*vScale[index].x;
	}
	dataSetName = "/residualX";
	dimensions[0] = (hsize_t)velocityPointCount;
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE,
			data.getData());
	if((status < 0))
	{
		fprintf(stderr, "ACCIVPass::appendOneSigmaUncertainty: Could not make dataset %s in file %s\n",
				  (const char *)dataSetName, (const char *)inOutScatteredVelocityFileName);
		return false;
	}

	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		data[index] = residuals[index].y*vScale[index].y;
	}
	dataSetName = "/residualY";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE,
			data.getData());
	if((status < 0))
	{
		fprintf(stderr, "CCIVPass::appendOneSigmaUncertainty: Could not make dataset %s in file %s\n",
			(const char *)dataSetName, (const char *)inOutScatteredVelocityFileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

bool ACCIVPass::appendScaledScatteredVelocity(const UString & inOutScatteredVelocityFileName)
{
	ScatteredVectorData2D scatteredVelocity;
	if(!scatteredVelocity.read(inOutScatteredVelocityFileName))
		return false;

	UArray<vector2d> & vScatter = scatteredVelocity.getData();

	SInt32 velocityPointCount = vScatter.getSize();

	UArray<vector2d> vScale(velocityPointCount);
	// interpolate to get a velocity scaling at each scattered point
	gridVelocityScaling.interpolate(scatteredVelocity.getX(),scatteredVelocity.getY(), vScale);

	hid_t fileID = H5Fopen(inOutScatteredVelocityFileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ACCIVPass::appendScaledVelocities: Could not open file %s\n",
				  (const char *)inOutScatteredVelocityFileName);
		return false;
	}

	UArray<double> data(velocityPointCount);
	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		data[index] = vScatter[index].x*vScale[index].x;
	}
	UString dataSetName = "/vx";
	hsize_t dimensions[1];
	dimensions[0] = (hsize_t)velocityPointCount;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE,
			data.getData());
	if((status < 0))
	{
		fprintf(stderr, "ACCIVPass::appendScaledVelocities: Could not make dataset %s in file %s\n",
				  (const char *)dataSetName, (const char *)inOutScatteredVelocityFileName);
		return false;
	}

	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		data[index] = vScatter[index].y*vScale[index].y;
	}
	dataSetName = "/vy";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE,
			data.getData());
	if((status < 0))
	{
		fprintf(stderr, "CCIVPass::appendScaledVelocities: Could not make dataset %s in file %s\n",
			(const char *)dataSetName, (const char *)inOutScatteredVelocityFileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

bool ACCIVPass::appendScaledGridVelocity(const UString & inOutGridVelocityFileName)
{

	VectorField2D gridVelocity;
	if(!gridVelocity.read(inOutGridVelocityFileName))
		return false;
	UArray<vector2d> & vGrid = gridVelocity.getData();

	UArray<vector2d> & vScale = gridVelocityScaling.getData();

	hid_t fileID = H5Fopen(inOutGridVelocityFileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ACCIVPass::appendScaledVelocities: Could not open file %s\n",
				  (const char *)inOutGridVelocityFileName);
		return false;
	}

    UArray<double> data(vGrid.getSize());
	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		 data[index] = vGrid[index].x*vScale[index].x;
	}
	UString dataSetName = "/vx";
	hsize_t dimensions[2];
	dimensions[0] = (hsize_t)gridVelocity.getYSize();
	dimensions[1] = (hsize_t)gridVelocity.getXSize();
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 2, dimensions, H5T_NATIVE_DOUBLE,
			data.getData());
	if((status < 0))
	{
		fprintf(stderr, "ACCIVPass::appendScaledVelocities: Could not make dataset %s in file %s\n",
				  (const char *)dataSetName, (const char *)inOutGridVelocityFileName);
		return false;
	}

	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		 data[index] = vGrid[index].y*vScale[index].y;
	}
	dataSetName = "/vy";
	status = H5LTmake_dataset(fileID, dataSetName, 2, dimensions, H5T_NATIVE_DOUBLE,
			data.getData());
	if((status < 0))
	{
		fprintf(stderr, "CCIVPass::appendScaledVelocities: Could not make dataset %s in file %s\n",
			(const char *)dataSetName, (const char *)inOutGridVelocityFileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

bool ACCIVPass::appendCorrelationUncertainties()
{
	TiePointSet allTiePoints;
	if(!allTiePoints.read(combinedTiePointsFileName))
	{
		fprintf(stderr, "ACCIVPass::appendCorrelationUncertainties: Could not read tie points from file %s\n",
                  (const char *)combinedTiePointsFileName);
		return false;
	}

	UArray<double> & x1 = allTiePoints.getX1();
	UArray<double> & y1 = allTiePoints.getY1();
	UArray<double> & x2 = allTiePoints.getX2();
	UArray<double> & y2 = allTiePoints.getY2();
	UArray<double> & deltaTs = allTiePoints.getDeltaTs();
	UArray<SInt32> & lowers = allTiePoints.getLowerIndexDeltaT();
	UArray<SInt32> & uppers = allTiePoints.getUpperIndexDeltaT();

	SInt32 velocityPointCount = x1.getSize();

	UArray<double> x(velocityPointCount), y(velocityPointCount);
	for(SInt32 index = 0; index < velocityPointCount; index++)
	{
		x[index] = 0.5*(x1[index]+x2[index]);
		y[index] = 0.5*(y1[index]+y2[index]);
	}

	UArray<vector2d> vScale(velocityPointCount);
	// interpolate to get a velocity scaling at each scattered point
	gridVelocityScaling.interpolate(x,y, vScale);

	UArray<double> correlationLocationRes(velocityPointCount), correlationVelocityRes(velocityPointCount);
	double rmsLocationRes = 0, rmsVelocityRes = 0;
	SInt32 tIndex = 0;
	for(SInt32 index = 0; index < velocityPointCount; index++)
	{
		if(index > uppers[tIndex])
			tIndex++;
		double deltaX = (x1[index]-x2[index])*vScale[index].x;
		double deltaY = (y1[index]-y2[index])*vScale[index].y;
		double deltaT = deltaTs[tIndex];
		double resSquared = deltaX*deltaX + deltaY*deltaY;
		correlationLocationRes[index] = sqrt(resSquared);
		correlationVelocityRes[index] = correlationLocationRes[index]/deltaT;
		rmsLocationRes += resSquared;
		rmsVelocityRes += resSquared/(deltaT*deltaT);
	}

	rmsLocationRes = sqrt(rmsLocationRes/(double)velocityPointCount);
	rmsVelocityRes = sqrt(rmsVelocityRes/(double)velocityPointCount);
	printf("  RMS correlation location uncertainty: %g\n", rmsLocationRes);
	printf("  RMS correlation velocity uncertainty: %g\n", rmsVelocityRes);
	
	hid_t fileID = H5Fopen(combinedTiePointsFileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ACCIVPass::appendCorrelationUncertainties: Could not open file %s\n",
                  (const char *)combinedTiePointsFileName);
		return false;
	}

	UString dataSetName = "/rmsCorrelationLocationUncertainty";
	hsize_t dimensions[1];
	dimensions[0] = 1;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, &rmsLocationRes);
	if((status < 0))
	{
		fprintf(stderr, "ACCIVPass::appendCorrelationUncertainties: Could not make dataset %s in file %s\n",
                  (const char *)dataSetName, (const char *)combinedTiePointsFileName);
		return false;
	}

	dataSetName = "/rmsCorrelationVelocityUncertainty";
	dimensions[0] = 1;
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, &rmsVelocityRes);
	if((status < 0))
	{
		fprintf(stderr, "ACCIVPass::appendCorrelationUncertainties: Could not make dataset %s in file %s\n",
                  (const char *)dataSetName, (const char *)combinedTiePointsFileName);
		return false;
	}

	dimensions[0] = (hsize_t)velocityPointCount;
	dataSetName = "/correlationLocationResiduals";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, correlationLocationRes.getData());
	if((status < 0))
	{
		fprintf(stderr, "CCIVPass::appendCorrelationUncertainties: Could not make dataset %s in file %s\n", 
			(const char *)dataSetName, (const char *)combinedTiePointsFileName);
		return false;
	}

	dimensions[0] = (hsize_t)velocityPointCount;
	dataSetName = "/correlationVelocityResiduals";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, correlationVelocityRes.getData());
	if((status < 0))
	{
		fprintf(stderr, "CCIVPass::appendCorrelationUncertainties: Could not make dataset %s in file %s\n", 
			(const char *)dataSetName, (const char *)combinedTiePointsFileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}


bool ACCIVPass::constructVelocityAlongCurvedPaths(const UString & tiePointsFileName,
	const UString & gridVelocityFileName, const UString & scatteredVelocityFileName)
{

	if(!overwriteWorkFiles)
	{
		if(hdf5FileExists(scatteredVelocityFileName))
		{
			return true;
		}
	}

	VectorField2D gridVelocity;
	if(!gridVelocity.read(gridVelocityFileName))
		return false;

	ScatteredVectorData2D scatteredVelocity;

	TiePointSet tiePoints;
	if(!tiePoints.read(tiePointsFileName))
		return false;

	if(!tiePoints.computeCurvedPathVelocities(gridVelocity, pathVectorCount, advectionErrorTolerance,
		advectionMaxTimeStepCount, scatteredVelocity))
		return false;

	UArray<double> & coeffs = tiePoints.getCorrelationCoefficients();
	UArray<double> & deltaTs = tiePoints.getDeltaTs();
	UArray<SInt32> & lowers = tiePoints.getLowerIndexDeltaT();
	UArray<SInt32> & uppers = tiePoints.getUpperIndexDeltaT();
	for(SInt32 timeIndex = 0; timeIndex < deltaTs.getSize(); timeIndex++)
	{
		double deltaT = deltaTs[timeIndex];
		for(SInt32 tiePointIndex = lowers[timeIndex]; tiePointIndex <= uppers[timeIndex]; tiePointIndex++)
		{
			for(SInt32 pathIndex = 0; pathIndex < pathVectorCount; pathIndex++)
			{
				SInt32 vIndex = tiePointIndex*pathVectorCount + pathIndex;
				scatteredVelocity.getWeights(vIndex) = coeffs[tiePointIndex]*deltaT/maxSeparationTime;
			}
		}
	}


	if(!scatteredVelocity.write(scatteredVelocityFileName))
		return false;

	if(!appendScaledScatteredVelocity(scatteredVelocityFileName))
		return false;

	return true;
}

bool ACCIVPass::removeOutliers(const UArray<vector2d> & residuals,
								const UString & inTiePointsFileName,
								const UString & outTiePointsFileName)
{
	if(!overwriteWorkFiles)
	{
		if(hdf5FileExists(outTiePointsFileName))
		{
			return true;
		}
	}

	UArray<double> resMag(residuals.getSize());
	for(SInt32 index = 0; index < resMag.getSize(); index++)
	{
		resMag[index] = residuals[index].length();
	}


	double median;
	{
		// find the median residual value
		// quick sort the residuals
		UArray<double> residualCopy = resMag;
		residualCopy.quickSort();

		// the median is the value at the middle of the list
		median = residualCopy[residualCopy.getSize()/2];
	}

	// read in the old tie points
	TiePointSet oldTiePoints, newTiePoints;
	if(!oldTiePoints.read(inTiePointsFileName))
		return false;
	
	SInt32 removedCount = 0;
	// for each value of deltaT
	for(SInt32 timeIndex = 0; timeIndex < oldTiePoints.getDeltaTs().getSize(); timeIndex++)
	{
		double deltaT = oldTiePoints.getDeltaTs()[timeIndex];
		SInt32 lower = oldTiePoints.getLowerIndexDeltaT()[timeIndex];
		SInt32 upper = oldTiePoints.getUpperIndexDeltaT()[timeIndex];
		// add this deltaT to the new tie points
		newTiePoints.getDeltaTs().add(deltaT);
	    // add this lowerIndexDeltaT to the new tie points
		newTiePoints.getLowerIndexDeltaT().add(newTiePoints.getX1().getSize());
		// for each tie point
		for(SInt32 tiePointIndex = lower; tiePointIndex <= upper; tiePointIndex++)
		{
			bool addTiePoint = true;
			// for each velocity value corresponding to this tie point
			for(SInt32 pathIndex = 0; pathIndex < pathVectorCount; pathIndex++)
			{
				SInt32 velocityIndex = pathIndex + pathVectorCount*tiePointIndex;
				if(resMag[velocityIndex] > median*outlierThresholdConstant)
				{
					addTiePoint = false;
					removedCount++;
					break;
				}
			}

			// if any of the velocities on this path was an outlier, skip this tie point
			if(addTiePoint)
			{
				newTiePoints.getX1().add(oldTiePoints.getX1()[tiePointIndex]);
				newTiePoints.getY1().add(oldTiePoints.getY1()[tiePointIndex]);
				newTiePoints.getX2().add(oldTiePoints.getX2()[tiePointIndex]);
				newTiePoints.getY2().add(oldTiePoints.getY2()[tiePointIndex]);
				newTiePoints.getCorrelationCoefficients().add(oldTiePoints.getCorrelationCoefficients()[tiePointIndex]);
			}
		}
	    // add this upperIndexDeltaT to the new tie points
		newTiePoints.getUpperIndexDeltaT().add(newTiePoints.getX1().getSize()-1);
	}
	double removedFraction = removedCount/(double)oldTiePoints.getX1().getSize();
	printf("Outliers removed, amounting to %.2f%% of total tie points.\n", removedFraction*100);
	// save the new tie points
	if(!newTiePoints.write(outTiePointsFileName))
		return false;
	
	return true;
}


bool ACCIVPass::hdf5FileExists(const UString & fileName)
{
	/* Save old error handler */
	//H5E_auto2_t old_func;
	//herr_t (*old_func)(void*);
	//void *old_client_data;

	//H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

	/* Turn off error handling */
	//H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	hid_t fileID;
	H5E_BEGIN_TRY {
		fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	} H5E_END_TRY;

	/* Restore previous error handler */
	//H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

	if(fileID < 0)
	{
		return false;
	}

	herr_t status = H5Fclose(fileID);
	return true;
}

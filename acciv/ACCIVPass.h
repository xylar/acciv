#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>
#include "VectorField2D.h"

class ACCIVPass
{
public:
	ACCIVPass()
	{
	}

	~ACCIVPass(void)
	{
	}

	bool doPass(const UString & inFolder);



private:
	bool getParameters();
	bool setUpWorkDirectory();
	void setUpTiePointAndVelocityTemplates();
	bool getGridInfo();
	bool getImageTimes();
	bool getImagePairs();
	void computeAdvectionTime();
	bool doImageAdvection();
	bool correlateImages();
	bool undoTiePointAdvection();
	bool constructVelocities();
	bool combineTiePoints();
	bool constructNoCurvatureVelocity();
	bool constructGridVelocityFromScatteredVelocity(const UString & inOutScatteredVelocityFileName,
		const UString & outGridVelocityFileName,
		UInt32 minBinPoints, UInt32 pathVectorCount,
		UArray<vector2d> & outResiduals, double & meanRes);
	
	bool writeFinalGridVelocity(const UString & inGridVelocityFileName,
		const UString & inScatteredVelocityFileName,
		const UString & inFinalGridVelocityFileName);

	bool appendOneSigmaUncertainty(const UArray<vector2d> & residuals,
		double meanResidual,
		const UString & inOutScatteredVelocityFileName);
	bool appendCorrelationUncertainties();
	bool appendScaledScatteredVelocity(const UString & inOutScatteredVelocityFileName);
	bool appendScaledGridVelocity(const UString & inOutGridVelocityFileName);

	bool constructVelocityAlongCurvedPaths(const UString & tiePointsFileName,
		const UString & gridVelocityFileName, const UString &  scatteredVelocityFileName);
	bool removeOutliers(const UArray<vector2d> & residuals,
		const UString & inTiePointsFileName,
		const UString & outTiePointsFileName);

	bool hdf5FileExists(const UString & fileName);

	UString folder;

	UString gridGeometryFactorsFileName;
	UArray<UString> earlierImageFileNames;
	UArray<UString> laterImageFileNames;
	UString outGridVelocityFileName;
	UString outScatteredVelocityFileName;

	bool advectImages;
	UString advectionTimeFlag;
	UArray<UString> inGridVelocityFileNames;
	UArray<SInt32> earlierImageIndices;
	UArray<SInt32> laterImageIndices;

	UArray<SInt32> correlationBoxSize;
	UArray<SInt32> searchRange;
	UArray<SInt32> stride;
	double correlationTolerance;
	double maxSeparationTime;

	double minimumCorrelationCoefficient;

	double minimumTimeSeparation;
	double maximumTimeSeparation;

	bool removeScatteredOutliers;
	SInt32 streamlineFollowingSmoothFitIterationCount;
	SInt32 smoothFitMinControlPointScatteredNeighbors;
	SInt32 pathVectorCount;
	double advectionErrorTolerance;
	SInt32 advectionMaxTimeStepCount;
	double outlierThresholdConstant;

	bool maskFinalVelocity;

	double advectionTime;

	UString workDirectory;

	bool doImageAdvectionStep, doImageCorrelationStep, 
		doUndoTiePointAdvectionStep, doVelocityConstructionStep;

	bool overwriteWorkFiles;

	UArray<SInt32> imageSize;
	UArray<double> imageBounds;

	UArray<double> earlierImageTimes;
	UArray<double> laterImageTimes;
	UArray<SInt32> imagePairs;

	UString earlierAdvectedImageTemplate,
		laterAdvectedImageTemplate;

	UString correlationTiePointTemplate,
		deadvectedTiePointFileName,
		combinedCorrelationTiePointsFileName,
		combinedTiePointsFileName;

	UString noCurvedPathScatteredVelocityFileName, 
		noCurvedPathGridVelocityFileName;

	UString curvedPathScatteredVelocityTemplate,
		curvedPathGridVelocityTemplate,
		noOutlierScatteredVelocityTemplate,
		noOutlierGridVelocityTemplate,
		noOutlierTiePointsTemplate;

	VectorField2D gridVelocityScaling;
};

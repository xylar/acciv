#include "ScatteredVectorData2D.h"

//#include "VectorField2D.h"
//#include <tools/URandom.h>

#include <hdf5.h>
#include <hdf5_hl.h>

/*
#include <math.h>
#include <stdio.h>
#include <algorithm>
*/


bool ScatteredVectorData2D::read(UString fileName)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Could not open file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/dataX";
	hsize_t dimensions[1];
	H5T_class_t typeClass;
	size_t typeSize;
	herr_t status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8))
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	SInt32 size = (SInt32)dimensions[0];
	setSize(size);

	UArray<double> data(size);

	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, data.getData());
	if(status < 0)
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		getData(index).x = data[index];
	}

	dataSetName = "/dataY";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, data.getData());
	if(status < 0)
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		getData(index).y = data[index];
	}

	dataSetName = "/x";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, getX().getData());
	if(status < 0)
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, getY().getData());
	if(status < 0)
	{
		fprintf(stderr, "ScatteredVectorData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/weights";
	bool weightsFound = H5LTfind_dataset(fileID, dataSetName);
	if(!weightsFound)
	{
		dataSetName.makeUpper();
		weightsFound = H5LTfind_dataset(fileID, dataSetName);
	}
	if(weightsFound)
	{

		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
		if(status < 0)
		{
			dataSetName.makeUpper();
			status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
		}
		if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
		{
			fprintf(stderr, "ScatteredData2D::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
			return false;
		}
		status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, weights.getData());
		if(status < 0)
		{
			fprintf(stderr, "ScatteredData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
			return false;
		}
	}
	else
	{
		for(UInt32 index = 0; index < getSize(); index++)
		{
			weights[index] = 1.0;
		}
	}

	status = H5Fclose(fileID);

	return true;
}

bool ScatteredVectorData2D::write(UString fileName)
{
	hid_t fileID = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ScatteredVectorData2D::write: Could not create file %s\n", (const char *)fileName);
		return false;
	}

	UArray<double> data(size);
	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		data[index] = getData(index).x;
	}

	UString dataSetName = "/dataX";
	hsize_t dimensions[1];
	dimensions[0] = (hsize_t)getSize();
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, data.getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredVectorData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		data[index] = getData(index).y;
	}

	dataSetName = "/dataY";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, data.getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredVectorData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}


	dataSetName = "/x";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, getX().getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredVectorData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, getY().getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredVectorData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/weights";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, getWeights().getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredVectorData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

/*
void ScatteredVectorData2D::smoothFit(VectorField2D & vField, UInt32 minBinPoints, double minVariance, double & outMeanRes, 
									  UArray<double> & outSmoothnessX, UArray<double> & outSmoothnessY, double testDataFraction) const
{
	// field is the x component of the vector field, used to get the x and y info (the y component contains the same x and y info)
	//ScalarField2D & field = vField.getDataX();

	// the x and y scattered data (these are the same for the x and y components, we chose the x component here)
	const UArray<double> & x = getDataX().getX();
	const UArray<double> & y = getDataX().getY();
	// the bounds of the grid
	double xl = vField.getXLower();
	double xu = vField.getXUpper();
	double yl = vField.getYLower();
	double yu = vField.getYUpper();

	// the bounds of the scattered data
	double maxX = x[0];
	double minX = x[0];
	double maxY = y[0];
	double minY = y[0];
	for(UInt32 index = 1; index < getDataX().getSize(); index++)
	{
		maxX = std::max(maxX,x[index]);
		minX = std::min(minX,x[index]);
		maxY = std::max(maxY,y[index]);
		minY = std::min(minY,y[index]);
	}

	// the size of the grid in pixels
	UInt32 nx = vField.getXSize();
	UInt32 ny = vField.getYSize();
	// the number of scattered data points
	UInt32 nPoints = getDataX().getData().getSize();

	// the area of a rectangle containing all the scattered data
	double area = (maxX-minX)*(maxY-minY);
	// the physical area of the image
	double imageArea = fabs((xu - xl)*(yu-yl));
	// the area of the scattered data rectangle in pixels
	double pixelArea = (nx-1)*(ny-1)*area/imageArea;

	// x.getSize() is the number of scattered points.
	// x.getSize()/minBinPoints is an approximate guess at the number of bins needed to divide the scattered
	//  data into bins, each containing minBinPoints data points.
	// pixelArea/(x.getSize()/minBinPoints) is the approximate number of pixels in each such bin
	// sqrt(pixelArea/(x.getSize()/minBinPoints)) is the number of pixels on a side for a square bin of pixels
	// goalBinPixels 8 times this size on each side, 64 times the area, and should therefore contain,
	//  on average 64*minBinPoints data points, plenty for good statistics
	double goalBinCount = x.getSize()/(double)minBinPoints;
	double goalBinPixels = 8*sqrt(pixelArea/goalBinCount);
    // the number of bins in x and y, the width or height of the grid in pixels (nx or ny) divided by the number of pixels on a
	//  side for a square bin, rounded; at least 2 bins are required
	UInt32 nxBin = std::max((UInt32)(nx/goalBinPixels + 0.5),(UInt32)2);
	UInt32 nyBin = std::max((UInt32)(ny/goalBinPixels + 0.5),(UInt32)2);

	UArray<UInt32> testIndices((UInt32)(testDataFraction*nPoints+0.5));
	for(UInt32 index = 0; index < testIndices.getSize(); index++)
	{
		// Withhold a fraction of the scattered data, selected at random, to test the quality of the fit.
		// Terminate the fitting process when the fit to the test data stops improving
		testIndices[index] = (UInt32)(URand()*nPoints+0.5);
		testIndices[index] = std::min(testIndices[index], nPoints-1);
	}
	testIndices.quickSort();
	testIndices.uniq(); // remove any duplicates

	// fit indices are the scattered data indices not being withheld for testing
	UInt32 fitCount = nPoints-testIndices.getSize();
	UArray<UInt32> fitIndices(fitCount);
	UArray<double> fitX(fitCount),fitY(fitCount);
	UInt32 testIndex = 0;
	UInt32 index = 0;
	for(UInt32 fitIndex = 0; fitIndex < fitCount; fitIndex++)
	{
		while(testIndex < testIndices.getSize() && index == testIndices[testIndex])
		{
			testIndex++;
			index++;
		}
		fitIndices[fitIndex] = index;
		fitX[fitIndex] = x[index];
		fitY[fitIndex] = y[index];
		index++;
	}

	printf("  smooth fit using %.2f%% of data for testing and %.2f%% for fitting\n",
			100*testIndices.getSize()/(double)nPoints,100*fitCount/(double)nPoints);

	//printf("  %i %i %i\n", testIndices.getSize(), fitCount, nPoints);

	// the first level of the quad tree is used to divide the data into nxBins by nyBins bins with bounds xl <= x <= xu
	// and similarly for y.  The x component of the scattered data getDataX() is passed to the quad tree, since only
	// the x and y data are needed (the y component could also have been passed since thye contain the same x and y data)
	QuadTreeLevel level(fitX, fitY, nxBin, nyBin, xl, xu, yl, yu);

	// the maximum number of levels is the approximate number of times that a square bin of
	//  goalBinPixels on a side would need to be bisected until there were 2 pixels on a side;
	//  maxLevels must be at least 1
	UInt32 maxLevels = std::max((UInt32)(log(goalBinPixels)/log(2.0) + 0.5),(UInt32)1);
	//maxLevels = 2;

	// Initially, the smooth fit is all zeros and the "residual" that we would like to fit
	//  is just a copy of the scattered vector data x and y components, stored in fxRes and fyRes
	UArray<double> fxRes(getDataX().getData());
	UArray<double> fyRes(getDataY().getData());

	// the size of the map of smoothing lengths should be nx by ny
	outSmoothnessX.setSize(nx*ny);
	outSmoothnessY.setSize(nx*ny);
	for(UInt32 index = 0; index < nx*ny; index++)
	{
		// the initial smoothing length in each direction is the size of bins at the coarsest level in pixels
		outSmoothnessX[index] = (nx-1)/(double)level.nx;
		outSmoothnessY[index] = (ny-1)/(double)level.ny;
	}

	double testRes = 1e30;
	// each level of the quad tree is a factor of 2 in each direction
	//  higher resolution than the previous level
	for(UInt32 levelIndex = 0; levelIndex < maxLevels; levelIndex++)
	{
		printf("   smooth fit level: %i, %i x %i bins with pixel smoothness %.1f x %.1f\n", levelIndex,
				level.nx, level.ny, (nx-1)/(double)level.nx, (ny-1)/(double)level.ny);
		bool result = smoothFitLevel(level, minBinPoints,
			minVariance, nx, ny, nPoints,
			x, y, xl, xu, yl, yu,
			testIndices, fitIndices,
			fxRes, fyRes, vField,
			outSmoothnessX, outSmoothnessY,
			outMeanRes, testRes);

		if(!result)
			return;

		if(levelIndex == maxLevels-1)
		{
			printf("  smooth fit finished: maximum levels reached\n");
		}
		else
		{
			// if this was not the last level, bisect each bin to get a new set of bins with twice the resolution
			level.bisect();
			printf("   fit RMS residual at this level: %g\n", outMeanRes);
		}
	}
}

bool ScatteredVectorData2D::smoothFitLevel(const QuadTreeLevel & level, UInt32 minBinPoints,
		double minVariance, UInt32 nx, UInt32 ny, UInt32 nPoints,
		const UArray<double> & x, const UArray<double> & y, double xl, double xu, double yl, double yu,
		const UArray<UInt32> testIndices, const UArray<UInt32> fitIndices,
		UArray<double> & fxRes, UArray<double> & fyRes, VectorField2D & field,
		UArray<double> & outSmoothnessX, UArray<double> & outSmoothnessY,
		double & outMeanRes, double & inOutTestRes) const
{
	// have we found any fits at this level?
	bool fitFound = false;
	// have there been any bilinear fits with enough points to perform a fit at this level?
	bool enoughPoints = false;

	// storage for the bilinear coefficients at this level.  There are 4 coefficients for each vector
	//  component for each bilinear fit, and there will be nx-1 times ny-1 fits (the number of groups
	//  of 4 adjacent bins in the image) plus a border of zeros around the edge of the domain (the vector
	//  field is required to go to zero at the edges).
	UArray<double> coeffs(8*(level.nx+1)*(level.ny+1));
	// similar coefficients used to find the smoothing length maps
	UArray<double> levelCoeffs(coeffs.getSize());

	// loop over bins in y, not including the first and last
	for(UInt32 yIndex = 1; yIndex < level.ny; yIndex++)
	{
		// loop over bins in x, not including the first and last
		for(UInt32 xIndex = 1; xIndex < level.nx; xIndex++)
		{
			// there are 8 coefficient per bin, so the current set of 8 starts
			//  at index given by coeffOffset.
			UInt32 coeffOffset = 8*(xIndex + (level.nx+1)*yIndex);
			// zero out the coefficients, the default in case we don't find a fit
			for(UInt32 index = 0; index < 8; index++)
			{
				coeffs[index+coeffOffset] = 0.0;
			}
			// create empty lists of indices, x's, y's, fx's and fy's that will
			//  be used to store the data in the 4 adjacent bins surrounding the
			//  current location
			UArray<UInt32> indices;
			UArray<double> xLocal,yLocal,fxLocal,fyLocal;
			// loop over the 4 adjacent bins
			for(UInt32 yOffset = 0; yOffset < 2; yOffset++)
			{
				for(UInt32 xOffset = 0; xOffset < 2; xOffset++)
				{
					// the current bin index (the same as a node in the QuadTreeLevel)
					UInt32 nodeIndex = (xIndex+xOffset-1)+level.nx*(yIndex+yOffset-1);
					// the current bin (or node)
					const QuadTreeNode & node = level.nodes[nodeIndex];
					// add all the indices of scattered data from this bin to the list
					// indices.appendArray(node.indices);
					for(UInt32 index = 0; index < node.x.getSize(); index++)
					{
						// add the location of the points from this bin to the list
						// The points have values from -0.5 to 0.5 in the bin, but we want
						// them to have values from -0.5 to 0.5 in the group of 4 bins, so
						// we scale them appropriately
						xLocal.add(0.5*(node.x[index]+xOffset-1));
						yLocal.add(0.5*(node.y[index]+yOffset-1));
						indices.add(fitIndices[node.indices[index]]);
					}
				}
			}
			// do we have enough points?
			if(indices.getSize() < minBinPoints)
				// no, skip to the next bin
				continue;

			// yep, we've found at least one bin with enough points to do a fit
			enoughPoints = true;

			// copy the scattered vector data from this group of 4 bins
			fxLocal.setSize(indices.getSize());
			fyLocal.setSize(indices.getSize());
			for(UInt32 index = 0; index < indices.getSize(); index++)
			{
				fxLocal[index] = fxRes[indices[index]];
				fyLocal[index] = fyRes[indices[index]];
			}

			// do a bilinear fit to the data in this group of 4 bins
			if(getBilinearCoeffs(xLocal,yLocal,fxLocal,fyLocal,minVariance,
				coeffs.getData() + coeffOffset))
			{
				// we found a fit, so the coefficients for the smoothing length mask will be (1,0,0,0,0,0,0,0)
				levelCoeffs[coeffOffset] = 1.0;
				// we found at least one fit
				fitFound = true;
			}
		}
	}
	if(!enoughPoints)
	{
		// no bins had enough points to perform a fit at this level.  We're done.
		printf("  smooth fit finished: no fits were found at this level\n");
		return false; // no fits were found, so we're done
	}
	if(fitFound)
	{
		// we found a fit at least once, so we should modify the residuals and the fit function

		// We had a set of bilinear fit coefficients on a rectangular grid of size level.nx+1 by level.ny+1.
		// Conceptually, we can also think of these as being a set of coefficients for a smooth bicubic spline fit to
		// the scattered data residual.  Each bicubic patch has 4 corners, with coefficients found by a different, overlapping
		// bilinear fit.

		// First, we evaluate the cubic spline fit to the residual at the scattered data points
		UArray<double> fxScatteredFit(nPoints),fyScatteredFit(nPoints);
		evaluateScatteredCubic(coeffs,x,y,xl,xu,yl,yu,level.nx,level.ny,fxScatteredFit,fyScatteredFit);

		// subtract the cubic fit to the residual from the previous residual at the scattered points
		outMeanRes = 0.0;
		for(UInt32 index = 0; index < nPoints; index++)
		{
			fxScatteredFit[index] = fxRes[index] - fxScatteredFit[index];
			fyScatteredFit[index] = fyRes[index] - fyScatteredFit[index];
			outMeanRes += fxScatteredFit[index]*fxScatteredFit[index] + fyScatteredFit[index]*fyScatteredFit[index];
		}
		// compute the root-mean-squared residual at this level, overwriting the RMS residual from the last level
		outMeanRes = sqrt(outMeanRes/nPoints);

		if(testIndices.getSize() != 0)
		{
			double testRes = 0.0;
			for(UInt32 index = 0; index < testIndices.getSize(); index++)
			{
				UInt32 testIndex = testIndices[index];
				if(testIndex >= fxScatteredFit.getSize())
					double a = 5;
				testRes += fxScatteredFit[testIndex]*fxScatteredFit[testIndex] + fyScatteredFit[testIndex]*fyScatteredFit[testIndex];
			}
			// compute the root-mean-squared residual at this level, overwriting the RMS residual from the last level
			testRes = sqrt(testRes/testIndices.getSize());
			printf("   RMS residual of the fit evaluated at test points: %g\n", testRes);
			if(testRes >= inOutTestRes)
			{
				printf("  smooth fit finished: no improvement in the test-point RMS residual at this level\n");
				return false; // no improvement in the fit to the test data, so we're done
			}
			inOutTestRes = testRes;
		}

		// we're keeping the fit so copy it to the scattered residuals
		for(UInt32 index = 0; index < nPoints; index++)
		{
			fxRes[index] = fxScatteredFit[index];
			fyRes[index] = fyScatteredFit[index];
		}

		// Next, we evaluate the bicubic spline fit on the regular grid.  Since the fit is to the *residuals*
		// of the scattered data (below we will subtract the current fit from the scattered data to update the residuals),
		// we add the fit at the current level to the fits from all previous levels.  The reason that we take this approach
		// is that the fit to the residuals should be zero where there is not enough data to make a fit. (It would be tricky
		// to figure out what the fit should be in cases of not enough data if we were using the full function and not the
		// residual.)

		// temporary data to store the two vector components of the grid fit to the current residuals
		UArray<double> fxGridFit(nx*ny),fyGridFit(nx*ny);
		// evaluate the bicubic function with coefficients given in coeffs on a regular grid.
		// The coefficients are on a grid of size level.nx by level.ny, while the evaluation
		// grid is of size nx by ny
		evaluateGridCubic(coeffs,level.nx,level.ny,nx,ny,fxGridFit,fyGridFit);

		// fx and fy are references to the data for the smooth fit output.
		//  Initially, they are all zeros, and a new fit is added to them at each level.
		for(UInt32 index = 0; index < nx*ny; index++)
		{
			// add the fit to the scattered data residuals at the current level to the overall
			//  grid fit
			field[index].x += fxGridFit[index];
			field[index].y += fyGridFit[index];
		}
		// the mask for the smoothing lengths at this level is found from the levelCoeffs
		// The y component is irrelevant, and is just there because we wanted to reuse the same
		// evaluateGridCubic function as above
		// fxGridFit contains a mask between 0 and 1 that says what fraction of the fit at
		// this level came from an actual bilinear fit (as opposed to zeros because there was not
		// enough data of a fit).
		evaluateGridCubic(levelCoeffs,level.nx,level.ny,nx,ny,fxGridFit,fyGridFit);
		for(UInt32 index = 0; index < nx*ny; index++)
		{
			outSmoothnessX[index] = (1.0 - fxGridFit[index])*outSmoothnessX[index] + fxGridFit[index]*(nx-1)/(double)level.nx;
			outSmoothnessY[index] = (1.0 - fxGridFit[index])*outSmoothnessY[index] + fxGridFit[index]*(ny-1)/(double)level.ny;
		}

	}

	return true;
}

bool ScatteredVectorData2D::getBilinearCoeffs(const UArray<double> & xs, const UArray<double> & ys,
	const UArray<double> & fx, const UArray<double> & fy, double minVar, double * outCoeffs) const
{
	double A[4][4], bx[4], by[4];
	//const double maxDist = 0.2;

	for(UInt32 j = 0; j < 4; j++)
	{
		for(UInt32 i = 0; i < 4; i++)
		{
			A[i][j] = 0.0;
		}
		bx[j] = 0.0;
		by[j] = 0.0;
	}

	UInt32 nPoints = xs.getSize();

	double meanX = 0.0, meanY = 0.0;
	double varXX = 0.0, varXY = 0.0, varYY = 0.0;
	for(UInt32 index = 0; index < nPoints; index++)
	{
		double x = xs[index];
		double y = ys[index];
		meanX += x;
		meanY += y;
		varXX += x*x;
		varXY += x*y;
		varYY += y*y;
	}
	A[0][0] = nPoints;
	A[0][1] = meanX;
	A[0][2] = meanY;

	A[1][1] = varXX;
	A[0][3] = varXY;
	A[1][2] = varXY;
	A[2][2] = varYY;

	meanX /= nPoints;
	meanY /= nPoints;

	varXX = varXX/nPoints-meanX*meanX;
	varXY = varXY/nPoints-meanX*meanY;
	varYY = varYY/nPoints-meanY*meanY;

	// eigenvalues of the covariance matrix
	// [varXX-lambda varXY
	//  varXY        varYY-lambda] = 0
	// (varXX-lambda)*(varYY-lambda) - varXY^2 = 0
	// lambda^2 - (varXX+varYY)*lambda + (varXX*varYY-varXY^2) = 0
	//double a = 1.0;
	double b = -(varXX+varYY);
	double deltaVar = varXX-varYY;
	double radicand = sqrt(deltaVar*deltaVar+varXY*varXY);
	double eig1 = 0.5*(-b + radicand);
	double eig2 = 0.5*(-b - radicand);
	double minEig = std::min(eig1,eig2);
	if(minEig < minVar)
		return false; // the points appear to be colinear

	double fxMin = fx[0], fxMax = fx[0], fyMin = fy[0], fyMax = fy[0];
	for(UInt32 index = 0; index < nPoints; index++)
	{
		double x = xs[index];
		double y = ys[index];
		double fxLocal = fx[index];
		double fyLocal = fy[index];
		fxMin = std::min(fxMin,fxLocal);
		fxMax = std::max(fxMax,fxLocal);
		fyMin = std::min(fyMin,fyLocal);
		fyMax = std::max(fyMax,fyLocal);
		A[1][3] += x*x*y;
		A[2][3] += x*y*y;
		A[3][3] += x*x*y*y;
		bx[0] += fxLocal;
		bx[1] += x*fxLocal;
		bx[2] += y*fxLocal;
		bx[3] += x*y*fxLocal;
		by[0] += fyLocal;
		by[1] += x*fyLocal;
		by[2] += y*fyLocal;
		by[3] += x*y*fyLocal;
	}

	for(UInt32 j = 0; j < 4; j++)
	{
		for(UInt32 i = j+1; i < 4; i++)
		{
			A[i][j] = A[j][i];
		}
	}
	int n=4, indx[4];
	double ACopy[4][4];
		
	for(UInt32 j = 0; j < 4; j++)
		for(UInt32 i = 0; i < 4; i++)
			ACopy[i][j] = A[i][j];

	legs(A,n,bx,outCoeffs,indx);
	for(UInt32 i=0; i<4; i++)
		outCoeffs[i] = std::max(std::min(outCoeffs[i],fxMax),fxMin);
	legs(ACopy,n,by,outCoeffs+4,indx);
	for(UInt32 i=0; i<4; i++)
		outCoeffs[i+4] = std::max(std::min(outCoeffs[i+4],fyMax),fyMin);

	// Make sure the coefficients are in the bounds of fxMin and fxMax

	bool areNaNs = false;
	for(UInt32 i=0; i<8; i++)
	{
		if(isnan(outCoeffs[i]))
		{
			areNaNs = true;
			break;
		}
	}

	if(areNaNs)
		for(UInt32 i=0; i<8; i++)
			outCoeffs[i] = 0.0;

	return !areNaNs;
}

void ScatteredVectorData2D::evaluateGridCubic(const UArray<double> & coeffs, UInt32 nxBin, UInt32 nyBin, UInt32 nx, UInt32 ny,
		UArray<double> & outFx, UArray<double> & outFy) const
{
	double scaleX = (nxBin)/(double)(nx-1);
	double scaleY = (nyBin)/(double)(ny-1);
	for(UInt32 yIndex = 0; yIndex < ny; yIndex++)
	{
		for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
		{
			double x = scaleX*xIndex;
			double y = scaleY*yIndex;
			evaluateCubic(x,y,nxBin,nyBin,coeffs,outFx[xIndex+nx*yIndex],outFy[xIndex+nx*yIndex]);
		}
	}
}

void ScatteredVectorData2D::evaluateScatteredCubic(const UArray<double> & coeffs, const UArray<double> & xs,
		const UArray<double> & ys, double xl, double xu, double yl, double yu,
		UInt32 nxBin, UInt32 nyBin,  UArray<double> & outFx, UArray<double> & outFy) const
{
	double scaleX = (nxBin)/(xu-xl);
	double scaleY = (nyBin)/(yu-yl);
	UInt32 nPoints = xs.getSize();
	for(UInt32 pIndex = 0; pIndex < nPoints; pIndex++)
	{
		double x = scaleX*(xs[pIndex]-xl);
		double y = scaleY*(ys[pIndex]-yl);
		evaluateCubic(x,y,nxBin,nyBin,coeffs,outFx[pIndex],outFy[pIndex]);
	}
}

void ScatteredVectorData2D::evaluateCubic(double x, double y, UInt32 nxBin, UInt32 nyBin,
		const UArray<double> & coeffs, double & outFx, double & outFy) const
{
	UInt32 xIndexBin = std::min(nxBin-1,std::max((UInt32)0,(UInt32)(x)));
	UInt32 yIndexBin = std::min(nyBin-1,std::max((UInt32)0,(UInt32)(y)));
	x -= xIndexBin;
	y -= yIndexBin;
	outFx = 0.0;
	outFy = 0.0;
	double bx[2][2], by[2][2];
	bx[0][0] = x*x*(2*x-3)+1;
	bx[0][1] = x*(x*(x-2)+1);
	bx[1][0] = x*x*(-2*x+3);
	bx[1][1] = x*x*(x-1);
	by[0][0] = y*y*(2*y-3)+1;
	by[0][1] = y*(y*(y-2)+1);
	by[1][0] = y*y*(-2*y+3);
	by[1][1] = y*y*(y-1);
	for(UInt32 yo = 0; yo < 2; yo++)
	{
		for(UInt32 xo = 0; xo < 2; xo++)
		{
			for(UInt32 yi = 0; yi < 2; yi++)
			{
				for(UInt32 xi = 0; xi < 2; xi++)
				{
					UInt32 coeffIndex = xi+2*yi+8*(xIndexBin+xo+(nxBin+1)*(yIndexBin+yo));
					double coeff = coeffs[coeffIndex];
					double poly = bx[xo][xi]*by[yo][yi];
					outFx += coeff*poly;
					coeff = coeffs[coeffIndex+4];
					outFy += coeff*poly;
				}
			}
		}
	}
}

*/

#include "MaskedImage.h"

#include <hdf5.h>
#include <hdf5_hl.h>

#include "VectorField2D.h"
#include "ParticleIntegrator.h"

bool MaskedImage::readMaskHDF5File(const UString & fileName)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "MaskedImage::readMaskHDF5File: Could not open file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/mask";
	hsize_t dimensions[2];
	H5T_class_t typeClass;
	size_t typeSize;
	herr_t status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_INTEGER) || (typeSize != 1) || (getYSize() != (SInt32)dimensions[0])
		|| (getXSize() != (SInt32)dimensions[1]))
	{
		fprintf(stderr, "MaskedImage::readMaskHDF5File: Could not get dataset info for %s in file %s\n", 
                  (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	mask.setSize(getXSize()*getYSize());
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_UCHAR, mask.getData());
	if(status < 0)
	{
		fprintf(stderr, "MaskedImage::readMaskHDF5File: Could not read dataset %s in file %s\n", 
                  (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

bool MaskedImage::writeMaskHDF5File(const UString & fileName) const
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "MaskedImage::writeMaskHDF5File: Could not create file %s\n", 
                  (const char *)fileName);
		return false;
	}

	UString dataSetName = "/mask";
	hsize_t dimensions[2];
	dimensions[0] = (hsize_t)getYSize();
	dimensions[1] = (hsize_t)getXSize();
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 2, dimensions, H5T_NATIVE_UCHAR, mask.getData());
	if((status < 0))
	{
		fprintf(stderr, "MaskedImage::writeMaskHDF5File: Could not make dataset %s in file %s\n", 
                  (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

bool MaskedImage::readTimeHDF5File(const UString & fileName, double & outTime)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "MaskedImage::readTimeHDF5File: Could not open file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/time";
	hsize_t dimensions[1];
	H5T_class_t typeClass;
	size_t typeSize;
	dimensions[0] = 1;
	herr_t status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
        if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || (dimensions[0] != 1))
	{
		fprintf(stderr, "MaskedImage::readTimeHDF5File: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, &outTime);
	if(status < 0)
	{
		fprintf(stderr, "MaskedImage::readTimeHDF5File: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

bool MaskedImage::writeTimeHDF5File(const UString & fileName, double inTime)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "MaskedImage::writeMaskHDF5File: Could not create file %s\n",
                  (const char *)fileName);
		return false;
	}

	UString dataSetName = "/time";
	hsize_t dimensions[1];
	dimensions[0] = (hsize_t)1;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, &inTime);
	if((status < 0))
	{
		fprintf(stderr, "MaskedImage::writeMaskHDF5File: Could not make dataset %s in file %s\n",
                  (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

void MaskedImage::maskedInterpolate(const UArray<double> & xs, const UArray<double> & ys, UArray<double> & outValues, UArray<bool> & outMask) const
{
	U_ASSERT((xs.getSize() == ys.getSize()) && (xs.getSize() == outValues.getSize()) && (xs.getSize() == outMask.getSize()));

	static double cubicSplineMatrix[4*4] = 
		{0, 1, 0, 0,
		-0.5, 0, 0.5, 0,
		1, -2.5, 2, -0.5,
		-0.5, 1.5, -1.5, 0.5};

	double xPowers[4];
	double yPowers[4];
	xPowers[0] = 1.0;
	yPowers[0] = 1.0;

	for(SInt32 dataIndex = 0; dataIndex < xs.getSize(); dataIndex++)
	{
		double x = (xs[dataIndex]- getXLower())/getDeltaX();
		double y = (ys[dataIndex]- getYLower())/getDeltaY();
		SInt32 xIndex = (SInt32)x;
		SInt32 yIndex = (SInt32)y;

		// if the index is too close to the edge, extrapolate from the nearest location inbounds
		if((xIndex < 1) || (xIndex > getXSize()-3)
			|| (yIndex < 1) || (yIndex > getYSize()-3))
		{
			outMask[dataIndex] = false;
			continue;
		}

		outMask[dataIndex] = true;
		for(SInt32 j = 0; j < 4; j++)
		{
			for(SInt32 i = 0; i < 4; i++)
			{
				outMask[dataIndex] = outMask[dataIndex] && getMask(xIndex + i - 1,yIndex + j - 1);
			}
		}

		if(!outMask[dataIndex])
			continue;

		x = x - (double)xIndex;
		y = y - (double)yIndex;

		xPowers[1] = x;
		xPowers[2] = x*x;
		xPowers[3] = x*x*x;
		yPowers[1] = y;
		yPowers[2] = y*y;
		yPowers[3] = y*y*y;

		double result = 0;
		for(SInt32 k = 0; k < 4; k++)
		{
			double localResult = 0;
			for(SInt32 j = 0; j < 4; j++)
				for(SInt32 i = 0; i < 4; i++)
					localResult = localResult + xPowers[i]*cubicSplineMatrix[j+4*i]*(*this)(xIndex + j - 1, yIndex + k - 1);
			for(SInt32 i = 0; i < 4; i++)
				result = result + yPowers[i]*cubicSplineMatrix[k+4*i]*localResult;
		}

		outValues[dataIndex] = result;
	}	
}

void MaskedImage::maskedInterpolateAtOffsetPixels(double xLower, double yLower, SInt32 xCount, SInt32 yCount,
												  UArray<double> & outValues, UArray<bool> & outMask) const
{
	U_ASSERT((xCount*yCount == outValues.getSize()) && (xCount*yCount == outMask.getSize()));

	static double cubicSplineMatrix[4*4] = 
		{0, 1, 0, 0,
		-0.5, 0, 0.5, 0,
		1, -2.5, 2, -0.5,
		-0.5, 1.5, -1.5, 0.5};

	double xPowers[4];
	double yPowers[4];
	xPowers[0] = 1.0;
	yPowers[0] = 1.0;

	double x = (xLower - getXLower())/getDeltaX();
	double y = (yLower - getYLower())/getDeltaY();

	SInt32 xLowerIndex = (SInt32)x;
	SInt32 yLowerIndex = (SInt32)y;
	x = x - (double)xLowerIndex;
	y = y - (double)yLowerIndex;

	xPowers[1] = x;
	xPowers[2] = x*x;
	xPowers[3] = x*x*x;
	yPowers[1] = y;
	yPowers[2] = y*y;
	yPowers[3] = y*y*y;

	double coefficients[4*4];
	for(SInt32 i = 0; i < 16; i++)
		coefficients[i] = 0;
	for(SInt32 k = 0; k < 4; k++)
	{
		for(SInt32 j = 0; j < 4; j++)
		{
			double result = 0;
			for(SInt32 l = 0; l < 4; l++)
				for(SInt32 i = 0; i < 4; i++)
					result = result + xPowers[i]*cubicSplineMatrix[j+4*i]*yPowers[l]*cubicSplineMatrix[k+4*l];
			coefficients[j+4*k] = result;
		}
	}

	
	SInt32 dataIndex = 0;
	for(SInt32 yIndex = yLowerIndex; yIndex < yLowerIndex+yCount; yIndex++)
	{
		for(SInt32 xIndex = xLowerIndex; xIndex < xLowerIndex+xCount; xIndex++)
		{
			// if the index is too close to the edge
			if((xIndex < 1) || (xIndex > getXSize()-3)
				|| (yIndex < 1) || (yIndex > getYSize()-3))
			{
				outMask[dataIndex] = false;
				dataIndex++;
				continue;
			}

			outMask[dataIndex] = true;
			for(SInt32 j = 0; j < 4; j++)
			{
				for(SInt32 i = 0; i < 4; i++)
				{
					outMask[dataIndex] = outMask[dataIndex] && getMask(xIndex + i - 1,yIndex + j - 1);
				}
			}

			if(!outMask[dataIndex])
			{
				dataIndex++;
				continue;
			}

			double result = 0;
			for(SInt32 k = 0; k < 4; k++)
				for(SInt32 j = 0; j < 4; j++)
					result = result + coefficients[j+4*k]*(*this)(xIndex + j - 1, yIndex + k - 1);

			outValues[dataIndex] = result;
			dataIndex++;
		}
	}
}

bool MaskedImage::advect(const VectorField2D & velocity, double deltaT, 
	double errorTolerance, SInt32 maximumTimeStepCount)
{
	SInt32 imageSize = getXSize()*getYSize();
	UArray<double> x(imageSize), y(imageSize);
	for(SInt32 imageIndex = 0; imageIndex < imageSize; imageIndex++)
	{
		SInt32 yIndex = (imageIndex)/getXSize();
		SInt32 xIndex = (imageIndex) - getXSize()*yIndex;
		x[imageIndex] = getX()[xIndex];
		y[imageIndex] = getY()[yIndex];
	}

	return advect(velocity, deltaT, x, y, errorTolerance, maximumTimeStepCount);
}

bool MaskedImage::advect(const VectorField2D & velocity, double deltaT,
	UArray<double> & inOutX, UArray<double> & inOutY, double errorTolerance, SInt32 maximumTimeStepCount)
{
	if(deltaT == 0.0)
		return true;
	ParticleIntegrator integrator(velocity);
	double nextTimeStep = 0.1*deltaT;
	double minimumTimeStep = fabs(0.1*deltaT/maximumTimeStepCount);

	SInt32 imageSize = getXSize()*getYSize();
	if(inOutX.getSize() != imageSize || inOutY.getSize() != imageSize)
		return false;

	SInt32 maxCount = 1000000; // keeps us from allocating too much memory at a time

	UArray<double> newData(imageSize);
	UArray<bool> newMask(imageSize);

	for(SInt32 imageIndex = 0; imageIndex < imageSize; imageIndex += maxCount)
	{
		SInt32 localCount = std::min(maxCount, imageSize-imageIndex);
		// set up the points, one for each output pixel
		UArray<double> points(2*localCount);

		for(SInt32 index = 0; index < localCount; index++)
		{
			points[2*index] = inOutX[index + imageIndex];
			points[2*index+1] = inOutY[index + imageIndex];
		}

		// integrate the pixels backward in time to their locations in the original image
		if(!integrator.integrate(points, deltaT, 0, errorTolerance, nextTimeStep, minimumTimeStep, maximumTimeStepCount))
			return false;

		// interpolate the advected image at the points in the original image
		UArray<double> pixelX(localCount), pixelY(localCount), advectedPixels(localCount);
		UArray<bool> advectedMask(localCount);
		for(SInt32 index = 0; index < localCount; index++)
		{
			pixelX[index] = points[2*index];
			pixelY[index] = points[2*index+1];
			inOutX[index + imageIndex] = points[2*index];
			inOutY[index + imageIndex] = points[2*index+1];
		}

		maskedInterpolate(pixelX, pixelY, advectedPixels, advectedMask);

		// copy the advected pixels and mask back into this image
		for(SInt32 index = 0; index < localCount; index++)
		{
			newData[index + imageIndex] = advectedPixels[index];
			newMask[index + imageIndex] = advectedMask[index];
		}
	}
	getData().copy(newData);
	mask.copy(newMask);

	return true;
}

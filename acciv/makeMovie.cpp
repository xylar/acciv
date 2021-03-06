// makeMovie.cpp : generates a sequence of advected/alpha-blended images that can be used to make a movie
//

#include <stdlib.h>

#include "VectorField2D.h"
#include "MaskedImage.h"
#include "MakeDirectory.h"
#include <core/UString.h>

int main(int argc, char * argv[])
{
	if(argc != 6)
	{
		fprintf(stderr, "usage: %s <velocity_file> <earlier_image> <later_image> <delta_t> <output_dir>\n", argv[0]);
		return 0;
	}
	char * velocityFileName = argv[1];
	char * image1FileName = argv[2];
	char * image2FileName = argv[3];
	double deltaT = atoi(argv[4]);
	char * folderName = argv[5];

	const double advectionErrorTolerance = 1e-4;
	const UInt32 advectionMaxTimeStepCount = 10000;

	makeDirectory(folderName);

	VectorField2D velocity;
	if(!velocity.read(velocityFileName))
	{
		fprintf(stderr, "Could not read velocity file %s.\n", velocityFileName);
		return 0;
	}

	MaskedImage image1, image2;
	if(!image1.read(image1FileName, true))
	{
		fprintf(stderr, "Could not read image file %s.\n", image1FileName);
		return 0;
	}
	if(!image2.read(image2FileName, true))
	{
		fprintf(stderr, "Could not read image file %s.\n", image2FileName);
		return 0;
	}

	UInt32 stepCount = UInt32((image2.getTime()-image1.getTime())/deltaT + 0.5) + 1;
	deltaT = (image2.getTime()-image1.getTime())/(stepCount-1);

	UString fileName;
	fileName.format("%s/image_f_%03i.h5",folderName, 0);
	if(!image1.write(fileName))
	{
		fprintf(stderr, "Could not write image file %s.\n", (const char *)fileName);
		return 0;
	}

	SInt32 imageSize = image1.getXSize()*image1.getYSize();
	UArray<double> x(imageSize), y(imageSize);
	for(SInt32 imageIndex = 0; imageIndex < imageSize; imageIndex++)
	{
		SInt32 yIndex = (imageIndex)/image1.getXSize();
		SInt32 xIndex = (imageIndex) - image1.getXSize()*yIndex;
		x[imageIndex] = image1.getX()[xIndex];
		y[imageIndex] = image1.getY()[yIndex];
	}

	for(SInt32 tIndex = 1; tIndex < stepCount; tIndex++)
	{
		MaskedImage advectedImage = image1;
		if(!advectedImage.advect(velocity,deltaT,x,y,advectionErrorTolerance,advectionMaxTimeStepCount))
		{
			fprintf(stderr, "Forward advection failed at frame %i .\n", tIndex);
			return 0;
		}
		fileName.format("%s/image_f_%03i.h5",folderName, tIndex);
		if(!advectedImage.write(fileName))
		{
			fprintf(stderr, "Could not write image file %s.\n", (const char *)fileName);
			return 0;
		}
	}

	for(SInt32 imageIndex = 0; imageIndex < imageSize; imageIndex++)
	{
		SInt32 yIndex = (imageIndex)/image1.getXSize();
		SInt32 xIndex = (imageIndex) - image1.getXSize()*yIndex;
		x[imageIndex] = image1.getX()[xIndex];
		y[imageIndex] = image1.getY()[yIndex];
	}

	fileName.format("%s/image_b_%03i.h5",folderName, stepCount-1);
	if(!image2.write(fileName))
	{
		fprintf(stderr, "Could not write image file %s.\n", (const char *)fileName);
		return 0;
	}
	for(SInt32 tIndex = stepCount-2; tIndex >= 0; tIndex--)
	{
		MaskedImage advectedImage = image2;
		if(!advectedImage.advect(velocity,-deltaT,x,y,advectionErrorTolerance,advectionMaxTimeStepCount))
		{
			fprintf(stderr, "Backward advection failed at frame %i .\n", tIndex);
			return 0;
		}
		fileName.format("%s/image_b_%03i.h5",folderName, tIndex);
		if(!advectedImage.write(fileName))
		{
			fprintf(stderr, "Could not write image file %s.\n", (const char *)fileName);
			return 0;
		}
	}

	for(SInt32 tIndex = 0; tIndex < stepCount; tIndex++)
	{
		fileName.format("%s/image_f_%03i.h5",folderName, tIndex);
		if(!image1.read(fileName))
		{
			fprintf(stderr, "Could not read image file %s.\n", (const char *)fileName);
			return 0;
		}
		fileName.format("%s/image_b_%03i.h5",folderName, tIndex);
		if(!image2.read(fileName))
		{
			fprintf(stderr, "Could not read image file %s.\n", (const char *)fileName);
			return 0;
		}
		double alpha = tIndex/(double)(stepCount-1);
		for(SInt32 dataIndex = 0; dataIndex < image1.getDataSize(); dataIndex++)
		{
			image1.getData()[dataIndex] = (1.0 - alpha)*image1.getData()[dataIndex]
			    + alpha*image2.getData()[dataIndex];
		}
		fileName.format("%s/image_%03i.h5",folderName, tIndex);
		if(!image1.write(fileName))
		{
			fprintf(stderr, "Could not write image file %s.\n", (const char *)fileName);
			return 0;
		}

	}

	return 0;
}


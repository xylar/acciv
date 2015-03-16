#pragma once

#include "MaskedImage.h"
#include <fftw3.h>
#include <complex>

#include "TiePointSet.h"

using namespace std;

class Correlator
{
public:
	Correlator(double inTolerance = 3e-8)
		:fftwRealField(NULL), fftwComplexField(NULL), tolerance(inTolerance)
	{
	}

	~Correlator(void)
	{
	}

	bool correlateImages(const MaskedImage & image1, const MaskedImage & image2,
		UArray<double> & outX1, UArray<double> & outY1, 
		UArray<double> & outX2, UArray<double> & outY2,
		UArray<double> & outCorrelationCoefficients,
		SInt32 boxSizeX, SInt32 boxSizeY, SInt32 searchMinX, SInt32 searchMaxX,
		SInt32 searchMinY, SInt32 searchMaxY, SInt32 strideX = 1, SInt32 strideY = 1);

	bool correlateImages(const MaskedImage & image1, const MaskedImage & image2,
		TiePointSet & outTiePoints,
		SInt32 boxSizeX, SInt32 boxSizeY, SInt32 searchMinX, SInt32 searchMaxX,
		SInt32 searchMinY, SInt32 searchMaxY, SInt32 strideX = 1, SInt32 strideY = 1)
	{
		return correlateImages(image1, image2, outTiePoints.getX1(), outTiePoints.getY1(),
			outTiePoints.getX2(), outTiePoints.getY2(), outTiePoints.getCorrelationCoefficients(),
			boxSizeX, boxSizeY, searchMinX, searchMaxX, searchMinY, searchMaxY, strideX, strideY);
	}


private:
	void setUpFFTW(SInt32 xSize, SInt32 ySize);
	void cleanUpFFTW();

	bool correlateSubImages(const MaskedImage & image1, const MaskedImage & image2,
		SInt32 xIndex, SInt32 yIndex, SInt32 boxSizeX, SInt32 boxSizeY,
		SInt32 searchMinX, SInt32 searchMaxX, SInt32 searchMinY, SInt32 searchMaxY,
		double & outX1, double & outY1, double & outX2, double & outY2,
		double & outCorrelationCoefficient);

	bool findFastCorrelations(const UArray<double> & subImage1, const UArray<double> & mask1,
		const UArray<double> & subImage2, const UArray<double> & mask2, 
		UArray<double> & outCorrelations, UArray<bool> & outCorrelationsMask, 
		SInt32 subImageSizeX, SInt32 subImageSizeY,
		SInt32 correlationsSizeX, SInt32 correlationsSizeY, SInt32 maxCount);

	void transform(const UArray<double> & inField, UArray<complex<double> > & outField_f);
	void crossCorrelate(const UArray<complex<double> > & field1_f, const UArray<complex<double> > & field2_f,
	UArray<double> & result);

	double * fftwRealField;
	fftw_complex * fftwComplexField;

	fftw_plan realToComplexPlan;
	fftw_plan complexToRealPlan;

	SInt32 realFieldSize;
	SInt32 complexFieldSize;

	UArray<complex<double> > subImage1_f, subImage2_f,
		mask1_f, mask2_f, subImage1Squared_f, subImage2Squared_f, product_f;

	UArray<double> correlations, mean1, mean2,
		var1, var2, maskCount, subImage1Squared, subImage2Squared;

	double tolerance;
};

class SlowCorrelationFunctor
{
public:

	SlowCorrelationFunctor(const MaskedImage & inMaskedImage1, SInt32 inXIndex1, SInt32 inYIndex1, const MaskedImage & inMaskedImage2,
		SInt32 inBoxSizeX, SInt32 inBoxSizeY)
		:maskedImage2(inMaskedImage2),
		boxSizeX(inBoxSizeX), boxSizeY(inBoxSizeY),
		image1(inBoxSizeX*inBoxSizeY),mask1(inBoxSizeX*inBoxSizeY)
	{
		for(SInt32 yIndex = 0; yIndex < boxSizeY; yIndex++)
		{
			for(SInt32 xIndex = 0; xIndex < boxSizeX; xIndex++)
			{
				SInt32 subImageIndex = xIndex + boxSizeX*yIndex;
				mask1[subImageIndex] = (double)inMaskedImage1.getMask(inXIndex1 + xIndex,inYIndex1 + yIndex);
				image1[subImageIndex] = inMaskedImage1(inXIndex1 + xIndex,inYIndex1 + yIndex);
			}
		}
	}

	bool operator() (const UArray<double> & x, double & outFunctionValue);

private:
	UArray<double> image1, mask1;
	const MaskedImage & maskedImage2;
	SInt32 x1Lower, y1Lower, boxSizeX, boxSizeY;
};

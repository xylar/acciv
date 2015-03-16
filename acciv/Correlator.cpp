#include "Correlator.h"

#include "SquareBoundedMinimizer.h"


bool Correlator::correlateImages(const MaskedImage & image1, const MaskedImage & image2,
		UArray<double> & outX1, UArray<double> & outY1, 
		UArray<double> & outX2, UArray<double> & outY2,
		UArray<double> & outCorrelationCoefficients,
		SInt32 boxSizeX, SInt32 boxSizeY, SInt32 searchMinX, SInt32 searchMaxX,
		SInt32 searchMinY, SInt32 searchMaxY, SInt32 strideX, SInt32 strideY)
{
	// See if the masks overlap at all.  If not, it's not worth searching
	// (even though the masks *could* possibly overlap if they were shifted).
	bool overlap = false;
	for(SInt32 index=0; index < image1.getDataSize(); index++)
	{
		if((image1.getMask(index)) && (image2.getMask(index)))
		{
			overlap = true;
			break;
		}
	}
	if(!overlap)
	{
		printf("Images don't overlap; they will not be correlated.\n");
		return false;
	}

	SInt32 x1Lower = std::max(0, -searchMinX);
	SInt32 x1Upper = image1.getXSize() - boxSizeX - std::max(0, searchMaxX);
	SInt32 y1Lower = std::max(0, -searchMinY);
	SInt32 y1Upper = image1.getYSize() - boxSizeY - std::max(0, searchMaxY);

	setUpFFTW(2*boxSizeX + searchMaxX - searchMinX, 2*boxSizeY + searchMaxY - searchMinY);

	UArray<double> correlations(boxSizeX*boxSizeY);

	UInt32 percentage = 0;
	for(SInt32 yIndex = y1Lower; yIndex <= y1Upper; yIndex += strideY)
	{
		UInt32 newPercentage = 10*((10*(yIndex-y1Lower))/(y1Upper-y1Lower));
		if(newPercentage > percentage)
		{
			printf("%i%%\n",newPercentage);
			percentage = newPercentage;
		}
		for(SInt32 xIndex = x1Lower; xIndex <= x1Upper; xIndex += strideX)
		{
			double x1, y1, x2, y2, correlationCoefficient;
			bool result = correlateSubImages(image1, image2, xIndex, yIndex,
				boxSizeX, boxSizeY, searchMinX, searchMaxX, searchMinY,
				searchMaxY, x1, y1, x2, y2, correlationCoefficient);
			if(!result)
				continue;
			outX1.add(x1);
			outY1.add(y1);
			outX2.add(x2);
			outY2.add(y2);
			outCorrelationCoefficients.add(correlationCoefficient);

		}
	}

	cleanUpFFTW();

	if(outX1.getSize() == 0)
	{
		printf("No tie points were found.\n");
		return false;
	}

	return true;
}

void Correlator::setUpFFTW(SInt32 xSize, SInt32 ySize)
{
	realFieldSize = xSize*ySize;
	complexFieldSize = (xSize/2 + 1)*ySize;
	fftwRealField = (double *)fftw_malloc(realFieldSize*sizeof(double));
	fftwComplexField = (fftw_complex *)fftw_malloc(complexFieldSize*sizeof(fftw_complex));

	realToComplexPlan = fftw_plan_dft_r2c_2d(ySize, xSize,
		fftwRealField, fftwComplexField, FFTW_MEASURE);
	complexToRealPlan = fftw_plan_dft_c2r_2d(ySize, xSize,
		fftwComplexField, fftwRealField, FFTW_MEASURE);

	subImage1_f.setSize(complexFieldSize);
	subImage2_f.setSize(complexFieldSize);
	mask1_f.setSize(complexFieldSize);
	mask2_f.setSize(complexFieldSize);
	subImage1Squared_f.setSize(complexFieldSize);
	subImage2Squared_f.setSize(complexFieldSize);
	product_f.setSize(complexFieldSize);

	correlations.setSize(realFieldSize);
	mean1.setSize(realFieldSize);
	mean2.setSize(realFieldSize);
	var1.setSize(realFieldSize);
	var2.setSize(realFieldSize);
	maskCount.setSize(realFieldSize);
	subImage1Squared.setSize(realFieldSize);
	subImage2Squared.setSize(realFieldSize);
}

void Correlator::cleanUpFFTW()
{
	fftw_destroy_plan(realToComplexPlan);
	fftw_destroy_plan(complexToRealPlan);

	fftw_free(fftwRealField);
	fftw_free(fftwComplexField);
}


bool Correlator::correlateSubImages(const MaskedImage & image1, const MaskedImage & image2,
	SInt32 inXIndex, SInt32 inYIndex, SInt32 boxSizeX, SInt32 boxSizeY,
	SInt32 searchMinX, SInt32 searchMaxX, SInt32 searchMinY, SInt32 searchMaxY,
	double & outX1, double & outY1, double & outX2, double & outY2,
	double & outCorrelationCoefficient)
{

	SInt32 x1Lower = inXIndex;
	SInt32 x1Upper = inXIndex+boxSizeX-1;
	SInt32 y1Lower = inYIndex;
	SInt32 y1Upper = inYIndex+boxSizeY-1;
	SInt32 x2Lower = inXIndex+searchMinX;
	SInt32 x2Upper = inXIndex+searchMaxX+boxSizeX-1;
	SInt32 y2Lower = inYIndex+searchMinY;
	SInt32 y2Upper = inYIndex+searchMaxY+boxSizeY-1;

	outX1 = (x1Lower + 0.5*(boxSizeX-1))*image1.getDeltaX() + image1.getXLower();
	outY1 = (y1Lower + 0.5*(boxSizeY-1))*image1.getDeltaY() + image1.getYLower();

	SInt32 subImageSizeX = 2*boxSizeX + searchMaxX - searchMinX;
	SInt32 subImageSizeY = 2*boxSizeY + searchMaxY - searchMinY;
	SInt32 subImageSize = subImageSizeX*subImageSizeY;
	UArray<double> subImage1(subImageSize), subImage2(subImageSize), mask1(subImageSize), mask2(subImageSize);

	bool isCompleteSubImageMasked = true;
	SInt32 subYIndex = 0;
	for(SInt32 yIndex = y1Lower; yIndex <= y1Upper; yIndex++)
	{
		SInt32 subXIndex = 0;
		for(SInt32 xIndex = x1Lower; xIndex <= x1Upper; xIndex++)
		{
			SInt32 index = subXIndex + subImageSizeX*subYIndex;
			bool localMask = image1.getMask(xIndex,yIndex);
			isCompleteSubImageMasked = isCompleteSubImageMasked & !localMask;
			mask1[index] = (double)localMask;
			subImage1[index] = mask1[index]*image1(xIndex,yIndex);
			subXIndex++;
		}
		subYIndex++;
	}

	if(isCompleteSubImageMasked)
		return false;

	isCompleteSubImageMasked = true;
	subYIndex = 0;
	for(SInt32 yIndex = y2Lower; yIndex <= y2Upper; yIndex++)
	{
		SInt32 subXIndex = 0;
		for(SInt32 xIndex = x2Lower; xIndex <= x2Upper; xIndex++)
		{
			SInt32 index = subXIndex + subImageSizeX*subYIndex;
			bool localMask = image2.getMask(xIndex,yIndex);
			isCompleteSubImageMasked = isCompleteSubImageMasked & !localMask;
			mask2[index] = (double)localMask;
			subImage2[index] = mask2[index]*image2(xIndex,yIndex);
			subXIndex++;
		}
		subYIndex++;
	}
	if(isCompleteSubImageMasked)
		return false;

	SInt32 correlationsSizeX = searchMaxX - searchMinX + 1;
	SInt32 correlationsSizeY = searchMaxY - searchMinY + 1;

	UArray<double> correlations(correlationsSizeX*correlationsSizeY);
	UArray<bool> correlationsMask(correlationsSizeX*correlationsSizeY);

	bool result = findFastCorrelations(subImage1, mask1, subImage2, mask2, correlations, 
		correlationsMask, subImageSizeX, subImageSizeY, correlationsSizeX, correlationsSizeY,
		boxSizeX*boxSizeY);

	if(!result)
		return false;

	// find the xIndex and yIndex corresponding to the maximum value of the correlation function
	SInt32 maxXIndex = 0, maxYIndex = 0;
	double maxCorrelation = correlations[0];
	for(SInt32 yIndex = 0; yIndex < correlationsSizeY; yIndex++)
	{
		for(SInt32 xIndex = 0; xIndex < correlationsSizeX; xIndex++)
		{
			SInt32 index = xIndex + correlationsSizeX*yIndex;
			if(!correlationsMask[index])
				continue;
			double value = correlations[index];
			if(value > maxCorrelation)
			{
				maxCorrelation = value;
				maxXIndex = xIndex;
				maxYIndex = yIndex;
			}
		}
	}

	if((maxXIndex == 0) || (maxXIndex == correlationsSizeX-1)
		|| (maxYIndex == 0) || (maxYIndex == correlationsSizeY-1))
		return false;

	SInt32 xIndex2 = maxXIndex+inXIndex+searchMinX;
	SInt32 yIndex2 = maxYIndex+inYIndex+searchMinY;

	// use the slow correlation function (with bicubic interpolation) within the minimization routine to find the location of maximum correlation
	SlowCorrelationFunctor functor(image1, inXIndex, inYIndex, image2, boxSizeX, boxSizeY);
	SquareBoundedMinimizer<SlowCorrelationFunctor> minimizer(functor, (double)(xIndex2-1), (double)(xIndex2+1),
		(double)(yIndex2-1), (double)(yIndex2+1), tolerance);
	UArray<double> point(2);
	point[0] = xIndex2;
	point[1] = yIndex2;
	result = minimizer.minimize(point, maxCorrelation);
	if(!result)
		return false;
	maxCorrelation = -maxCorrelation;

	// figure out x2 and y2 from the sub-pixel indices
	outX2 = (point[0] + 0.5*(boxSizeX-1))*image1.getDeltaX() + image1.getXLower();
	outY2 = (point[1] + 0.5*(boxSizeY-1))*image1.getDeltaY() + image1.getYLower();
	outCorrelationCoefficient = maxCorrelation;

	return true;
}

bool Correlator::findFastCorrelations(const UArray<double> & subImage1, const UArray<double> & mask1,
	const UArray<double> & subImage2, const UArray<double> & mask2, 
	UArray<double> & outCorrelations, UArray<bool> & outCorrelationsMask, 
	SInt32 subImageSizeX, SInt32 subImageSizeY,
	SInt32 correlationsSizeX, SInt32 correlationsSizeY, SInt32 maxCount)
{

	UArray<SInt32> subImageIndices(correlationsSizeX*correlationsSizeY);

	for(SInt32 yIndex = 0; yIndex < correlationsSizeY; yIndex++)
	{
		SInt32 subImageIndexY = (subImageSizeY - yIndex)%subImageSizeY;
		for(SInt32 xIndex = 0; xIndex < correlationsSizeX; xIndex++)
		{
			SInt32 subImageIndexX = (subImageSizeX - xIndex)%subImageSizeX;
			subImageIndices[xIndex + correlationsSizeX*yIndex] =
				subImageIndexX + subImageSizeX*subImageIndexY;
		}
	}

	for(SInt32 index = 0; index < realFieldSize; index++)
	{
		subImage1Squared[index] = subImage1[index]*subImage1[index];
		subImage2Squared[index] = subImage2[index]*subImage2[index];
	}

	transform(subImage1,subImage1_f);
	transform(subImage2,subImage2_f);
	transform(mask1,mask1_f);
	transform(mask2,mask2_f);
	transform(subImage1Squared,subImage1Squared_f);
	transform(subImage2Squared,subImage2Squared_f);

	crossCorrelate(subImage1_f, subImage2_f, correlations);
	crossCorrelate(subImage1_f, mask2_f, mean1);
	crossCorrelate(mask1_f, subImage2_f, mean2);
	crossCorrelate(subImage1Squared_f, mask2_f, var1);
	crossCorrelate(mask1_f, subImage2Squared_f, var2);
	crossCorrelate(mask1_f, mask2_f, maskCount);


	bool result = false;
	for(SInt32 correlationIndex = 0; correlationIndex < correlationsSizeX*correlationsSizeY; correlationIndex++)
	{
		SInt32 subImageIndex = subImageIndices[correlationIndex];
		double count = (double)(SInt32)(maskCount[subImageIndex] + 0.5);

		// if more than one quarter of the correlation window is masked out,
		// consider the correlation to be too unreliable
		// (this is why we don't correlate the edges of the image)
		if(count < maxCount)
		{
			outCorrelations[correlationIndex] = -1.0;
			outCorrelationsMask[correlationIndex] = false;
			continue;
		}
		double invCount = 1.0/count;
		double localMean1 = mean1[subImageIndex]*invCount;
		double localMean2 = mean2[subImageIndex]*invCount;
		double localVar1 = var1[subImageIndex]*invCount - localMean1*localMean1;
		double localVar2 = var2[subImageIndex]*invCount - localMean2*localMean2;
		if((localVar1 == 0.0) || (localVar2 == 0.0))
		{
			outCorrelations[correlationIndex] = -1.0;
			outCorrelationsMask[correlationIndex] = false;
			continue;
		}
		outCorrelations[correlationIndex] = (correlations[subImageIndex]*invCount - localMean1*localMean2)/sqrt(localVar1*localVar2);
		outCorrelationsMask[correlationIndex] = true;
		result = true;
	}

	return result;
}

void Correlator::transform(const UArray<double> & inField, UArray<complex<double> > & outField_f)
{
	U_ASSERT((inField.getSize() == realFieldSize) & (outField_f.getSize() == complexFieldSize));

	memcpy(fftwRealField, inField.getData(), realFieldSize*sizeof(double));
	fftw_execute(realToComplexPlan);
	memcpy(outField_f.getData(), fftwComplexField, complexFieldSize*sizeof(fftw_complex));
}

void Correlator::crossCorrelate(const UArray<complex<double> > & field1_f, const UArray<complex<double> > & field2_f,
	UArray<double> & result)
{
	U_ASSERT((field1_f.getSize() == complexFieldSize) & (field1_f.getSize() == complexFieldSize)
		& (result.getSize() == realFieldSize));

	for(SInt32 index = 0; index < complexFieldSize; index++)
		product_f[index] = field1_f[index]*conj(field2_f[index]);
	memcpy(fftwComplexField, product_f.getData(), complexFieldSize*sizeof(fftw_complex));
	fftw_execute(complexToRealPlan);
	memcpy(result.getData(), fftwRealField, realFieldSize*sizeof(double));
	for(SInt32 index=0; index<realFieldSize; index++)
		result[index] = result[index]/(double)realFieldSize;
}

bool SlowCorrelationFunctor::operator() (const UArray<double> & x, double & outFunctionValue)
{
	double inXIndex2 = x[0];
	double inYIndex2 = x[1];

	SInt32 boxSize = boxSizeX*boxSizeY;
	UArray<double> image2(boxSize), mask2(boxSize);
	if((inXIndex2 == (SInt32)inXIndex2) && (inYIndex2 == (SInt32)inYIndex2))
	{
		SInt32 xIndex2 = (SInt32)inXIndex2;
		SInt32 yIndex2 = (SInt32)inYIndex2;
		for(SInt32 yIndex = 0; yIndex < boxSizeY; yIndex++)
		{
			for(SInt32 xIndex = 0; xIndex < boxSizeX; xIndex++)
			{
				SInt32 subImageIndex = xIndex + boxSizeX*yIndex;
				mask2[subImageIndex] = (double)maskedImage2.getMask(xIndex2 + xIndex,yIndex2 + yIndex);
				image2[subImageIndex] = maskedImage2(xIndex2 + xIndex,yIndex2 + yIndex);
			}
		}
	}
	else
	{
		UArray<bool> tempMask2(boxSize);
		double xLower = inXIndex2*maskedImage2.getDeltaX() + maskedImage2.getXLower();
		double yLower = inYIndex2*maskedImage2.getDeltaY() + maskedImage2.getYLower();

		maskedImage2.maskedInterpolateAtOffsetPixels(xLower, yLower, boxSizeX, boxSizeY, image2, tempMask2);
		for(SInt32 index = 0; index < boxSize; index++)
		{
			mask2[index] = (double)tempMask2[index];
		}
	}

	double mean1 = 0.0, mean2 = 0.0, var1 = 0.0, var2 = 0.0, count = 0.0, coefficient = 0.0;

	for(SInt32 index = 0; index < boxSize; index++)
	{
		double maskFactor = mask1[index]*mask2[index];
		mean1 += maskFactor*image1[index];
		mean2 += maskFactor*image2[index];
		var1 += maskFactor*image1[index]*image1[index];
		var2 += maskFactor*image2[index]*image2[index];
		count += maskFactor*1.0;
		coefficient += maskFactor*image1[index]*image2[index];
	}

	// if more than one quarter of the correlation window is masked out,
	// consider the correlation to be too unreliable
	// (this is why we don't correlate the edges of the image)
	if(count < boxSize)
		return false;


	double radicand = (var1*count - mean1*mean1)*(var2*count - mean2*mean2);

	if(radicand <= 0.0) // the denominator is somehow invalid
	{
		printf("warning: bad denominator in correlation, possible bug in the code\n");
		return false;
	}

	coefficient = (coefficient*count - mean1*mean2)/sqrt(radicand);
	outFunctionValue = -coefficient; // we want to maximize the function, so minimize the negative function
	return true;
}

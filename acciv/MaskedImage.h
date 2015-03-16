#pragma once
#include "ScalarField2D.h"

class VectorField2D;

class MaskedImage :
	public ScalarField2D
{
public:
	MaskedImage(SInt32 inXSize, SInt32 inYSize, double inXLower, double inXUpper,
		double inYLower, double inYUpper)
		: ScalarField2D(inXSize, inYSize, inXLower, inXUpper, inYLower, inYUpper)
	{
		mask.setSize(getXSize()*getYSize());
		for(SInt32 index = 0; index < mask.getSize(); index++)
			mask[index] = true;
	}

	MaskedImage()
		: ScalarField2D()
	{
	}

	virtual ~MaskedImage(void)
	{
	}

	const bool& getMask(SInt32 index) const;
	bool& getMask(SInt32 index);
	const bool& getMask(SInt32 xIndex, SInt32 yIndex) const;
	bool& getMask(SInt32 xIndex, SInt32 yIndex);

	double getTime() const
	{
		return time;
	}

	void setTime(double inTime)
	{
		time = inTime;
	}

	void maskedInterpolate(const UArray<double> & xs, const UArray<double> & ys, UArray<double> & outValues, UArray<bool> & outMask) const;
	void maskedInterpolateAtOffsetPixels(double xLower, double yLower, SInt32 xCount, SInt32 yCount, UArray<double> & outValues, UArray<bool> & outMask) const;

	bool write(const UString & fileName, bool writeTime = false) const
	{
		if(!ScalarField2D::write(fileName))
			return false;
		if(!writeMaskHDF5File(fileName))
			return false;
		if(writeTime && !MaskedImage::writeTimeHDF5File(fileName, time))
			return false;
		return true;
	}

	bool read(const UString & fileName, bool readTime = false)
	{
		if(!ScalarField2D::read(fileName))
		{
			initialize(0,0,0.0,0.0,0.0,0.0);
			return false;
		}
		if(!readMaskHDF5File(fileName))
		{
			initialize(0,0,0.0,0.0,0.0,0.0);
			return false;
		}
		if(readTime && !MaskedImage::readTimeHDF5File(fileName,time))
		{
			initialize(0,0,0.0,0.0,0.0,0.0);
			return false;
		}
		return true;
	}

	bool advect(const VectorField2D & velocity, double deltaT, 
		double errorTolerance, SInt32 maximumTimeStepCount);
	bool advect(const VectorField2D & velocity, double deltaT,
		UArray<double> & inOutX, UArray<double> & inOutY, double errorTolerance, SInt32 maximumTimeStepCount);

	static bool readTimeHDF5File(const UString & fileName, double & outTime);
	static bool writeTimeHDF5File(const UString & fileName, double inTime);

private:
	bool readMaskHDF5File(const UString & fileName);
	bool writeMaskHDF5File(const UString & fileName) const;

	UArray<bool> mask;
	double time;
};


u_inline const bool& MaskedImage::getMask(SInt32 index) const
{
	return mask[index];
}

u_inline bool& MaskedImage::getMask(SInt32 index)
{
	return mask[index];
}

u_inline const bool& MaskedImage::getMask(SInt32 xIndex, SInt32 yIndex) const
{
	U_ASSERT((xIndex >= 0) & (xIndex < getXSize()) & (yIndex >= 0) & (yIndex < getYSize()));
	return mask[xIndex + getXSize()*yIndex];
}

u_inline bool& MaskedImage::getMask(SInt32 xIndex, SInt32 yIndex)
{
	U_ASSERT((xIndex >= 0) & (xIndex < getXSize()) & (yIndex >= 0) & (yIndex < getYSize()));
	return mask[xIndex + getXSize()*yIndex];
}


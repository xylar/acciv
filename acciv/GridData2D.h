#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>

template<class T>
class GridData2D
{
public:

	GridData2D(SInt32 inXSize, SInt32 inYSize, double inXLower, double inXUpper,
		double inYLower, double inYUpper)
	{
		initialize(inXSize,inYSize,inXLower,inXUpper,inYLower,inYUpper);
	}

	GridData2D()
	{
		initialize(0,0,0.0,0.0,0.0,0.0);
	}

	virtual ~GridData2D(void)
	{
	}

	// overloaded operator helpers
	const T& operator[](SInt32 index) const;
	T& operator[](SInt32 index);

	const T& operator()(SInt32 xIndex, SInt32 yIndex) const;
	T& operator()(SInt32 xIndex, SInt32 yIndex);

	virtual bool write(const UString & fileName) const = 0;

	virtual bool read(const UString & fileName) = 0;


	double getXLower() const
	{
		return xLower;
	}
	double getXUpper() const
	{
		return xUpper;
	}

	double getYLower() const
	{
		return yLower;
	}
	double getYUpper() const
	{
		return yUpper;
	}

	SInt32 getXSize() const
	{
		return xSize;
	}

	SInt32 getYSize() const
	{
		return ySize;
	}

	SInt32 getDataSize() const
	{
		return xSize*ySize;
	}

	double getDeltaX() const
	{
		return deltaX;
	}

	double getDeltaY() const
	{
		return deltaY;
	}

	const UArray<double>& getX() const
	{
		return x;
	}

	UArray<double>& getX()
	{
		return x;
	}

	const UArray<double>& getY() const
	{
		return y;
	}

	UArray<double>& getY()
	{
		return y;
	}

	const UArray<T>& getData() const
	{
		return data;
	}

	UArray<T>& getData()
	{
		return data;
	}

	// interpolate the values of data at the scattered data points (ioScatteredData.x and ioScatteredData.y)
	// using cubic splines
	// ioScatteredData.data will contain the interpolated values
	void interpolate(const UArray<double> & xs, const UArray<double> & ys, UArray<T> & outValues) const;

	void initialize(SInt32 inXSize, SInt32 inYSize, double inXLower, double inXUpper,
		double inYLower, double inYUpper)
	{
		xSize = inXSize;
		ySize = inYSize;
		xLower = inXLower;
		xUpper = inXUpper;
		yLower = inYLower;
		yUpper = inYUpper;
		x.setSize(inXSize);
		y.setSize(inYSize);
		data.setSize(xSize*ySize);
		setXandY();
	}

protected:
	void setXandY()
	{
		if((xSize == 0) || (ySize == 0))
			return;

		deltaX = (xUpper-xLower)/double(xSize-1);
		deltaY = (yUpper-yLower)/double(ySize-1);
		for(SInt32 xIndex = 0; xIndex < xSize; xIndex++)
			x[xIndex] = xIndex*deltaX + xLower;
		for(SInt32 yIndex = 0; yIndex < ySize; yIndex++)
			y[yIndex] = yIndex*deltaY + yLower;
	}

	double xLower, xUpper, yLower, yUpper;
	double deltaX, deltaY;
	
	SInt32 xSize, ySize;

	UArray<double> x;
	UArray<double> y;
	UArray<T> data;
};

template<class T>
u_inline const T& GridData2D<T>::operator[](SInt32 index) const
{
	return data[index];
}

template<class T>
u_inline T& GridData2D<T>::operator[](SInt32 index)
{
	return data[index];
}

template<class T>
u_inline const T& GridData2D<T>::operator()(SInt32 xIndex, SInt32 yIndex) const
{
	U_ASSERT((xIndex >= 0) & (xIndex < xSize) & (yIndex >= 0) & (yIndex < ySize));
	return data[xIndex + xSize*yIndex];
}

template<class T>
u_inline T& GridData2D<T>::operator()(SInt32 xIndex, SInt32 yIndex)
{
	U_ASSERT((xIndex >= 0) & (xIndex < xSize) & (yIndex >= 0) & (yIndex < ySize));
	return data[xIndex + xSize*yIndex];
}

template<class T>
void GridData2D<T>::interpolate(const UArray<double> & xs, const UArray<double> & ys, UArray<T> & outValues) const
{

	U_ASSERT((xs.getSize() == ys.getSize()) && (xs.getSize() == outValues.getSize()));

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
		double x = (xs[dataIndex]- xLower)/deltaX;
		double y = (ys[dataIndex]- yLower)/deltaY;
		SInt32 xIndex = (SInt32)x;
		SInt32 yIndex = (SInt32)y;

		// if the index is too close to the edge, extrapolate from the nearest location inbounds
		if(xIndex < 1)
			xIndex = 1;
		if(xIndex > xSize-3)
			xIndex = xSize-3;

		if(yIndex < 1)
			yIndex = 1;
		if(yIndex > ySize-3)
			yIndex = ySize-3;

		//if((xIndex-1 < 0) || (xIndex+2 > xSize-1)
		//	|| (yIndex-1 < 0) || (yIndex+2 > ySize-1))
		//{
		//	//To do: extrapolate/interpolate near the edges
		//	//some of the neighbors are out of range so we won't bother interpolating
		//	ioScatteredData.getData(dataIndex) = nan;
		//	continue;
		//}

		x = x - (double)xIndex;
		y = y - (double)yIndex;

		//U_ASSERT((x >= 0) & (x < 1) & (y >= 0) & (y < 1));

		xPowers[1] = x;
		xPowers[2] = x*x;
		xPowers[3] = x*x*x;
		yPowers[1] = y;
		yPowers[2] = y*y;
		yPowers[3] = y*y*y;

		T result = T(0.0);
		for(SInt32 k = 0; k < 4; k++)
		{
			T localResult = T(0.0);
			for(SInt32 j = 0; j < 4; j++)
				for(SInt32 i = 0; i < 4; i++)
					localResult = localResult + xPowers[i]*cubicSplineMatrix[j+4*i]
					              *data[(xIndex + j - 1) + xSize*(yIndex + k - 1)];
			for(SInt32 i = 0; i < 4; i++)
				result = result + yPowers[i]*cubicSplineMatrix[k+4*i]*localResult;
		}

		outValues[dataIndex] = result;

	}
}

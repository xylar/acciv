#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>

#include <algorithm>
#include "GridData2D.h"
#include <stdio.h>


template<class T>
class ScatteredData2D
{
public:
	ScatteredData2D(SInt32 inSize) : size(inSize), x(inSize), y(inSize), data(inSize), weights(inSize)
	{
		for(UInt32 index = 0; index < weights.getSize(); index++)
			weights[index] = 1.0;
	}

	ScatteredData2D() : size(0)
	{
	}

	virtual ~ScatteredData2D(void)
	{
	}

	SInt32 getSize() const
	{
		return size;
	}

	void setSize(SInt32 inSize)
	{
		size = inSize;
		x.setSize(inSize);
		y.setSize(inSize);
		data.setSize(inSize);
		bool initWeights = weights.getSize() < inSize;
		weights.setSize(inSize);
		if(initWeights)
			for(UInt32 index = 0; index < weights.getSize(); index++)
				weights[index] = 1.0;
	}

	const double& getX(SInt32 index) const;
	double& getX(SInt32 index);

	const double& getY(SInt32 index) const;
	double& getY(SInt32 index);

	const T& getData(SInt32 index) const;
	T& getData(SInt32 index);

	const double& getWeights(SInt32 index) const;
	double& getWeights(SInt32 index);

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

	const UArray<double>& getWeights() const
	{
		return weights;
	}

	UArray<double>& getWeights()
	{
		return weights;
	}

	virtual bool write(UString fileName) = 0; // must be implemented by the derived class
	virtual bool read(UString fileName) = 0; // must be implemented by the derived class

	void smoothFit(GridData2D<T> & field, UInt32 minControlPointScatteredNeighbors,
			UArray<T> & outResiduals, double & outMeanRes) const;

protected:

	// B-spline functions
	void getBSplineWeights(double x, double y, UInt32 nxLattice, UInt32 nyLattice,
			UInt32 & outXIndexLattice, UInt32 & outYIndexLattice, UArray<double> & outWeights) const;
	T evaluateBSpline(const UArray<T> & coeffs, const UArray<double> & weights, UInt32 nxLattice, UInt32 nyLattice,
			UInt32 xIndexLattice, UInt32 yIndexLattice) const;
	void evaluateGridBSpline(const UArray<T> & coeffs, UInt32 nxLattice, UInt32 nyLattice,
			UInt32 nx, UInt32 ny, double scale, UArray<T> & outF) const;
	void evaluateScatteredBSpline(const UArray<T> & coeffs,
			double xl, double xu, double yl, double yu,
			UInt32 nxLattice, UInt32 nyLattice,  UArray<T> & outF) const;

	void getPointBSplineCoeffs(const T & residual,
			double x, double y, UInt32 nxLattice, UInt32 nyLattice,
			UArray<T> & outNumerators, UArray<double> & outDenominators) const;

	void bSplineApproximation(const UArray<T> & residuals,
			double xl, double xu, double yl, double yu,
			UInt32 nxLattice, UInt32 nyLattice, UArray<T> & outControlPoints,
			UInt32 & outControlPointCount, UInt32 minimumPointsPerLatticeNode = 1) const;

	void subdivideControlPoints(
			UInt32 & inOutNxLattice, UInt32 & inOutNyLattice,
			const UArray<T> & inControlPoints, UArray<T> & outControlPoints) const;

	void multilevelBSplineApproximation(
			double xl, double xu, double yl, double yu,
			UInt32 maxLevelCount, UInt32 & inOutNxLattice, UInt32 & inOutNyLattice,
			UArray<T> & inOutResiduals, UArray<T> & outControlPoints,
			double & outScale,  double & outMeanResidual, UInt32 minimumPointsPerLatticeNode = 1) const;


	SInt32 size;
	UArray<double> x, y;
	UArray<T> data;
	UArray<double> weights;

};

template<class T>
u_inline const double& ScatteredData2D<T>::getX(SInt32 index) const
{
	return x[index];
}

template<class T>
u_inline double& ScatteredData2D<T>::getX(SInt32 index)
{
	return x[index];
}

template<class T>
u_inline const double& ScatteredData2D<T>::getY(SInt32 index) const
{
	return y[index];
}

template<class T>
u_inline double& ScatteredData2D<T>::getY(SInt32 index)
{
	return y[index];
}

template<class T>
u_inline const T& ScatteredData2D<T>::getData(SInt32 index) const
{
	return data[index];
}

template<class T>
u_inline T& ScatteredData2D<T>::getData(SInt32 index)
{
	return data[index];
}

template<class T>
u_inline const double& ScatteredData2D<T>::getWeights(SInt32 index) const
{
	return weights[index];
}

template<class T>
u_inline double& ScatteredData2D<T>::getWeights(SInt32 index)
{
	return weights[index];
}


template<class T>
void ScatteredData2D<T>::getBSplineWeights(double x, double y, UInt32 nxLattice, UInt32 nyLattice,
		UInt32 & outXIndex, UInt32 & outYIndex, UArray<double> & outWeights) const
{
	outXIndex = std::min(nxLattice-4,std::max((UInt32)0,(UInt32)(x)));
	outYIndex = std::min(nyLattice-4,std::max((UInt32)0,(UInt32)(y)));
	x -= outXIndex;
	y -= outYIndex;
	double bx[4], by[4];
	bx[0] = ((1-x)*(1-x)*(1-x));
	bx[1] = ((3*x-6)*x*x+4);
	bx[2] = (((-3*x+3)*x+3)*x+1);
	bx[3] = x*x*x;
	by[0] = ((1-y)*(1-y)*(1-y));
	by[1] = ((3*y-6)*y*y+4);
	by[2] = (((-3*y+3)*y+3)*y+1);
	by[3] = y*y*y;
	for(UInt32 yi = 0; yi < 4; yi++)
	{
		for(UInt32 xi = 0; xi < 4; xi++)
		{
			outWeights[xi + 4*yi] = bx[xi]*by[yi]/36.0;
		}
	}
}

template<class T>
u_inline T ScatteredData2D<T>::evaluateBSpline(const UArray<T> & controlPoints, const UArray<double> & weights, UInt32 nxLattice, UInt32 nyLattice,
		UInt32 xIndexLattice, UInt32 yIndexLattice) const
{
	T result = T(0.0);
	for(UInt32 yi = 0; yi < 4; yi++)
	{
		for(UInt32 xi = 0; xi < 4; xi++)
		{
			UInt32 cIndex = (xIndexLattice+xi) + nxLattice*(yIndexLattice+yi);
			result += weights[xi+4*yi]*controlPoints[cIndex];
		}
	}
	return result;
}

template<class T>
void ScatteredData2D<T>::evaluateGridBSpline(const UArray<T> & controlPoints, UInt32 nxLattice, UInt32 nyLattice,
		UInt32 nx, UInt32 ny, double scale, UArray<T> & outF) const
{
	for(UInt32 yIndex = 0; yIndex < ny; yIndex++)
	{
		for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
		{
			double x = scale*xIndex;
			double y = scale*yIndex;
			UInt32 outIndex = xIndex + nx*yIndex;
			UArray<double> weights(16);
			UInt32 xIndexLattice, yIndexLattice;
			getBSplineWeights(x, y, nxLattice, nyLattice, xIndexLattice, yIndexLattice, weights);
			outF[outIndex] = evaluateBSpline(controlPoints, weights, nxLattice, nyLattice,
					xIndexLattice, yIndexLattice);
		}
	}
}

template<class T>
void ScatteredData2D<T>::evaluateScatteredBSpline(const UArray<T> & controlPoints,
		double xl, double xu, double yl, double yu,
		UInt32 nxLattice, UInt32 nyLattice,  UArray<T> & outF) const
{
	double scaleX = (nxLattice-3)/(xu-xl);
	double scaleY = (nyLattice-3)/(yu-yl);
	UInt32 nPoints = getSize();
	outF.setSize(nPoints);
	for(UInt32 pIndex = 0; pIndex < nPoints; pIndex++)
	{
		double x = scaleX*(getX(pIndex)-xl);
		double y = scaleY*(getY(pIndex)-yl);
		UArray<double> weights(16);
		UInt32 xIndexLattice, yIndexLattice;
		getBSplineWeights(x, y, nxLattice, nyLattice, xIndexLattice, yIndexLattice, weights);
		outF[pIndex] = evaluateBSpline(controlPoints, weights, nxLattice, nyLattice,
				xIndexLattice, yIndexLattice);
	}
}


// The B-spline approximation (BA) algorithm from Lee, S., Wolberg, G., & Shin, S. (1997).
//   Scattered data interpolation with multilevel B-splines.
template<class T>
void ScatteredData2D<T>::bSplineApproximation(const UArray<T> & residuals,
		double xl, double xu, double yl, double yu,
		UInt32 nxLattice, UInt32 nyLattice, UArray<T> & outControlPoints,
		UInt32 & outControlPointCount, UInt32 minimumPointsPerLatticeNode) const
{

	outControlPointCount = 0;
	double scaleX = (nxLattice-3)/(xu-xl);
	double scaleY = (nyLattice-3)/(yu-yl);
	UInt32 nPoints = residuals.getSize();

	UArray<T> numerator(outControlPoints.getSize());
	UArray<double> denominator(outControlPoints.getSize());
	UArray<UInt32> count(outControlPoints.getSize());

	for(UInt32 index = 0; index < numerator.getSize(); index++)
	{
		numerator[index] = T(0.0);
		denominator[index] = 0.0;
		count[index] = 0;
	}

	for(UInt32 pIndex = 0; pIndex < nPoints; pIndex++)
	{
		double x = scaleX*(getX(pIndex)-xl);
		double y = scaleY*(getY(pIndex)-yl);
		UArray<double> bSplineWeights(16);
		UInt32 xIndexLattice, yIndexLattice;
		getBSplineWeights(x, y, nxLattice, nyLattice, xIndexLattice, yIndexLattice, bSplineWeights);
		double weightSquaredSum = 0.0;
		for(UInt32 wIndex = 0; wIndex < 16; wIndex++)
		{
			weightSquaredSum += bSplineWeights[wIndex]*bSplineWeights[wIndex];
		}
		double pointWeight = weights[pIndex];
		T residual = residuals[pIndex];
		for(UInt32 yi = 0; yi < 4; yi++)
		{
			for(UInt32 xi = 0; xi < 4; xi++)
			{
				double weight = bSplineWeights[xi+4*yi];
				double weightSquared = weight*weight;
				T phi = residual*weight/weightSquaredSum;
				UInt32 cIndex = (xIndexLattice+xi) + nxLattice*(yIndexLattice+yi);
				numerator[cIndex] += weightSquared*pointWeight*phi;
				denominator[cIndex] += weightSquared*pointWeight;
				count[cIndex]++;
			}
		}
	}

	UInt32 countSum = 0, maxCount = 0;
	for(UInt32 index = 0; index < outControlPoints.getSize(); index++)
	{
		countSum += count[index];
		maxCount = std::max(count[index],maxCount);
		if((count[index] >= minimumPointsPerLatticeNode) && (denominator[index] > 0.0))
		{
			outControlPoints[index] = numerator[index]/denominator[index];
			outControlPointCount++;
		}
		else
		{
			outControlPoints[index] = T(0.0);
		}
	}
}

template<class T>
void ScatteredData2D<T>::subdivideControlPoints(
		UInt32 & inOutNxLattice, UInt32 & inOutNyLattice,
		const UArray<T> & inControlPoints, UArray<T> & outControlPoints) const
{
	UInt32 inNx = inOutNxLattice;
	UInt32 inNy = inOutNyLattice;
	UInt32 outNx = 2*inNx-3;
	UInt32 outNy = 2*inNy-3;
	inOutNxLattice = outNx;
	inOutNyLattice = outNy;

	// subdivide the control points
	outControlPoints.setSize(0);
	outControlPoints.setSize(outNx*outNy);
	// control points 2i,2j in Lee et al. (1997) notation
	for(UInt32 yIndex = 1; yIndex < inNy-1; yIndex++)
	{
		for(UInt32 xIndex = 1; xIndex < inNx-1; xIndex++)
		{
			UInt32 xIndex2 = 2*xIndex-1;
			UInt32 yIndex2 = 2*yIndex-1;
			UInt32 coarseIndex = xIndex + inNx*yIndex;
			UInt32 fineIndex = xIndex2 + outNx*(yIndex2);
			outControlPoints[fineIndex] = (
				  inControlPoints[coarseIndex-1-inNx]
				+ inControlPoints[coarseIndex-1+inNx]
				+ inControlPoints[coarseIndex+1-inNx]
				+ inControlPoints[coarseIndex+1+inNx]
				+ 6.0*(inControlPoints[coarseIndex-1]
				     + inControlPoints[coarseIndex-inNx]
				     + inControlPoints[coarseIndex+1]
				     + inControlPoints[coarseIndex+inNx])
				+ 36.0*inControlPoints[coarseIndex])/64.0;
		}
	}
	// control points 2i,2j+1 in Lee et al. (1997) notation
	for(UInt32 yIndex = 0; yIndex < inNy-1; yIndex++)
	{
		for(UInt32 xIndex = 1; xIndex < inNx-1; xIndex++)
		{
			UInt32 xIndex2 = 2*xIndex-1;
			UInt32 yIndex2 = 2*yIndex;
			UInt32 coarseIndex = xIndex + inNx*yIndex;
			UInt32 fineIndex = xIndex2 + outNx*yIndex2;
			outControlPoints[fineIndex] = (
				  inControlPoints[coarseIndex-1]
			    + inControlPoints[coarseIndex-1+inNx]
			    + inControlPoints[coarseIndex+1]
			    + inControlPoints[coarseIndex+1+inNx]
			    + 6.0*(inControlPoints[coarseIndex]
			         + inControlPoints[coarseIndex+inNx]))/16.0;
		}
	}
	// control points 2i+1,2j in Lee et al. (1997) notation
	for(UInt32 yIndex = 1; yIndex < inNy-1; yIndex++)
	{
		for(UInt32 xIndex = 0; xIndex < inNx-1; xIndex++)
		{
			UInt32 xIndex2 = 2*xIndex;
			UInt32 yIndex2 = 2*yIndex-1;
			UInt32 coarseIndex = xIndex + inNx*yIndex;
			UInt32 fineIndex = xIndex2 + outNx*yIndex2;
			outControlPoints[fineIndex] = (
				  inControlPoints[coarseIndex-inNx]
			    + inControlPoints[coarseIndex+inNx]
			    + inControlPoints[coarseIndex+1-inNx]
			    + inControlPoints[coarseIndex+1+inNx]
			    + 6.0*(inControlPoints[coarseIndex]
			         + inControlPoints[coarseIndex+1]))/16.0;
		}
	}
	// control points 2i+1,2j+1 in Lee et al. (1997) notation
	for(UInt32 yIndex = 0; yIndex < inNy-1; yIndex++)
	{
		for(UInt32 xIndex = 0; xIndex < inNx-1; xIndex++)
		{
			UInt32 xIndex2 = 2*xIndex;
			UInt32 yIndex2 = 2*yIndex;
			UInt32 coarseIndex = xIndex + inNx*yIndex;
			UInt32 fineIndex = xIndex2 + outNx*yIndex2;
			outControlPoints[fineIndex] = (
				  inControlPoints[coarseIndex]
			    + inControlPoints[coarseIndex+inNx]
			    + inControlPoints[coarseIndex+1]
			    + inControlPoints[coarseIndex+1+inNx])/4.0;
		}
	}
}

// The multilevel B-spline approximation (BA) algorithm from Lee, S., Wolberg, G., & Shin, S. (1997).
//   Scattered data interpolation with multilevel B-splines.
template<class T>
void ScatteredData2D<T>::multilevelBSplineApproximation(
		double xl, double xu, double yl, double yu,
		UInt32 maxLevelCount, UInt32 & inOutNxLattice, UInt32 & inOutNyLattice,
		UArray<T> & inOutResiduals, UArray<T> & outControlPoints,
		double & inOutScale, double & outMeanResidual, UInt32 minimumPointsPerLatticeNode) const
{
	UInt32 nxLattice = inOutNxLattice;
	UInt32 nyLattice = inOutNyLattice;
	UInt32 nPoints = inOutResiduals.getSize();
	UArray<T> controlPoints(nxLattice*nyLattice);
	for(UInt32 index = 0; index < controlPoints.getSize(); index++)
		controlPoints[index] = T(0.0);
	outMeanResidual = 0.0;
	for(UInt32 pIndex = 0; pIndex < inOutResiduals.getSize(); pIndex++)
	{
		outMeanResidual += inOutResiduals[pIndex]*inOutResiduals[pIndex];
	}
	outMeanResidual = sqrt(outMeanResidual/inOutResiduals.getSize());
	printf("rms residual before fit %g\n", outMeanResidual);


	UInt32 lastLevelIndex = 0;
	UInt32 controlPointCount = 0;
	for(UInt32 levelIndex = 0; levelIndex < maxLevelCount; levelIndex++)
	{
		UArray<T> newControlPoints(nxLattice*nyLattice);
		bSplineApproximation(inOutResiduals, xl, xu, yl, yu,
				nxLattice, nyLattice, newControlPoints, controlPointCount,
				minimumPointsPerLatticeNode);

		printf("level %i: %i/%i nonzero control points found.\n", levelIndex,
				controlPointCount,newControlPoints.getSize());
		lastLevelIndex = levelIndex;

		if(controlPointCount == 0)
			break;


		UArray<T> scatteredResidual(nPoints);
		evaluateScatteredBSpline(newControlPoints, xl, xu, yl, yu,
				nxLattice, nyLattice,  scatteredResidual);

		outMeanResidual = 0.0;
		for(UInt32 pIndex = 0; pIndex < inOutResiduals.getSize(); pIndex++)
		{
			inOutResiduals[pIndex] -= scatteredResidual[pIndex];
			outMeanResidual += inOutResiduals[pIndex]*inOutResiduals[pIndex];
		}
		outMeanResidual = sqrt(outMeanResidual/inOutResiduals.getSize());
		printf("         rms residual after fit %g\n", outMeanResidual);


		for(UInt32 index = 0; index < controlPoints.getSize(); index++)
		{
			newControlPoints[index] += controlPoints[index];
		}

		subdivideControlPoints(nxLattice, nyLattice, newControlPoints, controlPoints);

		inOutScale *= 2;

	}
	printf("multilevelBSplineApproximation finished at level %i\n",lastLevelIndex);

	outControlPoints = controlPoints;
	inOutNxLattice = nxLattice;
	inOutNyLattice = nyLattice;

}

template<class T>
void ScatteredData2D<T>::smoothFit(GridData2D<T> & field,
		UInt32 minControlPointScatteredNeighbors,
		UArray<T> & outResiduals, double & outMeanRes) const
{

	// perform a linear fit to all data

	U_ASSERT(weights.getSize() == getSize);

	double A[3][3];
	T b[3];
	for(UInt32 j = 0; j<3; j++)
	{
		for(UInt32 i = 0; i<3; i++)
		{
			A[i][j] = 0.0;
		}
		b[j] = T(0.0);
	}
	for(UInt32 pIndex = 0; pIndex < getSize(); pIndex++)
	{
		double x = this->x[pIndex];
		double y = this->y[pIndex];
		double w = weights[pIndex];
		T z = getData(pIndex);
		A[0][0] += w*x*x;
		A[0][1] += w*x*y;
		A[0][2] += w*x;
		A[1][0] += w*y*x;
		A[1][1] += w*y*y;
		A[1][2] += w*y;
		A[2][0] += w*x;
		A[2][1] += w*y;
		A[2][2] += w*1;

		b[0] += w*x*z;
		b[1] += w*y*z;
		b[2] += w*z;
	}

	T linearCoefficients[3];
	double det =  A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])
                 -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                 +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);

	if (det == 0.0)
	{
		linearCoefficients[0] = T(0.0);
		linearCoefficients[1] = T(0.0);
		linearCoefficients[2] = T(0.0);
	}
	else
	{
		// Cramer's Rule (compute the determinant computed by replacing a given column of A by b
		// and divide by the determinant of A):
		linearCoefficients[0] = (b[0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])
								-A[0][1]*(b[1]*A[2][2]-A[1][2]*b[2])
								+A[0][2]*(b[1]*A[2][1]-A[1][1]*b[2]))/det;
		linearCoefficients[1] = (A[0][0]*(b[1]*A[2][2]-b[2]*A[1][2])
								-b[0]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
								+A[0][2]*(A[1][0]*b[2]-b[1]*A[2][0]))/det;
		linearCoefficients[2] = (A[0][0]*(A[1][1]*b[2]-A[2][1]*b[1])
                                -A[0][1]*(A[1][0]*b[2]-b[1]*A[2][0])
                                +b[0]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))/det;
	}

	// subtract the linear fit from the data to get the residuals for the first level
	outResiduals.setSize(0);
	outResiduals.setSize(getSize());

	for(UInt32 pIndex = 0; pIndex < getSize(); pIndex++)
	{
		outResiduals[pIndex] = getData(pIndex)
		    - linearCoefficients[0]*x[pIndex]
		    - linearCoefficients[1]*y[pIndex]
		    - linearCoefficients[2];
	}

	// determine the size of the coarsest lattice and the maximum number of levels

	// find the smallest power of 2 larger than each image dimension
	UInt32 xPower = 0, yPower = 0, nx = 1, ny = 1;
	while(nx < field.getXSize()-1)
	{
		nx *= 2;
		xPower++;
	}
	while(ny < field.getYSize()-1)
	{
		ny *= 2;
		yPower++;
	}
	UInt32 maxLevelCount = std::min(xPower,yPower);

	double xl = field.getXLower();
	double xu = nx*(field.getXUpper()-xl)/(field.getXSize()-1)+xl;
	double yl = field.getYLower();
	double yu = ny*(field.getYUpper()-yl)/(field.getYSize()-1)+yl;

	double scale = 1.0;


	for(UInt32 index = 0; index < maxLevelCount; index++)
	{
		nx /= 2;
		ny /= 2;
		scale /= 2.0;
	}

	nx += 3;
	ny += 3;

	// perform the MBA algorithm on the residual after the linear fit

	UArray<T> controlPoints;
	multilevelBSplineApproximation(
			xl, xu, yl, yu,
			maxLevelCount, nx, ny, outResiduals, controlPoints,
			scale, outMeanRes, minControlPointScatteredNeighbors);
	// interpolate using the resulting control points
	evaluateGridBSpline(controlPoints, nx, ny,
			field.getXSize(), field.getYSize(), scale, field.getData());

	// add the linear fit
	for(UInt32 yIndex = 0; yIndex < field.getYSize(); yIndex++)
	{
		for(UInt32 xIndex = 0; xIndex < field.getXSize(); xIndex++)
		{
			field(xIndex,yIndex) += linearCoefficients[0]*field.getX()[xIndex]
						  + linearCoefficients[1]*field.getY()[yIndex]
				          + linearCoefficients[2];
		}
	}
}

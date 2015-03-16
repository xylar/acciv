#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>

#include "ScatteredData2D.h"
#include "vector2.h"
//#include "QuadTree.h"

class VectorField2D;

class ScatteredVectorData2D : public ScatteredData2D<vector2d>
{
public:
	ScatteredVectorData2D(SInt32 inSize) : ScatteredData2D<vector2d>(inSize)
	{

	}

	ScatteredVectorData2D() : ScatteredData2D<vector2d>()
	{
	}

	virtual ~ScatteredVectorData2D(void)
	{
	}

	virtual bool write(UString fileName);
	virtual bool read(UString fileName);


/*
	void smoothFit(VectorField2D & field, UInt32 minBinPoints, double minVariance, double & outMeanRes,
		UArray<double> & outSmoothnessX, UArray<double> & outSmoothnessY, double testDataFraction = 0.0) const;
*/


private:

/*	bool smoothFitLevel(const QuadTreeLevel & level, UInt32 minBinPoints,
			double minVariance, UInt32 nx, UInt32 ny, UInt32 nPoints,
			const UArray<double> & x, const UArray<double> & y, double xl, double xu, double yl, double yu,
			const UArray<UInt32> testIndices, const UArray<UInt32> fitIndices,
			UArray<double> & fxRes, UArray<double> & fyRes, VectorField2D & field,
			UArray<double> & outSmoothnessX, UArray<double> & outSmoothnessY,
			double & outMeanRes, double & outTestRes) const;

	bool getBilinearCoeffs(const UArray<double> & x, const UArray<double> & y,
		const UArray<double> & fx, const UArray<double> & fy, double minVar, double * outCoeffs) const;
	void evaluateGridCubic(const UArray<double> & coeffs, UInt32 nxBin, UInt32 nyBin, UInt32 nx, UInt32 ny,
			UArray<double> & outFx, UArray<double> & outFy) const;
	void evaluateScatteredCubic(const UArray<double> & coeffs, const UArray<double> & x,
			const UArray<double> & y, double xl, double xu, double yl, double yu,
			UInt32 nxBin, UInt32 nyBin,  UArray<double> & outFx, UArray<double> & outFy) const;
	void evaluateCubic(double x, double y, UInt32 nxBin, UInt32 nyBin,
			const UArray<double> & coeffs, double & outFx, double & outFy) const;*/


};


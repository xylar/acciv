#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>

class VectorField2D;
class ScatteredVectorData2D;

class TiePointSet
{
public:
	TiePointSet(void)
	{
	}

	TiePointSet(const UString & fileName)
	{
		read(fileName);
	}
	~TiePointSet(void)
	{
	}

	UArray<double> & getX1()
	{
		return x1;
	}
	UArray<double> & getY1()
	{
		return y1;
	}
	UArray<double> & getX2()
	{
		return x2;
	}
	UArray<double> & getY2()
	{
		return y2;
	}
	UArray<double> & getCorrelationCoefficients()
	{
		return correlationCoefficients;
	}

	UArray<double> & getDeltaTs()
	{
		return deltaTs;
	}
	UArray<SInt32> & getLowerIndexDeltaT()
	{
		return lowerIndexDeltaT;
	}
	UArray<SInt32> & getUpperIndexDeltaT()
	{
		return upperIndexDeltaT;
	}


	bool read(const UString & fileName);
	bool write(const UString & fileName) const;

	bool unadvect(VectorField2D & velocity1, VectorField2D & velocity2, 
		const UArray<double> & deltaT1s, const UArray<double> &  deltaT2s,
		double errorTolerance, SInt32 maximumTimeStepCount);

	bool computeCurvedPathVelocities(const VectorField2D & velocity, SInt32 pathVectorCount,
		double errorTolerance, SInt32 maximumTimeStepCount,
		ScatteredVectorData2D & outVelocities) const;

private:
	UArray<double> x1, y1, x2, y2, correlationCoefficients;
	UArray<double> deltaTs;
	UArray<SInt32> lowerIndexDeltaT, upperIndexDeltaT;
};

#pragma once

#include "ScatteredData2D.h"

class ScatteredScalarData2D : public ScatteredData2D<double>
{
public:
	ScatteredScalarData2D(SInt32 inSize) : ScatteredData2D<double>(inSize)
	{

	}

	ScatteredScalarData2D() : ScatteredData2D<double>()
	{
	}

	virtual ~ScatteredScalarData2D(void)
	{
	}

	virtual bool write(UString fileName);
	virtual bool read(UString fileName);

};

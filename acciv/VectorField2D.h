#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>

#include "GridData2D.h"
#include "vector2.h"

class VectorField2D : public GridData2D<vector2d>
{
public:
	VectorField2D(SInt32 inXSize, SInt32 inYSize, double inXLower, double inXUpper,
		double inYLower, double inYUpper)
		: GridData2D<vector2d>(inXSize,inYSize,inXLower,inXUpper,inYLower,inYUpper)
	{
	}

	VectorField2D()
	{
	}

	virtual ~VectorField2D(void)
	{
	}

	virtual bool write(const UString & fileName) const;

	virtual bool read(const UString & fileName);

};

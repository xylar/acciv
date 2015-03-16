#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>

#include "GridData2D.h"

class ScalarField2D : public GridData2D<double>
{
public:

	ScalarField2D(SInt32 inXSize, SInt32 inYSize, double inXLower, double inXUpper,
		double inYLower, double inYUpper)
	    : GridData2D<double>(inXSize,inYSize,inXLower,inXUpper,inYLower,inYUpper)
	{
	}

	ScalarField2D()
		: GridData2D<double>()
	{
	}

	virtual ~ScalarField2D(void)
	{
	}


	virtual bool write(const UString & fileName) const;

	virtual bool read(const UString & fileName);

};

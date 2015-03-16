#pragma once

#include "PowelMinimizer.h"

template<class T>
class SquareBoundedMinimizer : public PowelMinimizer<T>
{
public:
	SquareBoundedMinimizer(T & inFunction, double inXLower, double inXUpper, double inYLower, double inYUpper,
		double inTolerance = 3.0e-8)
		:PowelMinimizer<T>(inFunction, inTolerance), xLower(inXLower), xUpper(inXUpper),
		yLower(inYLower), yUpper(inYUpper)
	{
	}

	virtual void get1DBounds(const UArray<double> & p, const UArray<double> & xi,  
		double & a, double & b, double & c)
	{
		const double kTiny = 1e-16;
		U_ASSERT(p.getSize() == 2);

		if(xi[0] == 0.0)
		{
			if(xi[1] > 0.0)
			{
				a = (yLower - p[1])/xi[1];
				c = (yUpper - p[1])/xi[1];
			}
			else
			{
				a = (yUpper - p[1])/xi[1];
				c = (yLower - p[1])/xi[1];
			}
		}
		else if(xi[1] == 0.0)
		{
			if(xi[0] > 0.0)
			{
				a = (xLower - p[0])/xi[0];
				c = (xUpper - p[0])/xi[0];
			}
			else
			{
				a = (xUpper - p[0])/xi[0];
				c = (xLower - p[0])/xi[0];
			}
		}
		else
		{
			UArray<double> as(4);
			as[0] = (xLower - p[0])/xi[0];
			as[1] = (xUpper - p[0])/xi[0];
			as[2] = (yLower - p[1])/xi[1];
			as[3] = (yUpper - p[1])/xi[1];
			as.quickSort();
			a = as[1];
			c = as[2];
			
		}
		b = 0;
	}
private:

	double xLower, xUpper, yLower, yUpper;
};


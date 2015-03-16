#pragma once

#include "BrentMinimizer.h"
#include <containers/UArray.h>

template<class T>
class PowelFunctor
{
public:
	PowelFunctor(const UArray<double> & inP, const UArray<double> & inXi, T & inFunction)
		:p(inP), xi(inXi), function(inFunction)
	{
		arraySize = p.getSize();
		U_ASSERT(xi.getSize() == arraySize);
		xt.setSize(arraySize);
	}

	bool operator() (const double x, double & outFunctionValue)
	{
		for(SInt32 index = 0; index < arraySize; index++)
			xt[index] = p[index] + x*xi[index];
		 return function(xt, outFunctionValue);
	}

private:
	const UArray<double> & p, & xi;
	SInt32 arraySize;
	T & function;
	UArray<double> xt;
};

template<class T>
class PowelMinimizer
{
public:
	PowelMinimizer(T & inFunction, double inTolerance = 3.0e-8)
		: function(inFunction), tolerance(inTolerance)
	{
	}

	virtual ~PowelMinimizer()
	{
	}

	bool lineMinimize(double & functionMinimum)
	{
		SInt32 arraySize = p.getSize();
		PowelFunctor<T> functor(p,xi,function);
		BrentMinimizer<PowelFunctor<T> > brentMinimizer(functor, tolerance);
		// bracket the function
		double a = 0.0, b = 1.0;
		bool result = brentMinimizer.bracket(a, b);
		if(!result)
			return false;
		//double a,b,c;
		//get1DBounds(p,xi,a,b,c);
		//bool result = brentMinimizer.setBoundingPoints(a,b,c);
		double fMin, xMin;
		result = brentMinimizer.minimize(xMin, fMin);
		if(!result)
			return false;
		for(SInt32 index = 0; index < arraySize; index++)
		{
			xi[index] *= xMin;
			p[index] += xi[index];
		}
		functionMinimum = fMin;
		return true;
	}

	bool minimize(UArray<double> & ioPoint, double & outFunctionValue)
	{
		const SInt32 kMaximumIterationCount = 200;
		SInt32 arraySize = ioPoint.getSize();
		UArray<double> xiMatrix(arraySize*arraySize);
		for(SInt32 index = 0; index < arraySize; index++)
			xiMatrix[index + arraySize*index] = 1.0; // make xiMatrix an identity matrix;

		p = ioPoint;
		UArray<double> pt(arraySize), ptt(arraySize);

		xi.setSize(arraySize);
		double fN;
		bool result = function(p, fN);
		if(!result)
			return false;
		pt = p;

		double fE;
		for(SInt32 iterationIndex = 0; iterationIndex < kMaximumIterationCount; iterationIndex++)
		{
			double f0 = fN;
			SInt32 bigDirection = 0;
			double deltaF = 0.0;
			for(SInt32 outerIndex = 0; outerIndex < arraySize; outerIndex++)
			{
				for(SInt32 innerIndex = 0; innerIndex < arraySize; innerIndex++)
					xi[innerIndex] = xiMatrix[innerIndex + arraySize*outerIndex];
				fE = fN;
				if(!lineMinimize(fN))
					return false;
				if(fE - fN > deltaF)
				{
					deltaF = fE - fN;
					bigDirection = outerIndex + 1;
				}
			}
			if((f0 - fN) <= tolerance)
			{
				ioPoint = p;
				outFunctionValue = fN;
				return true;
			}
			for(SInt32 index = 0; index < arraySize; index++)
			{
				ptt[index] = 2.0*p[index] - pt[index];
				xi[index] = p[index] - pt[index];
				pt[index] = p[index];
			}
			result = function(ptt,fE);
			if(!result)
				return false;
			if(fE < f0)
			{
				double t = 2.0*(f0 - 2.0*fN + fE)*(f0 - fN - deltaF)*(f0 - fN - deltaF)
					- deltaF*(f0 - fE)*(f0 - fE);
				if(t < 0.0)
				{
					if(!lineMinimize(fN))
						return false;
					double xiMax = 0.0;
					for(SInt32 index = 0; index < arraySize; index++)
					{
						xiMax = std::max(xiMax, xi[index]);
					}
					if(xiMax != 0.0)
					{
						for(SInt32 index = 0; index < arraySize; index++)
						{
							xiMatrix[index + arraySize*(bigDirection-1)] = xiMatrix[index + arraySize*(arraySize-1)];
							xiMatrix[index + arraySize*(arraySize-1)] = xi[index];
						}
					}
				}
			}


		}
		return false;


	}

	//virtual void get1DBounds(const UArray<double> & p, const UArray<double> & xi,  
	//	double & a, double & b, double & c) = 0;

private:
	UArray<double> p, xi;
	T & function;
	double tolerance;
};

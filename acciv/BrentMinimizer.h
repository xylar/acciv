#pragma once

#include <core/UTypes.h>

template<class T>
class BrentMinimizer
{
public:
	BrentMinimizer(T & inFunction, double inTolerance = 3.0e-8)
		: boundA(0), boundB(0), boundC(0), boundFa(0), boundFb(0), boundFc(0), function(inFunction), initialized(false),
		tolerance(inTolerance)
	{
	}

	virtual ~BrentMinimizer()
	{
	}

	u_inline void setBounds(double inA, double inB, double inC,
		double inFa, double inFb, double inFc)
	{
		if(inA > inC)
		{
			swap(inA,inC);
			swap(inFa, inFc);
		}
		U_ASSERT((inA <= inB) && (inB <= inC) && (inFb <= inFa) && (inFb <= inFc));

		boundA = inA;
		boundB = inB;
		boundC = inC;
		boundFa = inFa;
		boundFb = inFb;
		boundFc = inFc;
		initialized = true;
	}

	u_inline bool setBoundingPoints(double inA, double inB, double inC)
	{
		U_ASSERT((inA < inB) && (inB < inC));
		boundA = inA;
		boundB = inB;
		boundC = inC;
		bool result = function(boundA, boundFa);
		if(!result)
			return false;
		result = function(boundB, boundFb);
		if(!result)
			return false;
		result = function(boundC, boundFc);
		if(!result)
			return false;

		U_ASSERT((boundFb < boundFa) && (boundFb < boundFc));
		initialized = true;
	}

	u_inline bool bracket(double inA, double inB)
	{
		const double kGoldenRatio = 0.3819660;
		const double kMagnificationLimit = 100.0;
		const double kTiny = 1e-20; // used to prevent division by zero;

		double ax = inA;
		double bx = inB;
		double cx, fa, fb, fc, fu;
		bool result;
		result = function(ax, fa);
		if(!result)
			return false;
		result = function(bx, fb);
		if(!result)
			return false;

		if(fb > fa)
		{
			// swap the points so that the direction from a to b is downhill
			swap(ax,bx);
			swap(fa,fb);
		}

		cx = bx + kGoldenRatio*(bx - ax); // first guess for c
		result = function(cx, fc);
		if(!result)
			return false;

		while(fb > fc)
		{
			double r = (bx-ax)*(fb-fc); // compute u by parabolic extrapolation from a, b and c
			double q = (bx-cx)*(fb-fa);
			double denom;
			if(abs(q-r) < kTiny) 
			{
				if(q > r)
					denom = 2.0*kTiny;
				else
					denom = -2.0*kTiny;
			}
			else
			{
				denom = 2.0*(q-r);
			}

			double u = bx - ((bx-cx)*q - (bx-ax)*r)/denom;

			double uLimit = bx + kMagnificationLimit*(cx-bx);

			if((bx-u)*(u-cx) > 0.0) // parabolic u is between b and c
			{
				result = function(u, fu);
				if(!result)
					return false;
				if(fu < fc) // got a minimum between b and c
				{
					setBounds(bx,u,cx,fb,fu,fc);
					return true;
				}
				else if(fu > fb) // got a minimum between a and u
				{
					setBounds(ax,bx,u,fa,fb,fu);
					return true;
				}
				u = cx + kGoldenRatio*(cx - bx); // parabolic fit was not used.  Use default magnification.
				result = function(u, fu);
				if(!result)
					return false;
			}
			else if((cx-u)*(u-uLimit) > 0.0) // parabolic fit between c and its allowed limit
			{
				result = function(u, fu);
				if(!result)
					return false;
				if(fu < fc)
				{
					shift3(bx, cx, u, u+kGoldenRatio*(u-cx));
					shift2(fb, fc, fu);
					result = function(u, fu);
					if(!result)
						return false;
				}
			}
			else if((u-uLimit)*(uLimit-cx) >= 0.0) // limit parabloic u to maximum allowed value
			{
				u = uLimit;
				result = function(u, fu);
				if(!result)
					return false;
			}
			else // reject parabolic u, use default magnification
			{
				u = cx + kGoldenRatio*(cx-bx);
				result = function(u, fu);
				if(!result)
					return false;
			}
			// eliminate oldest point and continue
			shift3(ax,bx,cx,u);
			shift3(fa,fb,fc,fu);
		}
		setBounds(ax,bx,cx,fa,fb,fc);
		return true;
	}


	bool minimize(double & outXMin, double & outFMin)
	{
		U_ASSERT(initialized);
		const SInt32 kMaximumIterationCount = 100;
		const double kGoldenRatio = 0.3819660;

		double u;
		double d = 0.0;
		double e = 0.0;

		double a = boundA;
		double b = boundC;
		double x = boundB;
		double w = boundB;
		double v = boundB;
		double fx;
		bool result = function(x, fx);
		if(!result)
			return false;
		double fw = fx;
		double fv = fx;
		for(SInt32 iterationIndex = 0; iterationIndex < kMaximumIterationCount; iterationIndex++)
		{
			double xMidPoint = 0.5*(a+b);
			if((abs(x - xMidPoint) < tolerance) && (b - a) < 4*tolerance)
			{
				outFMin = fx;
				outXMin = x;
				return true;
			}
			bool didParabolicStep = false;
			if(abs(e) > tolerance) // try the parabolic fit
			{
				double r = (x-w)*(fx-fv);
				double q = (x-v)*(fx-fw);
				double p = (x-v)*q - (x-w)*r;
				q = 2.0*(q-r);
				if(q > 0.0)
					p = -p;
				else
					q = -q;
				double eTemp = e;
				e = d;
				if((abs(p) < abs(0.5*q*eTemp)) && (p > q*(a-x))
					&& (p < q*(b-x)))
				{
					d = p/q;
					u = x + d;
					if((u-a < 2.0*tolerance) || (b-u < 2.0*tolerance))
					{
						if(x < xMidPoint)
							d = tolerance;
						else
							d = -tolerance;
					}
					didParabolicStep = true;
				}
			}

			if(!didParabolicStep) // do the golden section step
			{
				if(x >= xMidPoint)
					e = a - x; 
				else
					e = b - x;
				d = kGoldenRatio*e;
			}

			if(abs(d) >= tolerance)
				u = x+d;
			else if(d > 0)
				u = x + tolerance;
			else
				u = x - tolerance;
			double fu;
			result = function(u, fu);
			if(!result)
				return false;

			if(fu <= fx)
			{
				if(u >= x)
					a = x;
				else
					b = x;
				shift3(v,w,x,u);
				shift3(fv,fw,fx,fu);
			}
			else
			{
				if(u < x)
					a = u;
				else
					b = u;
				if((fu <= fw) || (w == x))
				{
					shift2(v,w,u);
					shift2(fv,fw,fu);
				}
				else if((fu < fv) || (v == x) || (v == w))
				{
					v = u;
					fv = fu;
				}
			}
		}
		// too many iterations
		return false;
	}


private:

	u_inline void swap(double & a, double & b)
	{
		double temp = a;
		a = b;
		b = temp;
	}

	u_inline void shift3(double & a, double & b, double & c, double d)
	{
		a = b;
		b = c;
		c = d;
	}

	u_inline void shift2(double & a, double & b, double c)
	{
		a = b;
		b = c;
	}

	double boundA,boundB,boundC, boundFa,boundFb,boundFc;

	bool initialized;

	T & function;

	double tolerance;



};


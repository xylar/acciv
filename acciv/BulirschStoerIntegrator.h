#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>
#include<algorithm>

using namespace std;


#if	defined(_MSC_VER)
#pragma	warning( disable: 4786 ) //	identifier was truncated to	'255' characters in	the	debug information

#pragma	inline_depth( 255 )
#pragma	inline_recursion( on )
#pragma	auto_inline( on	)

#define	inline __forceinline 
#endif
const double kSafetyFactor1 = 0.25;
const double kSafetyFactor2 = 0.7;
const double kMaximumStepReduction = 1e-5;
const double kMinimumStepReduction = 0.7;
const double kMaximumScale = 0.1;

template<class ChildClass, int kMaxSubdivisionPairs=8>
class BulirschStoerIntegrator
{
public:

	inline BulirschStoerIntegrator()
		:initialized(false)
	{
	}

	inline ~BulirschStoerIntegrator(void)
	{
	}

	inline void initialize(SInt32 inVectorSize, double inErrorTolerance = 1e-6)
	{
		vectorSize = inVectorSize;

		polynomialTableau.setSize(vectorSize*kMaxSubdivisionPairs);

		errorTolerance = inErrorTolerance;
		xNew = -1e29;
		double epsilon1 = kSafetyFactor1*errorTolerance;

		workCost[0] = 3.0;
		for(SInt32 k = 1; k<=kMaxSubdivisionPairs; k++)
		{
			workCost[k] = workCost[k-1] + 2.0*(k+1);
		}

		for(SInt32 q = 1; q<kMaxSubdivisionPairs; q++)
		{
			for(SInt32 k = 0; k<q; k++)
			{
				double power = (workCost[k+1] - workCost[q+1])/((2.0*(k+1)+1.0)*(workCost[q+1]-workCost[0]+1.0));
				alphas[k][q] = pow(epsilon1, power);
			}
			alphas[q][q] = 1.0;
		}

		optimalIndex = 1;
		while((optimalIndex < kMaxSubdivisionPairs-1) 
			&& (workCost[optimalIndex+1] <= workCost[optimalIndex]*alphas[optimalIndex-1][optimalIndex]))
		{
			optimalIndex++;
		}
		maximumIndex = optimalIndex;
		initialized = true;
	}

	inline bool takeStep(UArray<double> & ys, UArray<double> & dydxs, double & x, double inTrialStepSize,
		const UArray<double> & yScalings, double & outActualStepSize, double & outNextStepSize)
	{
		if(!initialized)
			initialize(ys.getSize());

		double stepSize = inTrialStepSize;
//printf("initial step size: %g\n", stepSize);

		UArray<double> ySaves(vectorSize);

		for(SInt32 index=0; index<vectorSize; index++)
		{
			ySaves[index] = ys[index];
		}

		bool first = (stepSize != outNextStepSize) || (x != xNew); // sloppy
		if(first)
		{
			optimalIndex = maximumIndex;
		}

		bool reduced = false;

		bool converged = false;

		double errors[kMaxSubdivisionPairs];

		SInt32 previousIndex = 0;

		while(!converged)
		{

			double reduction;

			for(SInt32 index=0; index<=maximumIndex; index++)
			{
				xNew = x+stepSize;
//printf("%i %g %g %g\n", index, xNew, x, stepSize);
				if(xNew == x)
				{
					fprintf(stderr, "Step size underflow in BulirschStoerIntegrater::takeStep().\n");
					return false;
				}
				UArray<double> yErrors(vectorSize);
				UArray<double> yTemps(vectorSize);
				takeModifiedMidpointStep(ySaves, dydxs, x, stepSize, 2*(index+1), yTemps);
				double xEstimate = (stepSize*stepSize)/(4.0*(index+1)*(index+1));
				polynomialZeroExptrapolate(index, xEstimate, yTemps, ys, yErrors);
//printf("%i %g %g %g %g\n", index, xNew, yTemps[0], ys[0], yErrors[0]);

				if(index != 0)
				{
					previousIndex = index-1;
					double maximumError = fabs(yErrors[0]/yScalings[0]);
					SInt32 maxErrorIndex = 0;
					for(SInt32 i=1; i<vectorSize; i++)
					{
						double newError = fabs(yErrors[i]/yScalings[i]);
						if(newError > maximumError)
						{
							maximumError = newError;
							maxErrorIndex = i;
//printf("error: %g %g\n", yErrors[i], yScalings[i]);
						}
					}
//printf("max error: %i, %g\n", maxErrorIndex, maximumError);
					maximumError = maximumError/errorTolerance;
					errors[previousIndex] = pow(maximumError/kSafetyFactor1, 1.0/(2*index + 1));


					if((index >= optimalIndex) || first)
					{
						if(maximumError < 1.0)
						{
//printf("converged, error: %g\n", maximumError);
							converged = true;
							break;
						}
						if((index == maximumIndex) || (index == optimalIndex+1))
						{
							reduction = kSafetyFactor2/errors[previousIndex];
//printf("reduction type 1: %g\n", reduction);
							break;
						}
						else if(index == optimalIndex)
						{
							if(alphas[optimalIndex-1][optimalIndex] < errors[previousIndex])
							{
								reduction = 1.0/errors[previousIndex];
//printf("reduction type 2: %g\n", reduction);
								break;
							}
						}
						else if(optimalIndex == maximumIndex)
						{
							if(alphas[previousIndex][maximumIndex-1] < errors[previousIndex])
							{
								reduction = alphas[previousIndex][maximumIndex-1]*kSafetyFactor2/errors[previousIndex];
//printf("reduction type 3: %g\n", reduction);
								break;
							}
						}
						else if(alphas[previousIndex][optimalIndex] < errors[previousIndex])
						{
							reduction = alphas[previousIndex][optimalIndex-1]/errors[previousIndex];
//printf("reduction type 4: %g\n", reduction);
							break;
						}
					}
				}
			}

			if(!converged)
			{
				reduction = std::min(reduction, kMinimumStepReduction);
				reduction = std::max(reduction, kMaximumStepReduction);
				stepSize *= reduction;
				reduced = true;
//printf("step size reduced to: %g %g\n", stepSize, reduction);
			}
		}

		x = xNew;
//printf("x: %g\n", x);

		outActualStepSize = stepSize;
//printf("outActualStepSize: %g\n", outActualStepSize);

		first = false;

		SInt32 finalIndex = previousIndex + 1;
//printf("finalIndex: %i\n", finalIndex);

		double scale = std::max(errors[0], kMaximumScale);
//printf("scale: %g, %g, %g\n", scale, errors[0], kMaximumScale);
		double minimumWork = scale*workCost[1];
//printf("minimumWork: %g, %g\n", minimumWork, workCost[1]);
		optimalIndex = 1;
		for(SInt32 index = 1; index <= previousIndex; index++)
		{
//printf("index: %i\n", index);
			double factor = std::max(errors[index], kMaximumScale);
//printf("factor: %g, %g, %g\n", factor, errors[index], kMaximumScale);
			double work = factor*workCost[index+1];
//printf("work: %g, %g\n", work, workCost[index+1]);
			if(work < minimumWork)
			{
				scale = factor;
				minimumWork = work;
				optimalIndex = index + 1;
//printf("new minimumWork: %g, %g %i\n", minimumWork, scale, optimalIndex);
			}
		}

		outNextStepSize = stepSize/scale;
//printf("outNextStepSize %g %g %g\n", outNextStepSize, stepSize, scale);

		if((optimalIndex >= finalIndex) && (optimalIndex != maximumIndex) && (!reduced))
		{
			double factor = std::max(scale/alphas[optimalIndex-1][optimalIndex], kMaximumScale);
			if(workCost[optimalIndex + 1]*factor <= minimumWork)
			{
				outNextStepSize = stepSize/factor;
				optimalIndex++;
			}
		}

//printf("next step size: %g\n", outNextStepSize);
		return true;
	}

	inline void takeModifiedMidpointStep(const UArray<double> & ys, const UArray<double> & dydxs, double inX, double inTotalStepSize,
		SInt32 stepCount, UArray<double> & outYs)
	{
		ChildClass * child = static_cast<ChildClass*>(this);
		double stepSize = inTotalStepSize/stepCount;
		UArray<double> previousZs(vectorSize);
		UArray<double> currentZs(vectorSize);
		UArray<double> tempDydxs(vectorSize);
		for(SInt32 index = 0; index<vectorSize; index++)
		{
			previousZs[index] = ys[index];
			currentZs[index] = ys[index] + stepSize*dydxs[index];
		}
		double x = inX + stepSize;

		child->computeDerivatives(x, currentZs, tempDydxs);

		for(SInt32 n=1; n<stepCount; n++)
		{
			for(SInt32 index=0; index<vectorSize; index++)
			{
				double swap = previousZs[index] + 2.0*stepSize*tempDydxs[index];
				previousZs[index] = currentZs[index];
				currentZs[index] = swap;
			}
			x += stepSize;
			child->computeDerivatives(x, currentZs, tempDydxs);
		}

		for(SInt32 index = 0; index<vectorSize; index++)
		{
			outYs[index] = 0.5*(previousZs[index]+currentZs[index]+stepSize*tempDydxs[index]);
		}
 
	}

	inline void polynomialZeroExptrapolate(SInt32 callSequenceIndex, double currentX, UArray<double> & currentYs,
		UArray<double> & outYs, UArray<double> & yErrors)
	{
		polynomialXs[callSequenceIndex] = currentX;
		for(SInt32 index=0; index < vectorSize; index++)
		{
			yErrors[index] = currentYs[index];
			outYs[index] = currentYs[index];
		}
		if(callSequenceIndex == 0)
		{
			for(SInt32 index=0; index < vectorSize; index++)
			{
				polynomialTableau[index*kMaxSubdivisionPairs] = currentYs[index];
			}
		}
		else
		{
			UArray<double>  d(vectorSize);
			for(SInt32 index=0; index < vectorSize; index++)
			{
				d[index] = currentYs[index];
			}
			for(SInt32 outerIndex=0; outerIndex<callSequenceIndex; outerIndex++)
			{
				double delta = 1.0/(polynomialXs[callSequenceIndex-outerIndex-1]-currentX);
				double f1 = currentX*delta;
				double f2 = polynomialXs[callSequenceIndex-outerIndex-1]*delta;
				for(SInt32 innerIndex=0; innerIndex<vectorSize; innerIndex++)
				{
					double tableauValue = polynomialTableau[innerIndex*kMaxSubdivisionPairs + outerIndex];
					polynomialTableau[innerIndex*kMaxSubdivisionPairs + outerIndex] = yErrors[innerIndex];
					delta = d[innerIndex] - tableauValue;
					yErrors[innerIndex] = f1*delta;
					d[innerIndex] = f2*delta;
					outYs[innerIndex] += yErrors[innerIndex];
				}
			}
			for(SInt32 index=0; index < vectorSize; index++)
			{
				polynomialTableau[index*kMaxSubdivisionPairs + callSequenceIndex] = yErrors[index];
//printf("error[%i]: %g\n", index, yErrors[index]);
			}
		}
	}

private:
	bool initialized;
	double errorTolerance;
	double alphas[kMaxSubdivisionPairs][kMaxSubdivisionPairs];
	double workCost[kMaxSubdivisionPairs+1];

	SInt32 optimalIndex, maximumIndex;
	double polynomialXs[kMaxSubdivisionPairs];
	UArray<double> polynomialTableau;

	double xNew, nextTimeStep;
	SInt32 vectorSize;
};

#pragma once

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <core/UTypes.h>
#include <containers/UArray.h>
#include "BulirschStoerIntegrator.h"
#include "VectorField2D.h"

#if	defined(_MSC_VER)
#pragma	warning( disable: 4786 ) //	identifier was truncated to	'255' characters in	the	debug information

#pragma	inline_depth( 255 )
#pragma	inline_recursion( on )
#pragma	auto_inline( on	)

#define	inline __forceinline 
#endif

class ParticleIntegrator : public BulirschStoerIntegrator<ParticleIntegrator,8>
{
public:

	inline ParticleIntegrator(const VectorField2D & inVelocity)
		:velocity(inVelocity)
	{
	}

	inline ~ParticleIntegrator(void)
	{
	}

	inline bool integrate(UArray<double> & xs, double t1, double t2, 
		double errorTolerance, double & nextTimeStep, double minimumTimeStep, SInt32 maximumTimeStepCount)
	{
		const double kTiny = 1e-30;

		initialize(xs.getSize(), errorTolerance);

		minimumTimeStep = fabs(minimumTimeStep);

		double t = t1;
		double timeStep = (t2 > t1) ? fabs(nextTimeStep) : -fabs(nextTimeStep);
		SInt32 okStepCount = 0;
		SInt32 badStepCount = 0;

		UArray<double> dxdts(xs.getSize());
		UArray<double> xScalings(xs.getSize());
		
		double actualTimeStep = -1e29;
		nextTimeStep = -1e29;

		for(SInt32 stepIndex = 0; stepIndex < maximumTimeStepCount; stepIndex++)
		{
			computeDerivatives(t, xs, dxdts);
			for(SInt32 index=0; index<xs.getSize(); index++)
			{
				xScalings[index] = fabs(xs[index]) + fabs(timeStep*dxdts[index]) + kTiny;
//printf("scaling %i: %g, %g, %g\n", index, timeStep, xs[index], dxdts[index]);
//double temp1 = fabs(xs[index]);
//double temp2 = timeStep*dxdts[index];
//printf("%g, %g, %g\n", xScalings[index], temp1, temp2);
			}

			if((t + timeStep - t2)*(t + timeStep - t1) > 0) // if we'll go past t2 then decrease the timestep
                        {
				timeStep = t2 - t;
				if(fabs(timeStep) < minimumTimeStep) // are we done?
					return true;
			}
//printf("taking step, t: %g, timeStep: %g\n", t, timeStep);
			if(!takeStep(xs, dxdts, t, timeStep, xScalings, actualTimeStep, nextTimeStep))
			{
				return fabs(t2 - t) < minimumTimeStep;
			}
//printf("took step, actualTimeStep: %g, nextTimeStep: %g\n", actualTimeStep, nextTimeStep);

			if(timeStep == actualTimeStep)
				okStepCount++;
			else
				badStepCount++;

			if((t - t2)*(t2 - t1) >= -minimumTimeStep) // are we done?
			{
				//printf("Time steps taken: %i, %i \"good\", %i \"bad\"\n", stepIndex+1, okStepCount, badStepCount);
				return true;
			}
			if(fabs(nextTimeStep) < minimumTimeStep)
			{
				fprintf(stderr, "Time step %g is smaller than the minimum %g in ParticleIntegrator::integrate()\n", fabs(nextTimeStep), minimumTimeStep);
				return false;
			}
			timeStep = nextTimeStep;
		}
		fprintf(stderr, "Too many time steps in ParticleIntegrator::integrate()\n");
		return false;
	}

	inline void computeDerivatives(double t, const UArray<double> & xs, UArray<double> & outDxdts)
	{
		SInt32 vectorSize = xs.getSize()/2;
		UArray<double> x(vectorSize), y(vectorSize);
		UArray<vector2d> vels(vectorSize);
		for(SInt32 index = 0; index < vectorSize; index++)
		{
			x[index] = xs[2*index];
			y[index] = xs[2*index+1];
		}
		velocity.interpolate(x, y, vels);
		for(SInt32 index = 0; index < vectorSize; index++)
		{
			outDxdts[2*index] = vels[index].x;
			outDxdts[2*index+1] = vels[index].y;
		}
	}

private:

	const VectorField2D & velocity;

};

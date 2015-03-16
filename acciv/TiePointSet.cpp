#include "TiePointSet.h"

#include <hdf5.h>
#include <hdf5_hl.h>

#include "ParticleIntegrator.h"
#include "VectorField2D.h"
#include "ScatteredVectorData2D.h"

bool TiePointSet::read(const UString & fileName)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not open file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName;
	hsize_t dimensions[1];
	H5T_class_t typeClass;
	size_t typeSize;
	herr_t status;

	dataSetName = "/x1";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8))
	{
		fprintf(stderr, "TiePointSet::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	SInt32 size = (SInt32)dimensions[0];
	x1.setSize(size);
	y1.setSize(size);
	x2.setSize(size);
	y2.setSize(size);
	correlationCoefficients.setSize(size);

	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, x1.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y1";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "TiePointSet::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, y1.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	dataSetName = "/x2";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "TiePointSet::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, x2.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y2";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "TiePointSet::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, y2.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/correlationCoefficients";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "TiePointSet::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, correlationCoefficients.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/deltaTs";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8))
	{
		fprintf(stderr, "TiePointSet::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	SInt32 deltaTSize = (SInt32)dimensions[0];
	deltaTs.setSize(deltaTSize);
	lowerIndexDeltaT.setSize(deltaTSize);
	upperIndexDeltaT.setSize(deltaTSize);

	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, deltaTs.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/lowerIndexDeltaT";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_INTEGER) || (typeSize != 4) || ((SInt32)dimensions[0] != deltaTs.getSize()))
	{
		fprintf(stderr, "TiePointSet::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_INT, lowerIndexDeltaT.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/upperIndexDeltaT";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_INTEGER) || (typeSize != 4) || ((SInt32)dimensions[0] != deltaTs.getSize()))
	{
		fprintf(stderr, "TiePointSet::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_INT, upperIndexDeltaT.getData());
	if(status < 0)
	{
		fprintf(stderr, "TiePointSet::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5Fclose(fileID);
	return true;
}

bool TiePointSet::write(const UString & fileName) const
{
	if((x1.getSize() != y1.getSize()) || (x1.getSize() != x2.getSize())
		|| (x1.getSize() != y2.getSize()) || (x1.getSize() != correlationCoefficients.getSize()))
	{
		fprintf(stderr, "TiePointSet::write: data vectors are not the same size\n");
		return false;
	}
	hid_t fileID = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "TiePointSet::write: Could not create file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/x1";
	hsize_t dimensions[1];
	dimensions[0] = (hsize_t)x1.getSize();
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, x1.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y1";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, y1.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/x2";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, x2.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y2";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, y2.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/correlationCoefficients";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, correlationCoefficients.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/deltaTs";
	dimensions[0] = (hsize_t)deltaTs.getSize();
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, deltaTs.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/lowerIndexDeltaT";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_INT, lowerIndexDeltaT.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/upperIndexDeltaT";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_INT, upperIndexDeltaT.getData());
	if((status < 0))
	{
		fprintf(stderr, "TiePointSet::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

bool TiePointSet::unadvect(VectorField2D & velocity1, VectorField2D & velocity2, 
	const UArray<double> & deltaT1s, const UArray<double> &  deltaT2s,
	double errorTolerance, SInt32 maximumTimeStepCount)
{

	U_ASSERT(deltaT1s.getSize() == deltaTs.getSize())
	U_ASSERT(deltaT2s.getSize() == deltaTs.getSize())

	ParticleIntegrator integrator1(velocity1),
		integrator2(velocity2);

	for(SInt32 timeIndex = 0; timeIndex < deltaTs.getSize(); timeIndex++)
	{
		double deltaT1 = deltaT1s[timeIndex];
		double deltaT2 = deltaT2s[timeIndex];
		SInt32 lower = lowerIndexDeltaT[timeIndex];
		SInt32 upper = upperIndexDeltaT[timeIndex];
		double nextTimeStep = 0.1*deltaT1;
		double minimumTimeStep = fabs(0.1*deltaT1/maximumTimeStepCount);

		SInt32 size = upper-lower+1;

		UArray<double> points(2*size);
		if(deltaT1 != 0)
		{
			for(SInt32 index = 0; index < size; index++)
			{
				points[2*index] = getX1()[lower+index];
				points[2*index+1] = getY1()[lower+index];
			}
	
			// integrate the pixels backward in time to their locations in the original image
			if(!integrator1.integrate(points, 0, deltaT1, errorTolerance, nextTimeStep, minimumTimeStep, maximumTimeStepCount))
				return false;
	
			// copy the integrated points back
			for(SInt32 index = 0; index < size; index++)
			{
				getX1()[lower+index] = points[2*index];
				getY1()[lower+index] = points[2*index+1];
			}
                }

		if(deltaT2 != 0)
		{
			nextTimeStep = 0.1*deltaT2;
			minimumTimeStep = fabs(0.1*deltaT2/maximumTimeStepCount);
			for(SInt32 index = 0; index < size; index++)
			{
				points[2*index] = getX2()[lower+index];
				points[2*index+1] = getY2()[lower+index];
			}
	
			// integrate the pixels backward in time to their locations in the original image
			if(!integrator2.integrate(points, 0, deltaT2, errorTolerance, nextTimeStep, minimumTimeStep, maximumTimeStepCount))
				return false;
	
			// copy the integrated points back
			for(SInt32 index = 0; index < size; index++)
			{
				getX2()[lower+index] = points[2*index];
				getY2()[lower+index] = points[2*index+1];
			}
		}
	}
	return true;
}

bool TiePointSet::computeCurvedPathVelocities(const VectorField2D & velocity, SInt32 pathVectorCount,
		double errorTolerance, SInt32 maximumTimeStepCount,
		ScatteredVectorData2D & outVelocities) const
{
	ParticleIntegrator integrator(velocity);
	SInt32 size = x1.getSize();

	outVelocities.setSize(pathVectorCount*size);

		
	for(SInt32 timeIndex = 0; timeIndex < deltaTs.getSize(); timeIndex++)
	{
		double deltaT = deltaTs[timeIndex];
		SInt32 lower = lowerIndexDeltaT[timeIndex];
		SInt32 upper = upperIndexDeltaT[timeIndex];
		double timeStep = deltaT/(double)(pathVectorCount-1);
		double nextTimeStep = 0.1*timeStep;
		double minimumTimeStep = fabs(0.1*timeStep/maximumTimeStepCount);
		double t = 0.0;

		SInt32 size = upper-lower+1;
		if(size <= 0)
			continue;

		SInt32 pathSize = 2*pathVectorCount*size;
		UArray<double> points(2*size), v(2*size);
		UArray<double> path1(pathSize),path2(pathSize), pathV1(pathSize), pathV2(pathSize);

		for(SInt32 index = 0; index < size; index++)
		{
			path1[2*index] = points[2*index] = x1[index+lower];
			path1[2*index+1] = points[2*index+1] = y1[index+lower];
		}
		integrator.computeDerivatives(t, points, v);
		for(SInt32 index = 0; index < 2*size; index++)
		{
			pathV1[index] = v[index];
		}
		for(SInt32 pathIndex = 1; pathIndex < pathVectorCount; pathIndex++)
		{
			if(!integrator.integrate(points, t, t+timeStep, errorTolerance, nextTimeStep, minimumTimeStep, maximumTimeStepCount))
				return false;
			t += timeStep;
			integrator.computeDerivatives(t, points, v);
			for(SInt32 index = 0; index < size; index++)
			{
				SInt32 pIndex = 2*index + pathIndex*2*size;
				path1[pIndex] = points[2*index];
				path1[pIndex+1] = points[2*index+1];
				pathV1[pIndex] = v[2*index];
				pathV1[pIndex+1] = v[2*index+1];
			}
		}

		t = deltaT;
		for(SInt32 index = 0; index < size; index++)
		{
			SInt32 pIndex = 2*index + (pathVectorCount-1)*2*size;
			path2[pIndex] = points[2*index] = x2[index+lower];
			path2[pIndex+1] = points[2*index+1] = y2[index+lower];
		}
		integrator.computeDerivatives(t, points, v);
		for(SInt32 index = 0; index < 2*size; index++)
		{
			SInt32 pIndex = index + (pathVectorCount-1)*2*size;
			pathV2[pIndex] = v[index];
		}
		for(SInt32 pathIndex = pathVectorCount-2; pathIndex >= 0; pathIndex--)
		{
			double timeStep = deltaT/(double)(pathVectorCount-1);
			if(!integrator.integrate(points, t, t-timeStep, errorTolerance, nextTimeStep, minimumTimeStep, maximumTimeStepCount))
				return false;
			t -= timeStep;
			integrator.computeDerivatives(t, points, v);
			for(SInt32 index = 0; index < size; index++)
			{
				SInt32 pIndex = 2*index + pathIndex*2*size;
				path2[pIndex] = points[2*index];
				path2[pIndex+1] = points[2*index+1];
				pathV2[pIndex] = v[2*index];
				pathV2[pIndex+1] = v[2*index+1];
			}
		}

		for(SInt32 pathIndex = 0; pathIndex < pathVectorCount; pathIndex++)
		{
			double alpha = pathIndex/(double)(pathVectorCount-1);
			//double alpha = 0;
			for(SInt32 index = 0; index < size; index++)
			{
				SInt32 pIndex = 2*index + pathIndex*2*size;
				SInt32 vIndex = pathIndex + pathVectorCount*(index+lower);
				double x = (1-alpha)*path1[pIndex]
					+ alpha*path2[pIndex];
				double y = (1-alpha)*path1[pIndex+1]
					+ alpha*path2[pIndex+1];
				outVelocities.getX(vIndex) = x;
				outVelocities.getY(vIndex) = y;
				vector2d v;
				v.x = (1-alpha)*pathV1[pIndex]
					+ alpha*pathV2[pIndex]
					+ (path2[pIndex] - path1[pIndex])/deltaT;
				v.y = (1-alpha)*pathV1[pIndex+1]
					+ alpha*pathV2[pIndex+1]
					+ (path2[pIndex+1] - path1[pIndex+1])/deltaT;
				outVelocities.getData(vIndex) = v;
			}
		}
	}

	return true;
}

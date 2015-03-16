#include "ScalarField2D.h"
#include <platform/UStreamFile.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <math.h>
#include <algorithm>

//#include <fftw3.h>
//
//#include <complex>

using namespace std;

bool ScalarField2D::read(const UString & fileName)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ScalarField2D::read: Could not open file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/data";
	hsize_t dimensions[2];
	H5T_class_t typeClass;
	size_t typeSize;
	herr_t status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8))
	{
		fprintf(stderr, "ScalarField2D::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	ySize = (SInt32)dimensions[0];
	xSize = (SInt32)dimensions[1];

	data.setSize(xSize*ySize);
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, data.getData());
	if(status < 0)
	{
		fprintf(stderr, "ScalarField2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/bounds";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || (dimensions[0] != 4))
	{
		fprintf(stderr, "ScalarField2D::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	UArray<double> bounds(4);
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, bounds.getData());
	if(status < 0)
	{
		fprintf(stderr, "ScalarField2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	xLower = bounds[0];
	xUpper = bounds[1];
	yLower = bounds[2];
	yUpper = bounds[3];

	x.setSize(xSize);
	y.setSize(ySize);
	setXandY();
	status = H5Fclose(fileID);

	return true;
}

bool ScalarField2D::write(const UString & fileName) const
{
	hid_t fileID = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ScalarField2D::write: Could not create file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/data";
	hsize_t dimensions[2];
	dimensions[0] = (hsize_t)ySize;
	dimensions[1] = (hsize_t)xSize;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 2, dimensions, H5T_NATIVE_DOUBLE, data.getData());
	if((status < 0))
	{
		fprintf(stderr, "ScalarField2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/bounds";
	UArray<double> bounds(4);
	bounds[0] = xLower;
	bounds[1] = xUpper;
	bounds[2] = yLower;
	bounds[3] = yUpper;
	dimensions[0] = (hsize_t)4;
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, bounds.getData());
	if((status < 0))
	{
		fprintf(stderr, "ScalarField2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

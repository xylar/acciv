#include "VectorField2D.h"

#include <hdf5.h>
#include <hdf5_hl.h>

bool VectorField2D::read(const UString & fileName)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "VectorField2D::read: Could not open file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/dataX";
	hsize_t dimensions[2];
	H5T_class_t typeClass;
	size_t typeSize;
	herr_t status;

	dataSetName = "/dataX";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8))
	{
		fprintf(stderr, "VectorField2D::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	SInt32 ySize = (SInt32)dimensions[0];
	SInt32 xSize = (SInt32)dimensions[1];


	dataSetName = "/bounds";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || (dimensions[0] != 4))
	{
		fprintf(stderr, "VectorField2D::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	UArray<double> bounds(4);
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, bounds.getData());
	if(status < 0)
	{
		fprintf(stderr, "VectorField2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	initialize(xSize,ySize, bounds[0],bounds[1],bounds[2],bounds[3]);
    UArray<double> data(xSize*ySize);

	dataSetName = "/dataX";
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, data.getData());
	if(status < 0)
	{
		fprintf(stderr, "VectorField2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		(*this)[index].x = data[index];
	}

	dataSetName = "/dataY";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8))
	{
		fprintf(stderr, "VectorField2D::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	if(((SInt32)dimensions[0] != ySize) || ((SInt32)dimensions[1] != xSize))
	{
		fprintf(stderr, "VectorField2D::read: dataX is not the same size as dataY in file %s\n", (const char *)fileName);
		return false;
	}

	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, data.getData());
	if(status < 0)
	{
		fprintf(stderr, "VectorField2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5Fclose(fileID);

	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		(*this)[index].y = data[index];
	}


	return true;
}

bool VectorField2D::write(const UString & fileName) const
{
	hid_t fileID = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "VectorField2D::write: Could not create file %s\n", (const char *)fileName);
		return false;
	}

    UArray<double> data(xSize*ySize);
	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		 data[index] = (*this)[index].x;
	}

    UString dataSetName;
	hsize_t dimensions[2];
	dimensions[0] = (hsize_t)getYSize();
	dimensions[1] = (hsize_t)getXSize();
	herr_t status;
	dataSetName = "/dataX";
	status = H5LTmake_dataset(fileID, dataSetName, 2, dimensions, H5T_NATIVE_DOUBLE, data.getData());
	if((status < 0))
	{
		fprintf(stderr, "VectorField2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	for(UInt32 index = 0; index < data.getSize(); index++)
	{
		 data[index] = (*this)[index].y;
	}

	dataSetName = "/dataY";
	status = H5LTmake_dataset(fileID, dataSetName, 2, dimensions, H5T_NATIVE_DOUBLE, data.getData());
	if((status < 0))
	{
		fprintf(stderr, "VectorField2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/bounds";
	UArray<double> bounds(4);
	bounds[0] = getXLower();
	bounds[1] = getXUpper();
	bounds[2] = getYLower();
	bounds[3] = getYUpper();
	dimensions[0] = (hsize_t)4;
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, bounds.getData());
	if((status < 0))
	{
		fprintf(stderr, "VectorField2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

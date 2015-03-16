
#include "ScatteredScalarData2D.h"
#include <platform/UStreamFile.h>

#include <hdf5.h>
#include <hdf5_hl.h>

bool ScatteredScalarData2D::read(UString fileName)
{
	hid_t fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ScatteredData2D::read: Could not open file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/data";
	hsize_t dimensions[1];
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
		fprintf(stderr, "ScatteredData2D::read: Could not get dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	size = (SInt32)dimensions[0];
	x.setSize(size);
	y.setSize(size);
	data.setSize(size);

	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, data.getData());
	if(status < 0)
	{
		fprintf(stderr, "ScatteredData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/x";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "ScatteredData2D::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, x.getData());
	if(status < 0)
	{
		fprintf(stderr, "ScatteredData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y";
	status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	if(status < 0)
	{
		dataSetName.makeUpper();
		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
	}
	if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
	{
		fprintf(stderr, "ScatteredData2D::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}
	status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, y.getData());
	if(status < 0)
	{
		fprintf(stderr, "ScatteredData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/weights";
	bool weightsFound = H5LTfind_dataset(fileID, dataSetName);
	if(!weightsFound)
	{
		dataSetName.makeUpper();
		weightsFound = H5LTfind_dataset(fileID, dataSetName);
	}
	if(weightsFound)
	{

		status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
		if(status < 0)
		{
			dataSetName.makeUpper();
			status = H5LTget_dataset_info(fileID, dataSetName, dimensions, &typeClass, &typeSize);
		}
		if((status < 0) || (typeClass != H5T_FLOAT) || (typeSize != 8) || ((SInt32)dimensions[0] != size))
		{
			fprintf(stderr, "ScatteredData2D::read: Problem with dataset info for %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
			return false;
		}
		status = H5LTread_dataset(fileID, dataSetName, H5T_NATIVE_DOUBLE, weights.getData());
		if(status < 0)
		{
			fprintf(stderr, "ScatteredData2D::read: Could not read dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
			return false;
		}
	}
	else
	{
		for(UInt32 index = 0; index < getSize(); index++)
		{
			weights[index] = 1.0;
		}
	}

	status = H5Fclose(fileID);

	return true;
}

bool ScatteredScalarData2D::write(UString fileName)
{
	hid_t fileID = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "ScatteredData2D::write: Could not create file %s\n", (const char *)fileName);
		return false;
	}

	UString dataSetName = "/data";
	hsize_t dimensions[1];
	dimensions[0] = (hsize_t)size;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, data.getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/x";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, x.getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/y";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, y.getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	dataSetName = "/weights";
	status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, getWeights().getData());
	if((status < 0))
	{
		fprintf(stderr, "ScatteredVectorData2D::write: Could not make dataset %s in file %s\n", (const char *)dataSetName, (const char *)fileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}

// smoothFit.cpp : performs a smooth fit to scattered data
//

#include <core/UTypes.h>
#include <containers/UArray.h>
#include <core/UString.h>

#include "VectorField2D.h"
#include "ScatteredVectorData2D.h"
#include "ParameterFileReader.h"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <math.h>
#include <cstdlib>


bool getParameters();

bool getGridInfo();

bool constructGriddedFieldFromScatteredField();

UString inGridFileName;
UString inScatteredFieldFileName;
UString outGridFieldFileName;

UArray<SInt32> imageSize;
UArray<double> imageBounds;

SInt32 smoothFitMinControlPointScatteredNeighbors;


int main(int argc, char * argv[])
{

	if(!getParameters())
	{
		fprintf(stderr, "Could not read parameters\n");
		return -1;
	}

	if(!getGridInfo())
	{
		fprintf(stderr, "Could not read grid file\n");
		return -1;
	}

	if(!constructGriddedFieldFromScatteredField())
	{
		fprintf(stderr, "Could not perform smooth fit\n");
		return -1;
	}


	return 0;
}


bool getParameters()
{
	UString fileName = "parameters.ascii";
	ParameterFileReader reader(fileName);

	UString defaultParameterFile;
	UString tag;
	tag = "defaultParameterFile";
	if(reader.getString(tag, defaultParameterFile))
	{
		reader.addFile(defaultParameterFile);
	}

	tag = "inGridFileName";
	if(!reader.getString(tag, inGridFileName))
	{
		fprintf(stderr, "getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	tag = "inScatteredFieldFileName";
	if(!reader.getString(tag, inScatteredFieldFileName))
	{
		fprintf(stderr, "getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	tag = "outGridFieldFileName";
	if(!reader.getString(tag, outGridFieldFileName))
	{
		fprintf(stderr, "getParameters: Could not read tag %s\n", (const char *)tag);
		return false;
	}
	tag = "smoothFitMinControlPointScatteredNeighbors";
	if(!reader.getInteger(tag, smoothFitMinControlPointScatteredNeighbors))
	{
		smoothFitMinControlPointScatteredNeighbors = 16;
	}
	return true;
}

bool getGridInfo()
{

	VectorField2D grid;

	if(!grid.read(inGridFileName))
	{
		fprintf(stderr, "getGridInfo: Could not read grid velocity scaling from file %s\n", (const char *)inGridFileName);
		return false;
	}

	imageSize.setSize(2);
	imageSize[0] = grid.getXSize();
	imageSize[1] = grid.getYSize();

	imageBounds.setSize(4);
	imageBounds[0] = grid.getXLower();
	imageBounds[1] = grid.getXUpper();
	imageBounds[2] = grid.getYLower();
	imageBounds[3] = grid.getYUpper();

	return true;
}

bool constructGriddedFieldFromScatteredField()
{

	ScatteredVectorData2D scatteredField;
	if(!scatteredField.read(inScatteredFieldFileName))
		return false;

	VectorField2D griddedField(imageSize[0], imageSize[1], imageBounds[0], imageBounds[1], imageBounds[2], imageBounds[3]);

	double meanRes = 0.0;
	UArray<vector2d> residuals;

	scatteredField.smoothFit(griddedField, smoothFitMinControlPointScatteredNeighbors, residuals, meanRes);

	if(!griddedField.write(outGridFieldFileName))
		return false;

	// add the mean residual
	hid_t fileID = H5Fopen(outGridFieldFileName, H5F_ACC_RDWR, H5P_DEFAULT);
	if(fileID < 0)
	{
		fprintf(stderr, "constructGriddedFieldFromScatteredField: Could not open file %s\n",
			(const char *)outGridFieldFileName);
		return false;
	}

	UString dataSetName = "/meanResidual";
	hsize_t dimensions[2];
	dimensions[0] = 1;
	herr_t status = H5LTmake_dataset(fileID, dataSetName, 1, dimensions, H5T_NATIVE_DOUBLE, &meanRes);
	if((status < 0))
	{
		fprintf(stderr, "constructGriddedFieldFromScatteredField: Could not make dataset %s in file %s\n",
			(const char *)dataSetName, (const char *)outGridFieldFileName);
		return false;
	}

	status = H5Fclose(fileID);

	return true;
}


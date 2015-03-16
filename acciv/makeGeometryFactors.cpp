// makeGeometryFactors.cpp : generates a gridGeometryFactors.h5 file given the bounds and image size from an HDF5 image file
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "VectorField2D.h"
#include "MaskedImage.h"
#include <core/UString.h>

int main(int argc, char * argv[])
{
	if(argc < 3 || argc > 6)
	{
		fprintf(stderr, "example usage: %s image001.h5 gridGeometryFactors.h5 [centric] [Saturn] [flat] [Re=7e7 Rp=6e7]\n", argv[0]);
		fprintf(stderr, "default is planetographic latitude and Jovian radii\n");
		return 0;
	}
	char * imageFileName = argv[1];
	char * geometryFactorsFileName = argv[2];
	
	bool flat = false, centric = false, Saturn = false, radiiGiven = false;

	double Re = 7.1492e7; //71,492 km
	double Rp = 6.6854e7; //66,854 km

	for(UInt32 index = 3; index < argc; index++)
	{
		UString arg(argv[index]);
		if(arg.compareNoCase("flat") == 0)
		{
			flat = true;
			break;
		}
	}
	if(flat)
	{
		printf("Using flat geometry\n");
	}
	else
	{
		for(UInt32 index = 3; index < argc; index++)
		{
			UString arg(argv[index]);
			if(arg.compareNoCase("centric") == 0)
			{
				centric = true;
				break;
			}
		}
		if(centric)
			printf("Using planetocentric latitude\n");
		else
			printf("Using planetographic latitude\n");

		for(UInt32 index = 3; index < argc; index++)
		{
			UString line(argv[index]);
			SInt32 breakIndex = line.find('=');
			if(breakIndex != -1)
			{
				UString tag = line.left(breakIndex);
				while((tag.getLength() > 0) && (tag.getCharacterAt(tag.getLength()-1) == ' '))
					tag = tag.left(tag.getLength()-1);
				if(tag.compareNoCase("Re") == 0)
				{
					UString value = line.mid(breakIndex+1);
					while((value.getLength() > 0) && (value.getCharacterAt(0) == ' '))
						value = value.mid(1);
					Re = atof(value);
					radiiGiven = true;
				}
				else if(tag.compareNoCase("Rp") == 0)
				{
					UString value = line.mid(breakIndex+1);
					while((value.getLength() > 0) && (value.getCharacterAt(0) == ' '))
						value = value.mid(1);
					Rp = atof(value);
					radiiGiven = true;
				}
			}
		}
		if(radiiGiven)
		{
			printf("Using radii Re=%g and Rp=%g\n", Re, Rp);
		}
		else
		{
			for(UInt32 index = 3; index < argc; index++)
			{
				UString arg(argv[index]);
				if(arg.compareNoCase("Saturn") == 0)
				{
					Saturn = true;
					break;
				}
			}
			if(Saturn)
			{
				// from Wikipedia
				// Equatorial radius 	60,268 � 4 km
				// Polar radius 		54,364 � 10 km
				Re = 6.0268e7;
				Rp = 5.4364e7;
				printf("Using Saturn's radii\n");
			}
			else
			{
				printf("Using Jupiter's radii\n");
			}
		}
	}

	MaskedImage image;
	if(!image.read(imageFileName))
	{
		fprintf(stderr, "Could not read image file %s", imageFileName);
		return 0;
	}

	double epsilon = Rp/Re;
	double epsilon2 = epsilon*epsilon;
	const double pi = 3.14159265358979323846264338327950;
	const double piOver180 = pi/180.0;

	UInt32 nx = image.getXSize();
	UInt32 ny = image.getYSize();

	VectorField2D geometryFactors(nx, ny, image.getXLower(),
			image.getXUpper(), image.getYLower(), image.getYUpper());

	for(UInt32 yIndex = 0; yIndex < ny; yIndex++)
	{
		if(flat)
		{
			for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
			{
				UInt32 index = xIndex + nx*yIndex;
				geometryFactors[index].x = 1.0;
				geometryFactors[index].y = 1.0;
			}
		}
		else
		{
			double lat = piOver180*image.getY()[yIndex];
			// convert to planetographic latitude
			if(centric)
				lat = atan(tan(lat)/epsilon2);
			double secLat = 1.0/cos(lat);
			double tanLat = tan(lat);
			double t = atan(epsilon*tanLat);
			double dtdLatitude = piOver180*epsilon*(secLat*secLat/(1 + epsilon2*tanLat*tanLat));
			double dydLatitude = sqrt(Re*Re*sin(t)*sin(t) + Rp*Rp*cos(t)*cos(t))*dtdLatitude;
			double dxdLongitude = -Re*cos(t)*piOver180;
			for(UInt32 xIndex = 0; xIndex < nx; xIndex++)
			{
				UInt32 index = xIndex + nx*yIndex;
				geometryFactors[index].x = dxdLongitude;
				geometryFactors[index].y = dydLatitude;
			}
		}
	}
	if(!geometryFactors.write(geometryFactorsFileName))
	{
		fprintf(stderr, "Could not write geometry factors file %s", geometryFactorsFileName);
		return 0;
	}

	return 0;
}


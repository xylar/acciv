#pragma once

#include <core/UTypes.h>
#include <containers/UArray.h>

// a quad tree for use in computing a smooth fit to scattered data

class QuadTreeNode
{
public:
	UArray<UInt32> indices;
	UArray<double> x;
	UArray<double> y;
};

class QuadTreeLevel
{
public:
	QuadTreeLevel()
		:nx(0),ny(0)
	{
	}

	// Construct the base of a quad tree from scattered data
	QuadTreeLevel(const UArray<double> & x, const UArray<double> & y, UInt32 inNx, UInt32 inNy,
		double xLower, double xUpper, double yLower, double yUpper );

	// Construct this level of a quad tree by bisecting the previous level
	void bisect();

	UArray<QuadTreeNode> nodes;
	UInt32 nx, ny;
};

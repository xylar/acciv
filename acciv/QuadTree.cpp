#include "QuadTree.h"
#include <math.h>
#include <algorithm>

// Construct the base of a quad tree from scattered data
QuadTreeLevel::QuadTreeLevel(const UArray<double> & x, const UArray<double> & y, UInt32 inNx, UInt32 inNy,
		double xLower, double xUpper, double yLower, double yUpper )
	: nx(inNx), ny(inNy)
{
	nodes.setSize(nx*ny);
	double deltaX = (xUpper-xLower)/nx;
	double deltaY = (yUpper-yLower)/ny;
	for(UInt32 index = 0; index < x.getSize(); index++)
	{
		double xFrac = (x[index]-xLower)/deltaX;
		double yFrac = (y[index]-yLower)/deltaY;
		UInt32 xIndex = std::min(std::max((SInt32)(xFrac),0),(SInt32)nx-1);
		UInt32 yIndex = std::min(std::max((SInt32)(yFrac),0),(SInt32)ny-1);
		xFrac -= xIndex;
		yFrac -= yIndex;
		UInt32 nodeIndex = xIndex+nx*yIndex;
		nodes[nodeIndex].indices.add(index);
		nodes[nodeIndex].x.add(xFrac);
		nodes[nodeIndex].y.add(yFrac);

	}
}

// Construct this level of a quad tree by bisecting the previous level
void QuadTreeLevel::bisect()
{
	UArray<QuadTreeNode> nodes_old(nodes);
	nx = 2*nx;
	ny = 2*ny;
	nodes.setSize(0);
	nodes.setSize(nx*ny);
	for(UInt32 yIndex = 0; yIndex < ny/2; yIndex++)
	{
		for(UInt32 xIndex = 0; xIndex < nx/2; xIndex++)
		{
			UInt32 inBinIndex = xIndex+nx/2*yIndex;
			const UArray<UInt32> & indices = nodes_old[inBinIndex].indices;
			const UArray<double> & x = nodes_old[inBinIndex].x;
			const UArray<double> & y = nodes_old[inBinIndex].y;
			for(UInt32 dataIndex = 0; dataIndex < indices.getSize(); dataIndex++)
			{
				UInt32 xOffset = (x[dataIndex] > 0.5) ? 1 : 0;
				UInt32 yOffset = (y[dataIndex] > 0.5) ? 1 : 0;
				UInt32 binIndex = 2*xIndex+xOffset + nx*(2*yIndex+yOffset);
				nodes[binIndex].indices.add(indices[dataIndex]);
				nodes[binIndex].x.add(2*x[dataIndex]-xOffset);
				nodes[binIndex].y.add(2*y[dataIndex]-yOffset);
			}
		}
	}
}

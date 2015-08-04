#!/bin/csh

echo \# Parameters for ACCIV used construct a velocity field from a set of earlier 
echo \# and a set of later images \(which may be the same: an image will never be 
echo \# compared with itself and each image will only be used once\).
echo \# 
echo \# default parameters file, used to assign parameters that don\'t change between runs
echo defaultParameterFile = ..\/defaultParameters.ascii
echo 
echo \# the number of images in the earlier image set
echo earlierImageIndices = \[1,2,3,4\]
echo 
echo \# the number of images in the later image set
echo laterImageIndices = \[1,2,3,4\]
echo 
echo \# do we want to advect the images using a velocity field we already computed\?
echo \#   If this is the first ACCIV pass, there will be no velocity field to use yet,
echo \#   so this flag should be set to false.  Otherwise it should be set to true.
echo advectImages = $1
echo 
echo \# the velocity data \(if any\) to be used to advect the earlier images to a common time
echo inEarlierVelocityFileName = ..\/$2\/outGridVelocity.h5

echo \# the velocity data \(if any\) to be used to advect the later images to a common time
echo inLaterVelocityFileName = ..\/$2\/outGridVelocity.h5
echo 
echo \# the size of the box of pixels used to perform correlations: \[xSize, ysize\]
echo correlationBoxSize = \[$3, $3\]
echo 
echo \# the range of pixels over which to search for correlations: \[xMin, xMax, yMin, yMax\]
echo \# This search range will be scaled by the time separation between image
echo \# pairs \(normalized by the max time between image pairs\).  This insures that ACCIV
echo \# searches the same approximate range of velocities for each image pair
echo searchRange = \[-$4, $4, -$5, $5\]
echo 
echo \# the stride between correlations \(after finding each correlation, how many pixels
echo \# do we shift the correlation box by to find the next correlation\): \[xStride, yStride\]
echo stride = \[$6,$6\]
echo 
echo \# the minimum time between images, below which image pairs will not be correlated
echo minimumTimeSeparation = $7

echo \# the minimum number of scattered data points that neighbor a given b-spline control point, 
echo \# below which the control point is dropped from the fit
echo smoothFitMinControlPointScatteredNeighbors = $8

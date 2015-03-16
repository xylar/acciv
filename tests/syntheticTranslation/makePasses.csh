#!/bin/csh

# usage:
# ./writeParameterFile.csh writeImages previousPass correlationBoxSize searchRangeX searchRangeY stride minimumTimeSeparation smoothFitMinControlPointScatteredNeighbors > currentPass/parameters.ascii
mkdir pass01
../syntheticTestScripts/writeParameterFile.csh false pass01 80 110 55 4 0 16 > pass01/parameters.ascii
mkdir pass02
../syntheticTestScripts/writeParameterFile.csh true pass01 80 20 20 4 0 16 > pass02/parameters.ascii
mkdir pass03
../syntheticTestScripts/writeParameterFile.csh true pass02 80 20 20 4 0 16 > pass03/parameters.ascii
mkdir pass04
../syntheticTestScripts/writeParameterFile.csh true pass03 80 20 20 4 0 16 > pass04/parameters.ascii

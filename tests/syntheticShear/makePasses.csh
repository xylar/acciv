#!/bin/csh

# usage:
# ./writeParameterFile.csh writeImages previousPass correlationBoxSize searchRangeX searchRangeY stride minimumTimeSeparation smoothFitMinControlPointScatteredNeighbors > currentPass/parameters.ascii
mkdir -p pass01
../syntheticTestScripts/writeParameterFile.csh false pass01 120 75 10 4 0 32 > pass01/parameters.ascii
mkdir -p pass02
../syntheticTestScripts/writeParameterFile.csh true pass01 80 20 20 4 0 32 > pass02/parameters.ascii
mkdir -p pass03
../syntheticTestScripts/writeParameterFile.csh true pass02 80 20 20 4 0 32 > pass03/parameters.ascii
mkdir -p pass04
../syntheticTestScripts/writeParameterFile.csh true pass03 80 20 20 4 0 32 > pass04/parameters.ascii
mkdir -p pass05
../syntheticTestScripts/writeParameterFile.csh true pass04 80 20 20 4 0 32 > pass05/parameters.ascii
mkdir -p pass06
../syntheticTestScripts/writeParameterFile.csh true pass05 80 20 20 4 0 32 > pass06/parameters.ascii

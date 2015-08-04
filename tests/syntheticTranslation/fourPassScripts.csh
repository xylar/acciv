#!/bin/csh
../../acciv/makeGeometryFactors.exe image001.h5 gridGeometryFactors.h5 flat
../syntheticTestScripts/serialPasses.csh pass01 pass02 pass03 pass04

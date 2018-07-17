#!/bin/csh
../../acciv/makeGeometryFactors image001.h5 gridGeometryFactors.h5 flat
../syntheticTestScripts/serialPasses.csh day1/pass1 day1/pass2 day2/pass1 day2/pass2 days_1_2/pass1 days_1_2/pass2

#!/bin/csh
../../acciv/makeGeometryFactors.exe image001.h5 gridGeometryFactors.h5 centric
./serialPasses.csh later/pass1 later/pass2 later/pass3 full/pass1 full/pass2 full/pass3

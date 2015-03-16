#!/bin/csh
setenv acciv '../../../../acciv/acciv.exe'
../../acciv/makeGeometryFactors.exe image001.h5 gridGeometryFactors.h5 centric
cd later/pass1
date
echo  later/pass1
$acciv
cd ../pass2
date
echo  later/pass2
$acciv
cd ../pass3
date
echo  later/pass3
$acciv
cd ../../full/pass1
date
echo  full/pass1
$acciv
cd ../pass2
date
echo  full/pass2
$acciv
cd ../pass3
date
echo  full/pass3
$acciv
cd ../..
date

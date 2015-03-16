#!/bin/csh
setenv acciv '../../acciv/acciv.exe'
echo $1
stdbuf -oL $acciv $1 > $1/output.log
python ../syntheticTestScripts/plotVelocities.py --folder=$1 --imageFileName=image001.h5 --savePlots
python ../syntheticTestScripts/plotAdvectedImages.py --inFolder=$1/_work --outFolder=$1
../syntheticTestScripts/plotWorkVelocities.csh $1/_work

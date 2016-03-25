#!/bin/bash
shopt -s expand_aliases

which acciv
if [ "$?" != "0" ]; then
    alias acciv="../../../../acciv/acciv.exe"
fi

which acciv-makeGeometryFactors
if [ "$?" != "0" ]; then
    alias acciv-makeGeometryFactors="../../acciv/makeGeometryFactors.exe"
fi

set -e

acciv-makeGeometryFactors images/image000.h5 images/gridGeometryFactors.h5 flat

cd passes/pass1
acciv
cd ../..

python3 plotVelocities.py --folder passes/pass1 --imageFileName images/image000.h5

cd passes/pass2
acciv
cd ../..

python3 plotVelocities.py --folder passes/pass2 --imageFileName images/image000.h5 

#!/bin/bash

acciv-makeGeometryFactors images/image000.h5 images/gridGeometryFactors.h5 flat

cd passes/pass1
acciv
cd ../..

python3 plotVelocities.py --folder passes/pass1 --imageFileName images/image000.h5 --plotSettings plotSettings.json

cd passes/pass2
acciv
cd ../..

python3 plotVelocities.py --folder passes/pass2 --imageFileName images/image000.h5 --plotSettings plotSettings.json
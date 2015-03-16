#!/bin/csh
python ../syntheticTestScripts/plotResiduals.py --folder=$1 --scatterFileName=noCurvedPathScatteredVelocity.h5 --figurePrefix=noCurved --savePlots
foreach pass (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
  python ../syntheticTestScripts/plotResiduals.py --folder=$1 --scatterFileName=curvedPathScatteredVelocity_$pass.h5 --figurePrefix=curved_$pass --savePlots
end

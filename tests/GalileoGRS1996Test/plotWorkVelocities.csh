#!/bin/csh
python plotVelocities.py --folder=$1 --imageFileName=image001.h5 --scatterFileName=noCurvedPathScatteredVelocity.h5 --gridFileName=noCurvedPathGridVelocity.h5 --figurePrefix=noCurved --savePlots --tiePointsFolder=.
foreach pass (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
  python plotVelocities.py --folder=$1 --imageFileName=image001.h5 --scatterFileName=curvedPathScatteredVelocity_$pass.h5 --gridFileName=curvedPathGridVelocity_$pass.h5 --figurePrefix=curved_${pass}_ --savePlots --tiePointsFolder=.

end

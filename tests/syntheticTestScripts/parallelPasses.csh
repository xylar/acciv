#!/bin/csh
foreach pass ($*)
  echo $pass
  ../syntheticTestScripts/passScript.csh $pass > /dev/null &
end

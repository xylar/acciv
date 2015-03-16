#!/bin/csh
foreach pass ($*)
  ../syntheticTestScripts/passScript.csh $pass
end

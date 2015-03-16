#!/bin/csh
foreach pass ($*)
  ./passScript.csh $pass
end

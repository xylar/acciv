#!/bin/csh
foreach pass ($*)
  echo $pass
  ./passScript.csh $pass > /dev/null &
end

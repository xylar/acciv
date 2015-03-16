#!/bin/csh

foreach name ( *.ascii */*/parameters.ascii )
  echo $name
  cat $name | tr -d '\r' > temp
  mv temp $name
end

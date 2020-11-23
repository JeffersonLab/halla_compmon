#!/bin/bash

start=0
end=0
if [ "$1" -eq 1 ]; then
  start=1
  end=44
elif [ "$1" -eq 2 ]; then
  start=101
  end=223
else
  echo 
  echo "Please enter proper run code. PREX==1, CREX==2" 
  echo 
  exit 1;
fi

for ((i=$start; i<=$end; i++))
do
  ./aggregate.sh -s $i --nogrand
done

$COMPMON_GRAND/exptPlots.sh $1

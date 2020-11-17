#!/bin/bash

start=0
end=0
if [ "$1" -eq 1 ]; then
  start=4300
  end=4620
elif [ "$1" -eq 2]; then
  start=4960
  end=6240
else
  echo 
  echo "Please enter proper run code. PREX==1, CREX==2" 
  echo 
  exit 1;
fi

for i in {$start..$end..20}
do
  ./groupRunPlots.sh $i
  sleep 10m
done


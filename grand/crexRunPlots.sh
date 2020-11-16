#!/bin/bash

for i in {4960..6240..20}
do
  ./groupRunPlots.sh $i
  sleep 10m
done


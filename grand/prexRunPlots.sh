#!/bin/bash

for i in {4300..4620..20}
do
  ./groupRunPlots.sh $i
  sleep 10m
done


#!/bin/bash

for i in {4960..6240..20}
do
  ./group_online.sh $i
  sleep 5m
done

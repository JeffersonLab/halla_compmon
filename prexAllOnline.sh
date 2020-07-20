#!/bin/bash

for i in {4340..4620..20}
do
  ./group_online.sh $i
  sleep 20m
done

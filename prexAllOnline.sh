#!/bin/bash

for i in {4300..4620..20}
do
  ./group_online.sh $i
  sleep 30m
done

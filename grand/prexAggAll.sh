#!/bin/bash

for i in {1..44}
do
  ./aggregate.sh -s $i --nogrand
done

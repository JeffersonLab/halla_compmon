#!/bin/bash

for i in {101..223}
do
  ./aggregate.sh -s $i --nogrand
done

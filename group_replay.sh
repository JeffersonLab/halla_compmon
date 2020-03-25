#!/bin/bash

for ((run = $1; run <= $2; run++)); do
  ./compmon.sh -r $run &
done

#!/bin/bash

for ((run = $1; run <= $2; run++)); do
  ./online.sh -r $run --nopanguin --webupload &
done

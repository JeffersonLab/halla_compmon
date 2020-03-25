#!/bin/bash

add=20
end=$(( $1 + $add ))

for ((run = $1; run <= $end; run++)); do
  ./online.sh -r $run &
done

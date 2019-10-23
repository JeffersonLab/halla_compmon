#!/bin/bash

for ((snail = $1; snail <= $2; snail++)); do
  #root -l -b -q cyclePlot.C'("snail'$snail'")' &
  root -l -b -q aggregate.C'("snail'$snail'", 0, 1)' &
done

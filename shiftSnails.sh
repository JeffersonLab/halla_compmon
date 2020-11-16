#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Illegal number of parameters"
  echo "Usage: ./shiftSnails.sh <snail_min> <snail_max> <shift>"
fi

snail_min=$1
snail_max=$2
shift=$3

if [ $shift > 0 ]; then
  for ((snail=$snail_max; $snail >= $snail_min; snail--)); do
    echo "Moving file snail$snail.list to snail"$(( $snail + $shift ))".list"
    #mv snails/snail$snail.list snails/snail$(( $snail + $shift )).list
    mv $COMPMON_WEB/snails/snail$snail $COMPMON_WEB/snails/snail$(( $snail + $shift ))
  done
else
  for ((snail=$snail_min; $snail <= $snail_min; snail++)); do
    echo "Moving file snail$snail.list to snail"$(( $snail + $shift ))".list"
    #mv snails/snail$snail.list snails/snail$(( $snail + $shift )).list
    mv $COMPMON_WEB/snails/snail$snail $COMPMON_WEB/snails/snail$(( $snail + $shift ))
  done
fi

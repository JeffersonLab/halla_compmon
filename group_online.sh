#!/bin/bash

exclude_list=("4294" "4338" "4380" "4429" "4479" "4481" \
"4498" "4500" "4511" "4516" "4566" "4567" "4582" "4587" \
"4588" "4619" "4947" "4948" "5064" "5195" "5278" "5310" \
"5315" "5387" "5398" "5424" "5425" "5475" "5486" "5488" \
"5530" "5546" "5582" "5626" "5637" "5641" "5655" "5657" \
"5718" "5719" "5790")

add=19
end=$(( $1 + $add ))

for ((run = $1; run <= $end; run++)); do
  if [[ ! "${exclude_list[@]}" =~ "${run}" ]];
  then
    ./online.sh -r $run &
  fi
done

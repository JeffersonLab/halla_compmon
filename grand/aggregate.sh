#!/bin/bash

usage() { echo "Usage: $0 [-s <snail number>] [-h, print help info] [--nogrand Optional: skip the grand rootfile build] [--exptPlots make the snail/run/cyclewise plots for the entire experimental run]" 1>&2; exit 1; }

#while getopts ":h:s:" opt; do
#  case ${opt} in
#    h )
#      usage
#      ;;
#    s )
#      snail_num=$OPTARG
#      ;;
#  esac
#done

do_grand=1
do_expt_plots=0

while (( "$#" )); do
  case "$1" in
    -s|--snail)
      snail_num=$2
      shift 2
      ;;
    --nogrand)
      do_grand=0
      shift 1
      ;;
    --exptPlots)
      do_expt_plots=1
      shift 1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

if [ -z "${snail_num}" ]; then
  usage
fi

if [ ! -d $COMPMON_WEB/snails/snail$snail_num ]; then
  mkdir $COMPMON_WEB/snails/snail$snail_num
fi

run_code=2
if [ "$snail_num" -lt 100 ]; then
  run_code=1
fi

#root -l -b -q miniruns.C\(\"snail${snail_num}\"\)
#root -l -b -q laserCycles.C\(${snail_num}\)
#root -l -b -q aggregate.C\(\"snail${snail_num}\"\)
#root -l -b -q aggregate.C\(\"snail${snail_num}\",4\)
#root -l -b -q aggregate.C\(\"snail${snail_num}\",0,1\)
#root -l -b -q aggregate.C\(\"snail${snail_num}\",4,1\)

if [ $do_grand -eq 1 ]; then
  root -l -b -q $COMPMON_GRAND/grandConstruction/buildGrandRootfile.C\($run_code\)
fi
root -l -b -q $COMPMON_GRAND/grandOnline/aggregateGrand.C\(${snail_num}\)
if [ $do_expt_plots -eq 1 ]; then
  $COMPMON_GRAND/exptPlots.sh $run_code
fi
#root -l -b -q $COMPMON_GRAND/grandMacros/snailwisePlots.C\($run_code\)
#root -l -b -q $COMPMON_GRAND/grandMacros/runwisePlots.C\($run_code\)
#root -l -b -q $COMPMON_GRAND/grandMacros/cyclewisePlots.C\($run_code\)

date=`date +"%Y-%m-%d"`
time=`date +"%T"`
python $COMPMON_ONLINE/write_html.py $date $time index.html

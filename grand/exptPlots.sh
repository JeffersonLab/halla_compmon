#!/bin/bash

root -l -b -q $COMPMON_GRAND/grandOnline/snailwisePlots.C\($1\)
root -l -b -q $COMPMON_GRAND/grandOnline/runwisePlots.C\($1\)
root -l -b -q $COMPMON_GRAND/grandOnline/cyclewisePlots.C\($1\)

date=`date +"%Y-%m-%d"`
time=`date +"%T"`
python $COMPMON_ONLINE/write_html.py $date $time index.html

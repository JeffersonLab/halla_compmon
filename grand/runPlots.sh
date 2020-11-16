#!/bin/bash

root -l -b -q $COMPMON_LASERCYCLES/laserCycles.C\($1\)
root -l -b -q $COMPMON_GRAND/buildRunRootfile.C\($1\)

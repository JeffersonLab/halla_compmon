#!/bin/bash

for i in {4960..6259}
do
  rm -rf $COMPMON_WEB/runs/Run$i/
  mv $COMPMON_WEB/../../prex2/compton/runs/Run$i/ $COMPMON_WEB/runs/
done

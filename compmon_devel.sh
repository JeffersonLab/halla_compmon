#!/bin/bash
#source /raid1/cornejo/clmeg_repo/jc2_utils.sh
#source /raid1/cornejo/clmeg_repo/repo.sh jl_2.0
#bin/compmon -c compmon_cmu.config $*

configfile=compmon_prex.config
removeconflictingrootfiles=true

if $removeconflictingrootfiles;
then
  ## Get the runnumber, first split the arguments using -r as delimiter
  ## Then to get rid of any other command line options that may come after
  ## the run number, split by space and just get the first one
  runnum=`echo $*| awk -F-r '{print $2}' | awk '{print $1}'`

  if [ "x${runnum}" != "x" ];
  then
    ## If the run number was specified, then make sure we cleanup any old rootfiles
    ## that have an underscore on them. Otherwise ROOT would not overwrite them and
    ## instead append yet another underscore to the ROOT file.
    ## This is guaranteed to confuse us in the future
  
    outPath=`grep outPath ${configfile} | awk -F= '{print $2}'`
    outPrefix=`grep outPrefix ${configfile} | awk -F= '{print $2}'`
    fpat="${outPath}/${outPrefix}${runnum}_[1-9]*.root"
    nofilesfound=true
    echo "Looking for conflicting ROOT files"
    echo "Looking for: ${fpat}"
    for fname in `ls ${fpat} 2> /dev/null`
    do
      echo "Removing found ROOT file: ${fname}"
      rm "${fname}"
      nofilesfound=false
    done
    if ${nofilesfound}
    then
      echo "No conflicting ROOT files found"
    fi
  fi

fi

#bin/compmon -c ${configfile} $*
src/bin_test/compmon -c ${configfile} $*

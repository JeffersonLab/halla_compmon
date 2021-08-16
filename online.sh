#!/bin/bash

PARAMS=""
rootfile_keep=0
web_upload=1
run_num=-1
max_evt=0
replay=0
do_grand=1
do_cycles=1
do_cycleCut=1
DATE="`date +%F`"
TIME="`date +%T`"

# echo "AJZ Says: Do not run this script! It is being edited right now."
# echo "AJZ will inform youm when it is safe to run again"
# exit 0

print_help () {
  echo "Usage: ./online.sh [-r --run <number>] [-h --help] < --nopanguin --rootfile --webupload ... >"
  echo "Required Arguments:"
  echo "    -r, --run: Specify run number to analyze"
  echo "Optional Arguments:"
  echo "    -h, --help: Print this help message"
  echo "    --maxevent <number>: Specify a maximum event number to analyze up to"
  echo "    --nocycles: Don't recalculate laser cycles and bursts"
  echo "    --nogrand: Don't generate any grand or run rootfiles."
  echo "    --nowebupload: Don't generate webpage PDFs"
  echo "    --panguin: Enable panguin window"
  echo "    --replay: Force a compmon re-analysis of the run"
  echo "    --rootfile: Don't delete the plotfile generated at the end of the script"
}

while (( "$#" )); do
  case "$1" in
    -r|--run)
      run_num=$2
      shift 2
      ;;
    --maxevent)
      max_evt=$2
      shift 2
      ;;
    --rootfile)
      rootfile_keep=1
      shift 1
      ;;
    --nowebupload)
      web_upload=0
      shift 1
      ;;
    --replay)
      replay=1
      shift 1
      ;;
    --nogrand)
      do_grand=0
      shift 1
      ;;
    --nocycles)
      do_cycles=0
      shift 1
      ;;
    --nocyclecut)
      do_cycleCut=0
      shift 1
      ;;
    -h|--help)
      print_help
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

# set positional arguments in their proper place
eval set -- "$PARAMS"

if [ $run_num -lt 2239 ]; then
  echo
  echo "Please specify run number!"
  echo
  print_help
  exit 1;
fi

#echo $DATE $TIME
#echo "Upload setting: $web_upload"
#echo "Rootfile setting: $rootfile_keep"
#echo $PARAMS

if [ $web_upload -eq 1 ] && [ ! -d $COMPMON_WEB/runs/Run$run_num ]; then
  mkdir $COMPMON_WEB/runs/Run$run_num
fi

if [ ! -f $COMP_ROOTFILES/compmon_$run_num.root ] || [ $replay -eq 1 ]; then
  ./compmon.sh -r $run_num
fi

if [ ! -f $COMP_ROOTFILES/compmon_$run_num.root ]; then
  echo "Couldn't create or find rootfile. Exiting..."
  exit 1;
fi

if [ $do_cycles == 1 ]; then
  root -l -b -q $COMPMON_LASERCYCLES/laserCycles.C\($run_num\)
  root -l -b -q $COMPMON_LASERCYCLES/bursts.C\($run_num,300\)
fi

if [ $max_evt -gt 0 ]; then
  root -q -b -l $COMPMON_ONLINE/dataQualityCheck.C\($run_num,$max_evt\)
else
  root -q -b -l $COMPMON_ONLINE/dataQualityCheck.C\($run_num\)
fi

if [ $web_upload == 1 ]; then
  root -q -b -l $COMPMON_ONLINE/writeToPDF.C\($run_num\)
  root -l -b -q $COMPMON_ONLINE/plotAllCycles.C\($run_num\)
  root -l -b -q $COMPMON_ONLINE/plotAllBursts.C\($run_num\)
  python $COMPMON_ONLINE/write_html.py $DATE $TIME index.html
fi

if [ $rootfile_keep -eq 0 ]; then
  echo "Deleting rootfile..."
  rm -f $COMPMON_PLOTFILES/compton_online_run_$run_num.root
fi

if [ $do_grand -eq 1 ]; then
  python3 $COMPMON_ONLINE/write_position_file.py $run_num
  $COMPMON_GRAND/runPlots.sh $run_num
fi

if [ $do_cycleCut -eq 1 ]; then
  root -l -b -q $COMPMON_ONLINE/cycleCutOnlinePlots.C\($run_num\)
fi

#rm -f $COMPMON_DIR/core.*
#rm -f $COMPMON_ONLINE/core.*
#rm -f $COMPMON_GRAND/grandConstruction/core.*

echo "Online analysis finished for run $run_num"

exit 0;


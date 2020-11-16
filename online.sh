#!/bin/bash

PARAMS=""
panguin=0
rootfile_keep=0
web_upload=1
run_num=-1
max_evt=0
replay=0
do_grand=1
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
  echo "    --maxevent <number>: Specify a maximum event number to analyze up to"
  echo "    --panguin: Enable panguin window"
  echo "    --nowebupload: Don't generate webpage PDFs"
  echo "    --rootfile: Don't delete the plotfile generated at the end of the script"
  echo "    --replay: Force a compmon re-analysis of the run"
  echo "    -h, --help: Print this help message"
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
    --panguin)
      panguin=1
      shift 1
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

if [ $run_num -lt 0 ]; then
  echo
  echo "Please specify run number!"
  echo
  print_help
  exit 1;
elif [[ $run_num -ge 2239 && $run_num -le 2857 ]]; then
  tests=2
elif [[ $run_num -ge 3114 && $run_num -le 4231 ]]; then
  tests=1
elif [ $run_num -ge 4232 ]; then
  tests=0
else
  echo 
  echo "Please enter real run number."
  echo
  print_help
  exit 1;
fi

#echo $DATE $TIME
#echo "Panguin setting: $panguin"
#echo "Upload setting: $web_upload"
#echo "Rootfile setting: $rootfile_keep"
#echo "Tests setting: $tests"
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

root -l -b -q $COMPMON_LASERCYCLES/laserCycles.C\($run_num\)

if [ $max_evt -gt 0 ]; then
  root -q -b -l $COMPMON_DIR/dataQualityCheck.C\($run_num,$max_evt\)
else
  root -q -b -l $COMPMON_DIR/dataQualityCheck.C\($run_num\)
fi
#if [ $web_upload -eq 1 ]; then
#  root -l -b -q $COMPMON_LASERCYCLES/laserPatternWise.C\($run_num\)
#  cp -f $COMPMON_WEB/runs/Run$run_num/laserCycles_$run_num.dat $COMPMON_MINIRUNS/
#fi

if [ $panguin -eq 1 ]; then
  python $COMPMON_PANGUIN/macros/writeCFG.py $run_num
  $COMPMON_PANGUIN/build/panguin -f $COMPMON_PANGUIN/macros/prex_auto_runs.cfg -r $run_num
fi

if [ $web_upload == 1 ]; then
  root -q -b -l $COMPMON_DIR/writeToPDF.C\($run_num\)
  root -l -b -q $COMPMON_DIR/plotAllCycles.C\($run_num\)
  python $COMPMON_DIR/write_html.py $DATE $TIME index.html
fi

if [ $rootfile_keep -eq 0 ]; then
  echo "Deleting rootfile..."
  rm -f $COMPMON_PLOTFILES/compton_online_run_$run_num.root
fi

if [ $do_grand -eq 1 ]; then
  python3 $COMPMON_DIR/write_position_file.py $run_num
  root -l -q $COMPMON_GRAND/buildRunRootfile.C\($run_num\)
fi

echo "Online analysis finished for run $run_num"

exit 0;


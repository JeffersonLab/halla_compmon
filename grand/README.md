Code to build the compton grand rootfile. Will not work without the regular CompMon included and all environment variables set correctly.

Once downloaded the folder must be pointed to by $COMPMON_GRAND. Also the runPlots/ folder must be created in the main directory.

To run you need to run buildRunRootfile.C on every production run. This only has to be done once. Then you run buildGrandRootfile.C with either argument 1 or 2. (1 == PREX, 2 == CREX) and it will construct the experimental grand rootfile.

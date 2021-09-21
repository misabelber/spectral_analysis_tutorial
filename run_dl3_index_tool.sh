#!/bin/bash
BASEDIR="/fefs/aswg/workspace/maria.bernardos/LSTanalysis/data_MC_comp/data/DL3/real/"

DATES=( 20210904 )
TAILCUT=tailcut84

SUFFIX="_test"
CUT='hardcut'
VERS="v0.7.3"

SCRIPT=/fefs/aswg/workspace/maria.bernardos/GitHub/cta-lstchain/lstchain/tools/lstchain_create_dl3_index_files.py

for date in "${DATES[@]}"
do
    DIR=$BASEDIR/$date/$VERS$SUFFIX/$TAILCUT"_"$CUT/
    #srun -o out.txt python $SCRIPT -d $DIR --overwrite &
    python $SCRIPT -d $DIR --overwrite &
done

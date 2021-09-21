#!/bin/bash
#Script to Run the dl3_tool over a set of DL2 REAL data.
#Change paths according to your system

#Some variables for retrieving the input files
DATES=( 20210904 ) #Dates to analyze
VERS="v0.7.3" #lstchain version
TAILCUT="tailcut84" #tailcut
BASEDIR="/fefs/aswg/data/real/DL2/" #Base directory where to find the dl2 data

#Variables for the output directory
SUFFIX="_test" #Optional suffix
CUT="softcut" #Name of the cut applied to data ( nocuts, softcut, hardcut, hardestcut )

#IRF file
IRF="/fefs/aswg/data/mc/IRF/20200629_prod5_trans_80/zenith_20deg/south_pointing/20210416_v0.7.3_prod5_trans_80_local_taicut_8_4/off0.4deg/irf_20210416_v073_prod5_trans_80_local_taicut_8_4_gamma_point-like_off04deg.fits.gz"

#Configuration file)
CONFIG="/fefs/aswg/workspace/maria.bernardos/spectral_analysis_tutorial/config_dl3_tool_"$CUT.json

#Source to analyze
SRCNAME="Crab"
SRCRA="83.633deg"
SRCDEC="22.01deg"

#Script (in ./cta-lstchain/tools/lstchain_create_dl3_file.py)

SCRIPT="/fefs/aswg/workspace/maria.bernardos/GitHub/cta-lstchain/lstchain/tools/lstchain_create_dl3_file.py"

for date in "${DATES[@]}" #Run over dates
do
    #Define and create an output path for each date
    OUTPUT_PATH="/fefs/aswg/workspace/maria.bernardos/LSTanalysis/data_MC_comp/data/DL3/real/$date/$VERS$SUFFIX/"$TAILCUT"_"$CUT
    `mkdir -p $OUTPUT_PATH`

    FILES=`ls $BASEDIR/$date/$VERS/$TAILCUT/*.h5`
    # Run over the merged run files
    for f in $FILES
    do
        b=$(basename $f)
        run=${b:13:-3}
        if [[ $run != *"."* ]]; then
            FILE=$BASEDIR/$date/$VERS/$TAILCUT/dl2_LST-1.Run$run.h5
            #Execute the script. You can remove the "srun --mem=20g -o out.txt" part if you don't want to use the IT cluster
            #srun --mem=20g -o out.txt python $SCRIPT -d $FILE -o $OUTPUT_PATH --input-irf $IRF --source-name $SRCNAME --source-ra $SRCRA --source-dec $SRCDEC --config $CONFIG --overwrite &
            echo "File $b"
            rm -rf jobs/run_store_all_results_${run}_${CUT}.sh
            echo "#!/bin/bash
#SBATCH -N 3
#SBATCH --mem 100000
ulimit -l unlimited
ulimit -s unlimited
ulimit -a

            python $SCRIPT -d $FILE -o $OUTPUT_PATH --input-irf $IRF --source-name $SRCNAME --source-ra $SRCRA --source-dec $SRCDEC --config $CONFIG --overwrite " >> jobs/run_store_all_results_${run}_${CUT}.sh
            chmod gu+x jobs/run_store_all_results_*
        fi
    done
done

#!/bin/bash
DATASET=$1
OUTDIR=$2
MASS=$3
FILE=$4
#off=5
#num=$(($process+$off))

# CMSSW setup etc
export SCRAM_ARCH="slc5_amd64_gcc462" 
export VO_CMS_SW_DIR="/cms/base/cmssoft"
export COIN_FULL_INDIRECT_RENDERING=1
source $VO_CMS_SW_DIR/cmsset_default.sh


export MYSCRIPT=/users/h2/clseitz/MyCMS/RU/CMSSW_5_3_6/src/RUAnalysis/Ntupler/submit/ntupler_MC_Includes_cfg.py
DIR=/users/h2/clseitz/MyCMS/RU/CMSSW_5_3_6/src/RUAnalysis/Ntupler/
cd $DIR
eval `scramv1 runtime -sh`
#---------------------------------------------------------------

echo "Running script " $MYSCRIPT 

cd $DIR
cmsRun $MYSCRIPT $FILE $DATASET $MASS $OUTDIR > $OUTDIR/$DATASET/Log/Log_${DATASET}_${GLOBALTAG}_${FILE}.log
#stealth
#cmsRun $MYSCRIPT $FILE $MASS $FLAVOR $isr> $OUTDIR/StealthSbottomTo3b_MSb_$MASS/Log/Log_${FILE}.log
#cmsRun $MYSCRIPT $process >& /cms/data23/evah/2011Data/Muons_anaPlots_mFO_${process}.log


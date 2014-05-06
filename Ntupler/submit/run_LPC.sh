#!/bin/bash

cluster=$2
process=$1

# CMSSW setup etc
source /uscmst1/prod/sw/cms/shrc prod
export SCRAM_ARCH="slc5_amd64_gcc462"

export MYSCRIPT=/uscms/home/clseitz/MyCMS/RU/CMSSW_5_2_2/src/RUAnalysis/Ntupler/submit/ntupler_data2012_LPC_cfg.py
OUTDIR=/uscms_data/d3/clseitz/ThreeJet/RUNtuple/SingleMuDataset/BHT650_SHT750
DIR=/uscms/home/clseitz/MyCMS/RU/CMSSW_5_2_2/src/
cd $DIR
eval `scramv1 runtime -sh`
#---------------------------------------------------------------

echo "Running script " $MYSCRIPT 

cp $MYSCRIPT $OUTDIR

cmsRun $MYSCRIPT $process  >& $OUTDIR/Log/MCProd_${process}.log



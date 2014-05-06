#! /bin/env python
import os
import sys
#make all the jdl scripts to submit the data at Rutgers
name_list = os.popen("cat ../python/AllData.list | awk '{print $1}'").readlines()
globaltag_list = os.popen("cat ../python/AllData.list | awk '{print $2}'").readlines()
nfiles_list = os.popen("cat ../python/AllData.list | awk '{print $3}'").readlines()
nDataSet=len(name_list)
k=0
OutDir='/cms/clseitz/ThreeJet/RUNtuple/Jan29ReRun_53xJES/MultiJet53xv2_BaseMu40_SelQuad60Di20_plots'
cmd = 'mkdir -p '+OutDir+'/Log'
os.popen(cmd)
for k in range(0,nDataSet):
    file_list = os.popen("ls -1 /cms/dan/clseitz/ThreeJet/TLBSM53xv2/MultiJet/SoftLinks/"+name_list[k][0:-1]+"*/*root").readlines()
    nFiles = len(file_list)
    ncfg=int(nfiles_list[k][0:-1])
    fileCounter=0
    
    FILE = open("batch_"+name_list[k][0:-1]+"_"+globaltag_list[k][0:-1]+".jdl","w")
    FILE.write('OutDir = '+OutDir+' \n')
    FILE.write('universe = vanilla \n')
    FILE.write('Requirements = (Arch == "X86_64") \n')
    FILE.write('Executable = run_Full2012Data.sh \n')
    FILE.write('should_transfer_files = NO \n')
    
    FILE.write('Dataset = '+name_list[k][0:-1]+'\n')
    FILE.write('Globaltag= '+globaltag_list[k][0:-1]+'\n')
    FILE.write('Output =  $(OutDir)/Log/TLBSM_'+name_list[k][0:-1]+'_'+globaltag_list[k][0:-1]+'_$(Process).stdout \n')
    FILE.write('Error  =  $(OutDir)/Log/TLBSM_'+name_list[k][0:-1]+'_'+globaltag_list[k][0:-1]+'_$(Process).stderr \n')
    FILE.write('Log    =  $(OutDir)/Log/TLBSM_'+name_list[k][0:-1]+'_'+globaltag_list[k][0:-1]+'_$(Process).stdlog \n')
    FILE.write('Arguments = $(Dataset) $(OutDir) $(Globaltag) $(Process) \n')
    FILE.write('Queue '+nfiles_list[k][0:-1]+'\n')
    FILE.close()
    #cmd = 'condor_submit  -p batch_'+name_list[k][0:-1]+'_'+globaltag_list[k][0:-1]+'.jdl'
    #os.popen(cmd)
        

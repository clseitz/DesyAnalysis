#! /bin/env python
import os
import sys

mass = 1250
type = 'Decoup'
flavor = '113_223'

#file_list = os.popen("cat MultiJetData_16.128fb.list").readlines()
#file_list = os.popen("cat QCD_HT500-1000.list").readlines()
#file_list = os.popen("cat QCD_HT250-500.list").readlines()

#file_list = os.popen("cat Stealth"+str(mass)+"_ISR"+str(ISR)+".list").readlines()
#file_list = os.popen("ls -1 /cms/clseitz/ThreeJet/MC/AOD/RPV_M"+str(mass)+"_"+str(flavor)+"_"+str(type)+"/*/*root").readlines()
name_list = os.popen("cat AllData.list | awk '{print $1}'").readlines()
globaltag_list = os.popen("cat AllData.list | awk '{print $2}'").readlines()
nfiles_list = os.popen("cat AllData.list | awk '{print $3}'").readlines()
nDataSet=len(name_list)
k=0
for k in range(0,nDataSet):
    file_list = os.popen("ls -1 /cms/dan/clseitz/ThreeJet/TLBSM53xv2/MultiJet/SoftLinks/"+name_list[k][0:-1]+"*/*root").readlines()
    nFiles = len(file_list)
    ncfg=int(nfiles_list[k][0:-1])
    fileCounter=0
    for indexcount in range(0,ncfg):
        nFilesPercfg=nFiles/ncfg+1
        print nFilesPercfg, nFiles, ncfg
        fileCounter=indexcount*nFilesPercfg
        if(len(file_list) > fileCounter):
            FILE = open(name_list[k][0:-1]+"_"+globaltag_list[k][0:-1]+"_"+str(indexcount)+"_cfi.py","w")
            FILE.write("import FWCore.ParameterSet.Config as cms \n")
            FILE.write('source = cms.Source("PoolSource", \n')
            FILE.write("fileNames = cms.untracked.vstring( \n ")
            
            while fileCounter < (indexcount+1)*nFilesPercfg:
                if(len(file_list) > fileCounter):
                    this_fout = file_list[fileCounter][0:-1]
                    FILE.write("'file:"+this_fout+"', \n")
                fileCounter=fileCounter+1
            FILE.write( "), \n")
        #FILE.write( "duplicateCheckMode = cms.untracked.string('noDuplicateCheck') \n")
            FILE.write( ") \n")
            print fileCounter, ncfg
        
        
            FILE.close()
    k=k+1
        
        

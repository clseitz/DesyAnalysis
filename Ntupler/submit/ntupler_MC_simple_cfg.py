import FWCore.ParameterSet.Config as cms
import string
import os,sys

print sys.argv[2]
index = int (sys.argv[2])

############################################
#######Needs to be adjusted to number of jobs###########
newInd=index % 20
#####################

suffix_list=[]
output_list=[]
dir_list=[]

##################RPV#######################
suffix_list = os.popen('ls -1 /cms/clseitz/ThreeJet/TLBSM/TLBSM53xv2/PAT/RPV/RPV_M750_112/*root').readlines()
#this_fin_nocfi  = "file:"+suffix_list[newInd]

#this_fout0 = outdir+"/"+dataset+"/RPV_M"+mass+"_"+flavor+"/RPV_M"+mass+"_"+flavor+"_"+str(newInd)+"_plots.root"
#this_fout1 = outdir+"/"+dataset+"/RPV_M"+mass+"_"+flavor+"/RPV_M"+mass+"_"+flavor+"_"+str(newInd)+"_tree.root"
#############Stealth################
#suffix_list = os.popen('ls -1 /cms/clseitz-2/MakeDi_Example/PAT/*root').readlines()
this_fin_nocfi  = "file:"+suffix_list[newInd]

#this_fout0 = "/cms/clseitz-2/MakeDi_Example/PAT/tth_"+str(newInd)+"_plots.root"
#this_fout1 = "/cms/clseitz-2/MakeDi_Example/PAT/tth_"+str(newInd)+"_tree.root"
this_fout0 = "th_"+str(newInd)+"_plots.root"
this_fout1 = "tth_"+str(newInd)+"_tree.root"



#################Bprim############################
#suffix_list = os.popen('ls -1 /cms/clseitz/'+dataset+'/*root').readlines()
#this_fin_nocfi  = "file:"+suffix_list[newInd]

#this_fout0 = outdir+"/"+dataset+"/Bprime_M"+mass+"_"+flavor+"_"+str(newInd)+"_plots.root"
#this_fout1 = outdir+"/"+dataset+"/Bprime_M"+mass+"_"+flavor+"_"+str(newInd)+"_tree.root"


print this_fin_nocfi
print this_fout0
print this_fout1

process = cms.Process("data2")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string( 'START53_V15::All' )

# START53_V7E is the old tag used by TLBSM53x
# process.GlobalTag.globaltag = cms.string( 'START53_V7E::All' )
# For data use an AN3 GT like the following
# process.GlobalTag.globaltag = cms.string( 'FT_53_V10_AN3::All' )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# JetCorrectionService ak5PFLchs1FastL2L3 for MC
# JetCorrectionService ak5PFchsL1FastL2L3Residual  for data
# jetcorrserv = 'ak5PFchsL1FastL2L3Residual'
jetcorrserv = 'ak5PFchsL1FastL2L3'

## load jet corrections
# process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
process.load("RUAnalysis.Ntupler.ak5PFchsJetCorrectionServices_cff")
process.prefer(jetcorrserv)


process.maxEvents = cms.untracked.PSet(
#      input = cms.untracked.int32(200)
     input = cms.untracked.int32(-1)
     )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(this_fin_nocfi),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )


##process.source = cms.Source("PoolSource",
##  skipEvents = cms.untracked.uint32(0),
##  fileNames = cms.untracked.vstring(
##    '',
##    ''
##    ),
##  duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
##)

process.data2 = cms.EDAnalyzer('Ntupler',
                               sumPtMin       = cms.untracked.double(0.0),
                               debug          = cms.untracked.bool(False),
                               NtuplePlots = cms.untracked.string(this_fout0),
                               NtupleTree= cms.untracked.string(this_fout1),
                               #PatJetType     = cms.untracked.vstring('goodPatJetsCA8PF'),
                               #PatJetType     = cms.untracked.vstring('goodPatJetsPFlow'),
                               #PatJetType     = cms.untracked.vstring('goodPatJetsCA8PF'),
                               PatJetType = cms.untracked.vstring('goodPatJetsPFlow','goodPatJetsCA8PF','goodPatJetsCA8PrunedPF'),
                               PrimaryVertex = cms.untracked.string('goodOfflinePrimaryVertices'),
                               METtype = cms.untracked.string('patMETsPFlow'),
                               htTrigger      = cms.untracked.string('HLT_HT360_v2'),
                               IsoMu12_v1     = cms.untracked.string('HLT_IsoMu12_v1'),
                               IsoMu17_v6     = cms.untracked.string('HLT_IsoMu17_v6'),
                               eeta           = cms.untracked.double(2.1),
                               ept            = cms.untracked.double(20.0),
                               meta           =  cms.untracked.double(2.1),
                               mpt            = cms.untracked.double(20.0),
                               isData         =  cms.untracked.bool(False),
                               NjetsMin       = cms.untracked.int32(4),
                               NjetsMax       = cms.untracked.int32(1000),
                               jetptcut       = cms.untracked.double(20),
                               etacut         = cms.untracked.double(3.),
                               jecAdj = cms.untracked.string('none'), # Can also be "up" or "down
															 jetCorrectionService = cms.untracked.string(jetcorrserv),
                               #selection trigger
                               #Qua60_Di20
                               #TriggerNamesSel = cms.untracked.vstring('HLT_QuadJet60_DiJet20_v1','HLT_QuadJet60_DiJet20_v2','HLT_QuadJet60_DiJet20_v3','HLT_QuadJet60_DiJet20_v4',
                               #                                     'HLT_QuadJet60_DiJet20_v5','HLT_QuadJet60_DiJet20_v6','HLT_QuadJet60_DiJet20_v7'),
                               #Quad70
                               #TriggerNamesSel = cms.untracked.vstring('HLT_QuadJet80_v1','HLT_QuadJet80_v2','HLT_QuadJet80_v3','HLT_QuadJet80_v4',
                               #'HLT_QuadJet80_v5','HLT_QuadJet80_v6','HLT_QuadJet80_v7'),

                               TriggerNamesSel = cms.untracked.vstring('HLT_SixJet45_v1','HLT_SixJet45_v2','HLT_SixJet45_v3','HLT_SixJet45_v4'),
                               
                               #base trigger
                               #TriggerNamesBase = cms.untracked.vstring('HLT_QuadJet50_v9','HLT_QuadJet50_v8'),
                               #,'HLT_QuadJet50_v3','HLT_QuadJet50_v4',
                               #                                      'HLT_QuadJet50_v5','HLT_QuadJet50_v6','HLT_QuadJet50_v7'),
                               #TriggerNamesBase = cms.untracked.vstring('HLT_QuadJet70_v1','HLT_QuadJet70_v2','HLT_QuadJet70_v3','HLT_QuadJet70_v4',
                               #                                   'HLT_QuadJet70_v5','HLT_QuadJet70_v6','HLT_QuadJet70_v7'),
                               TriggerNamesBase = cms.untracked.vstring('HLT_SixJet35_v1','HLT_SixJet35_v2','HLT_SixJet35_v3','HLT_SixJet35_v4'),
                               #JSON file 900/pb May18
#                               JSONFilename = cms.untracked.string('/home/clseitz/MyCMS/RU/CMSSW_5_2_2/src/RUAnalysis/Ntupler/submit/Cert_190456-194076_8TeV_PromptReco_Collisions12_JSON.txt'),
                               #JSON File 2.42/fb Jun01
                               JSONFilename = cms.untracked.string('/home/clseitz/MyCMS/RU/CMSSW_5_2_2/src/RUAnalysis/Ntupler/submit/Cert_190456-195016_8TeV_PromptReco_Collisions12_JSON.txt'),
                               )


process.p = cms.Path(process.data2)




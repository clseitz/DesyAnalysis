import FWCore.ParameterSet.Config as cms
import string
import os,sys

#this_fin_nocfi  = "file:../../../PhysicsTools/PatAlgos/test/miniAOD/patTuple_mini.root"


#this_fin_nocfi  = "file:/nfs/dust/cms/user/clseitz/MiniAOD/Files/TT_Tune4C_13TeV-pythia8-tauola_PAT.root"
this_fin_nocfi  = "file:/nfs/dust/cms/user/clseitz/MiniAOD/Files/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT.root"
#this_fin_nocfi  = "file:/nfs/dust/cms/user/clseitz/MiniAOD/Files/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PAT.root"
#this_fin_nocfi  = "file:/nfs/dust/cms/user/clseitz/MiniAOD/CMSSW_7_0_4/src/PhysicsTools/PatAlgos/test/miniAOD/miniAOD-prod_PAT.root"

this_fout0 = "test_plots_t2tt.root"
this_fout1 = "test_tree_t2tt.root"



process = cms.Process("data2")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string( 'START53_V15::All' )

# START53_V7E is the old tag used by TLBSM53x
# process.GlobalTag.globaltag = cms.string( 'START53_V7E::All' )
# For data use an AN3 GT like the following
# process.GlobalTag.globaltag = cms.string( 'FT_53_V10_AN3::All' )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200

# JetCorrectionService ak5PFLchs1FastL2L3 for MC
# JetCorrectionService ak5PFchsL1FastL2L3Residual  for data
# jetcorrserv = 'ak5PFchsL1FastL2L3Residual'

jetcorrserv = 'ak5PFchsL1FastL2L3'


## load jet corrections
# process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
process.load("DesyAnalysis.Ntupler.ak5PFchsJetCorrectionServices_cff")
process.prefer(jetcorrserv)


process.maxEvents = cms.untracked.PSet(
#      input = cms.untracked.int32(200)
     input = cms.untracked.int32(-1)
     )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(this_fin_nocfi),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

inputJetCorrLabelAK7PFchs = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJetsPFlow = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag('packedPFCandidates'),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax = cms.double(6.0),
    Ghost_EtaMax = cms.double(7.0)
    )
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
#addJetCollection(process,
#                 labelName = 'slimmedJetsCA8',
#                 jetSource = cms.InputTag('ca8PFJetsPFlow'),
#                 trackSource = cms.InputTag('unpackedTracksAndVertices'),
#                 pvSource = cms.InputTag('unpackedTracksAndVertices'),
#                 jetCorrections = inputJetCorrLabelAK7PFchs
#                 )


from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True) # no idea while doArea is false by default, but it's True in RECO so we have to set it
process.ak5PFJetsCHS = ak5PFJets.clone(src = 'chs', doAreaFastjet = True) # no idea while doArea is false by default, but it's True in RECO so we have to set it
process.ak5GenJets = ak5GenJets.clone(src = 'packedGenParticles')


from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(
   process,
   postfix   = "",
   labelName = 'AK5PFCHS',
   jetSource = cms.InputTag('ak5PFJetsCHS'),
   trackSource = cms.InputTag('unpackedTracksAndVertices'),
   pvSource = cms.InputTag('unpackedTracksAndVertices'),
   jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
   btagDiscriminators = [      'combinedSecondaryVertexBJetTags'     ]
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
                               debug          = cms.untracked.bool(True),
                               NtuplePlots = cms.untracked.string(this_fout0),
                               NtupleTree= cms.untracked.string(this_fout1),
                               #PatJetType     = cms.untracked.vstring('goodPatJetsCA8PF'),
                               #PatJetType     = cms.untracked.vstring('goodPatJetsPFlow'),
                               #PatJetType     = cms.untracked.vstring('goodPatJetsCA8PF'),
                               PatJetType = cms.untracked.vstring('slimmedJets','AK5PFCHS','slimmedJetsCA8'),
                               PrimaryVertex = cms.untracked.string('offlineSlimmedPrimaryVertices'),
                               METtype = cms.untracked.string('slimmedMETs'),
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
                               rhoIsoInputTag          = cms.InputTag("fixedGridRhoFastjetAll"),
                               isoValInputTags = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                                               cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                                               cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
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


process.p = cms.Sequence(process.data2)


#process.endpath= cms.EndPath(process.OUT)
process.p0 = cms.Path( process.p)

################################################################################
#
# ak5PFchsJetCorrectionServices_cff.py
# -------------------------
#
# Jet corrections for ak5PFchs jets. CMSSW seems not to include these
# JECs, even though others like ak5PF are provided.
################################################################################

from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

ak5PFchsL2Relative = cms.ESProducer(
    'LXXXCorrectionESProducer',
    level     = cms.string('L2Relative'),
    algorithm = cms.string('AK5PFchs')
    )
ak5PFchsL3Absolute = cms.ESProducer(
    'LXXXCorrectionESProducer',
    level     = cms.string('L3Absolute'),
    algorithm = cms.string('AK5PFchs')
    )
ak5PFchsL2L3 = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak5PFchsL2Relative','ak5PFchsL3Absolute')
    )
ak5PFchsL1Fastjet = cms.ESProducer(
    'L1FastjetCorrectionESProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK5PFchs'),
    srcRho      = cms.InputTag('kt6PFJets','rho')
    )
ak5PFchsL1FastL2L3 = ak5PFchsL2L3.clone()
ak5PFchsL1FastL2L3.correctors.insert(0,'ak5PFchsL1Fastjet')

ak5PFchsResidual  = ak5CaloResidual.clone( algorithm = 'AK5PFchs' )

ak5PFchsL1FastL2L3Residual = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak5PFchsL1Fastjet','ak5PFchsL2Relative','ak5PFchsL3Absolute','ak5PFchsResidual')
    )

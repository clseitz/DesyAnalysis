 // -*- C++ -*-
//
// Package:    Ntupler
// Class:      Ntupler
// 
/**\class Ntupler Ntupler.cc DesyAnalysis/Ntupler/src/Ntupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Claudia Seitz
//         Created:  Mon Apr  9 12:14:40 EDT 2012
// $Id: Ntupler.cc,v 1.28 2013/08/08 08:26:13 clseitz Exp $
// Modified for use at Desy May 2014
//


// system include files
#include <memory>


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h" // based on DataFormats/Candidate/interface/Particle.h
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"  
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"


#include "TTree.h"

//NEW
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"



#include "DataFormats/PatCandidates/interface/Electron.h" 
#include "DataFormats/PatCandidates/interface/Muon.h" 
#include "DataFormats/PatCandidates/interface/MET.h" 
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
//filtering susy events
//#include "UserCode/ModelFilter/interface/ModelFilter.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
//pile upt stuff
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//Own libraries
#include "DesyAnalysis/Ntupler/interface/Ntupler.h"
#include "DesyAnalysis/Ntupler/interface/NtpReader.h"
#include <memory>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>

using namespace std;



// constants, enums and typedefs
//

typedef struct {
	reco::Candidate::LorentzVector adjJet, diffVec;
	const pat::Jet *origJet;
	double jecUnc;
} jetElem;

// Comparison function for jet list
// Want highest pt first
bool cmpJets(jetElem first, jetElem second)
{
	return (first.adjJet.Pt() > second.adjJet.Pt());
}


//
// static data member definitions
//

//
// constructors and destructor
//
Ntupler::Ntupler(const edm::ParameterSet& iConfig)

{
  _sumPtMin        = iConfig.getUntrackedParameter<double>("sumPtMin",         300.0);
  _debug           = iConfig.getUntrackedParameter<bool>  ("debug",            false);
  _ntuplePlots  = iConfig.getUntrackedParameter<string>("NtuplePlots",  "PatJeTs_test.root");
  _ntupleTree = iConfig.getUntrackedParameter<string>("NtupleTree", "PatJets_testTree.root");
  //  _patJetType      = iConfig.getUntrackedParameter<string>("PatJetType",      "selectedPatJets");
  _patJetType      = iConfig.getUntrackedParameter<std::vector<std::string> >("PatJetType", std::vector<std::string> ());
  _primaryVertex   = iConfig.getUntrackedParameter<string> ("PrimaryVertex","goodOfflinePrimaryVertices");
  _jecAdj   			 = iConfig.getUntrackedParameter<string> ("jecAdj", "none");
  _jetCorrectionService = iConfig.getUntrackedParameter<string> ("jetCorrectionService",
		"ak5PFL1L2L3");
  _METtype   = iConfig.getUntrackedParameter<string> ("METtype","patMETsPFlow");
  _njetsMin        = iConfig.getUntrackedParameter<int>   ("NjetsMin",         4);
  _njetsMax        = iConfig.getUntrackedParameter<int>   ("NjetsMax",         4);
  _etacut          = iConfig.getUntrackedParameter<double>("etacut",           3.0); 
  _jetptcut        = iConfig.getUntrackedParameter<double>("jetptcut",         20.0);

  _eeta            = iConfig.getUntrackedParameter<double>("eeta",           2.1); 
  _ept             = iConfig.getUntrackedParameter<double>("ept",         20.0);
  
  _meta            = iConfig.getUntrackedParameter<double>("meta",           2.1); 
  _mpt             = iConfig.getUntrackedParameter<double>("mpt",         20.0);
  
  _pheta           = iConfig.getUntrackedParameter<double>("pheta",           1.45); 
  _phpt            = iConfig.getUntrackedParameter<double>("phpt",         30.0);
  
  _nbTagsMin       = iConfig.getUntrackedParameter<int>   ("nbTagsMin",        0);
  _nbTagsMax       = iConfig.getUntrackedParameter<int>   ("nbTagsMax",        1000);
  _isData          = iConfig.getUntrackedParameter<bool>  ("isData",           true);
  _isSUSY          = iConfig.getUntrackedParameter<bool>  ("isSUSY",           false);
  _noTripletBtag   = iConfig.getUntrackedParameter<bool>  ("noTripletBtag",    false);
  _nTripletBtagsMin= iConfig.getUntrackedParameter<int>   ("nTripletBtagsMin", 0);
  _nTripletBtagsMax= iConfig.getUntrackedParameter<int>   ("nTripletBtagsMax", 1000);
  _doBtagEff       = iConfig.getUntrackedParameter<bool>   ("doBtagEff", true);
  _rhoIsoInputTag  = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
   _isoValInputTags        = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags"); 
  fTriggerNamesSel = iConfig.getUntrackedParameter<std::vector<std::string> >("TriggerNamesSel", std::vector<std::string>());
  for (std::vector<std::string>::iterator It = fTriggerNamesSel.begin(); It != fTriggerNamesSel.end(); ++It) {
    fTriggerMap[*It] = false;
  }
  fTriggerNamesBase = iConfig.getUntrackedParameter<std::vector<std::string> >("TriggerNamesBase", std::vector<std::string>());
  for (std::vector<std::string>::iterator It = fTriggerNamesBase.begin(); It != fTriggerNamesBase.end(); ++It) {
    fTriggerMap2[*It] = false;
  }
  fTriggerNamesBase2 = iConfig.getUntrackedParameter<std::vector<std::string> >("TriggerNamesBase2", std::vector<std::string>());
  for (std::vector<std::string>::iterator It = fTriggerNamesBase2.begin(); It != fTriggerNamesBase2.end(); ++It) {
    fTriggerMapBase2[*It] = false;
  }
  JSONFilename  = iConfig.getUntrackedParameter<string>("JSONFilename","Cert_160404-166502_7TeV_PromptReco_Collisions11_JSON.txt");


}


Ntupler::~Ntupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


static double getJERAdj(double recoPt, const pat::Jet &jet, bool down)
{
	const reco::GenJet *genJetPtr = jet.genJet();
	if (genJetPtr == NULL) {
		// cout << "Bad GenJet pointer\n";
		return (1.0);
	}
	double ptdiff = recoPt - genJetPtr->pt();
	if (down)
		ptdiff *= -1.0;
	double ptscale = ((ptdiff * 0.1) + recoPt) / recoPt;
	// cout << " JER scaling " << ptscale << " ";
	if (ptscale < 0.0)
		ptscale = 1.0;
	return (ptscale);
}


// ------------ method called for each event  ------------
void
Ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   /////////////////
   //GET EVT INFO
   ///////////////////////////////////////////////////////////
   run   = iEvent.id().run();
   event = iEvent.id().event();
   lumis = iEvent.id().luminosityBlock();
   /////////////////////////////////////////////////////////////
   ////CHECK JSON
   /////////////////////////////////////////////////////////////
   GoodRun=kFALSE;
   UseJson(GoodRuns,GoodLumiStart,GoodLumiEnd,nGoodRuns,run,lumis);
   ////////////////////////////////////////////////////////////
   //only use good runs from the JSON file
   if(GoodRun){
     //////////////////////
     ///// check triggers
     ////////////////////////////////////////////////////////// 
     // Look for one of the triggers we care about

     bool HasTrigger = false;
     bool HasTrigger2 = false;
     bool HasTriggerBase2 = false;
     if (_isData) {
       getTriggerDecision(iEvent, fTriggerMap);
       for (std::map<std::string, bool>::iterator It = fTriggerMap.begin(); It != fTriggerMap.end(); ++It) {
	 if (It->second) {
	   HasTrigger = true;
	   break;
	 }
       }//triggers

       getTriggerDecision(iEvent, fTriggerMap2);
       for (std::map<std::string, bool>::iterator It = fTriggerMap2.begin(); It != fTriggerMap2.end(); ++It) {
         if (It->second) {
           HasTrigger2 = true;
           break;
         }
       }//triggers
       getTriggerDecision(iEvent, fTriggerMapBase2);
       for (std::map<std::string, bool>::iterator It = fTriggerMapBase2.begin(); It != fTriggerMapBase2.end(); ++It) {
         if (It->second) {
           HasTriggerBase2 = true;
           break;
         }
       }//triggers

     }//isData
     else {
       HasTrigger= true;
       HasTrigger2= true;
       HasTriggerBase2 = true;
       if (_isSUSY){
	 	 GetSUSYpoint(iEvent,iSetup);

       }//isSUSY
     }//else isData
    /////////////////////////////////////////////////////
    ////CLEAN UP VARIABLES
    ////////////////////////////////////////////////////
     HasSelTrigger = HasTrigger;
     DataIs=_isData;
     HasBaseTrigger = HasTrigger2 || HasTriggerBase2;
     if (HasTrigger || HasTrigger2 || HasTriggerBase2){
       fGoodJets.clear(); fCleanJets.clear(); 
       nGoodJets=0; nCleanJets=0;

       fGoodPFJets.clear(); fCleanPFJets.clear(); 
       nGoodPFJets=0; nCleanPFJets=0;

       fGoodCA8PFJets.clear(); fCleanCA8PFJets.clear(); 
       nGoodCA8PFJets=0; nCleanCA8PFJets=0;

       fGoodElectrons.clear(); fCleanElectrons.clear();
       nGoodElectrons=0; nCleanElectrons=0; 
       fCleanElectronsPFrelIso.clear(); 
       fCleanElectronsPFabsIso.clear(); 
       fGoodMuons.clear(); fCleanMuons.clear(); 
       fCleanMuonsPFrelIso.clear(); 
       fCleanMuonsPFabsIso.clear(); 
       nGoodMuons=0; nCleanMuons=0;  
       fGoodPhotons.clear(); fCleanPhotons.clear();
       nGoodPhotons=0; nCleanPhotons=0;
       fGoodVtx.clear();
       nGoodVtx=0;
       nElectrons=0;
       nMuons=0;
       nPhotons=0;

       nGoodPFJets=0;       nCleanPFJets=0;
       nGoodCA8PFJets=0;       nCleanCA8PFJets=0;
       nGoodCA8PrunedPFJets=0;       nCleanCA8PrunedPFJets=0;
       nCA8PrunedPFJets=0;
       nCA8PFJets=0;
       nPFJets=0;
       /*
       for (int i=0; i<200; ++i)
	 {
	   pdgID[i] = -99;
	   MCpx[i] = -99;
	   MCpy[i] = -99;
	   MCpz[i] = -99;
	   MCe[i] = -99;
	 }
       */
       Triplet.clear();   
       sumScalarPtTriplet.clear();
       sumVectorPtTriplet.clear();
       massTriplet.clear();
       nTriplets=0; q=0; //basically just triplet counting
       IsVtxGood = 0; 
       ///////////////////////////////////////////////////////
       /////DO OBJECT ID       ///////////////////////////////////////////////////////
       //Select all the objects int the event (vertex function makes some plots)
       
       
       //JETS already have loose JetID applied
       //	 edm::Handle<edm::View<pat::Jet> > fGoodPFJets;
       // edm::Handle< std::vector<pat::Jet> > fGoodPFJets;
       //       edm::Handle< std::vector<pat::Jet> > fCleanPFJets;
       //iEvent.getByLabel(_patJetType[0], fCleanPFJets);
       // edm::Handle< std::vector<pat::Jet> > fGoodCA8PFJets;
       //edm::Handle< std::vector<pat::Jet> > fCleanCA8PFJets;
       //iEvent.getByLabel(_patJetType[1], fCleanCA8PFJets);
       //need as PAT to get btag stuff
       //   edm::Handle< std::vector<pat::Jet> > fGoodCA8PrunedPFJets;
       //edm::Handle< std::vector<pat::Jet> > fCleanCA8PrunedPFJets;
       //iEvent.getByLabel(_patJetType[2], fCleanCA8PrunedPFJets);
       
       
       ////////////////////
       if(!_isData) GetTruePileUp(iEvent);
       
       DoVertexID(iEvent);
       if(_debug) cout<<"before jetid"<<endl;
       DoJetID(iEvent, iSetup, _patJetType[0], 0);
       DoJetID(iEvent, iSetup, _patJetType[1], 1);
       if(_debug)cout<<"after jetid"<<endl;
       DoElectronID(iEvent);
       if(_debug)cout<<"after eleid"<<endl;
       DoMuonID(iEvent);
       if(_debug)cout<<"after muid"<<endl;
       //Photons are not implemented at the moment
       //DoPhotonID(iEvent);
       DoMETID(iEvent);
       if(_debug)  cout<<"after metid"<<endl;
       ///////REMOVE OVERLAP IN OCJECT COLLECTION
       //make some plots before cleanup 
       h_nGoodPFJets->Fill(nGoodPFJets);       
       h_nGoodJets->Fill(nGoodJets);
       h_nGoodCA8PFJets->Fill(nGoodCA8PFJets);
       h_nGoodCA8PrunedPFJets->Fill(nGoodCA8PrunedPFJets);

       h_nGoodElectrons->Fill(nGoodElectrons);  
       h_nGoodMuons->Fill(nGoodMuons);
       h_nGoodPhotons->Fill(nGoodPhotons);

       ///////////////////////////////////////////////////
       //cross clean the collections and also get isolation for e and mu
       DoCleanUp(iEvent,fGoodMuons,fGoodElectrons,fGoodPhotons,fGoodJets);
       if(_debug) cout<<"after cleanup"<<endl;
       /*fCleanMuons = fGoodMuons;
       fCleanElectrons = fGoodElectrons;
       fCleanPhotons= fGoodPhotons;
       nCleanMuons = nGoodMuons; 
       nCleanElectrons = nGoodElectrons;
       nCleanPhotons= nGoodPhotons;
       */

       ///////////////////////////////////////////////
       //make plots after the clean up
       
       h_nCleanElectrons->Fill(nCleanElectrons);
       h_nCleanMuons->Fill(nCleanMuons);
       h_nCleanPhotons->Fill(nCleanPhotons);
       h_nCleanJets->Fill(nCleanJets);
       h_nCleanPFJets->Fill(nCleanPFJets);
       h_nCleanCA8PFJets->Fill(nCleanCA8PFJets);
       h_nCleanCA8PrunedPFJets->Fill(nCleanCA8PrunedPFJets);
       /////////////////
       //////KINEMATIC PLOTS OF OBJECTS + STUFF FOR THE TREE
       ////////////////
       //make some kinematic plots and write out variables for the tree
       
       nPFJets=nCleanPFJets;
       //       std::auto_ptr<reco:GenParticleCollection> &parents(new reco::GenParticleCollection());
       


       int i=0;

       for (std::vector<pat::Jet>::const_iterator Jet = fCleanPFJets.begin(); Jet != fCleanPFJets.end(); ++Jet) { 
	 
	 h_chargedHadronEnergyFraction->Fill(Jet->correctedJet("Uncorrected").chargedHadronEnergyFraction());
	 h_neutralHadronEnergyFraction->Fill(Jet->correctedJet("Uncorrected").neutralHadronEnergyFraction());
	 h_neutralEmEnergyFraction->Fill(Jet->correctedJet("Uncorrected").neutralEmEnergyFraction());
	 h_chargedEmEnergyFraction->Fill(Jet->correctedJet("Uncorrected").chargedEmEnergyFraction());
	 h_numberOfDaughters->Fill(Jet->correctedJet("Uncorrected").numberOfDaughters());
	 h_chargedMultiplicity->Fill(Jet->correctedJet("Uncorrected").chargedMultiplicity());
	 if(i<6){
	 v_PFjet_pt[i]->Fill(Jet->pt()); 
	 v_PFjet_eta[i]->Fill(Jet->eta()); 
	 v_PFjet_phi[i]->Fill(Jet->phi()); 
	 v_PFjet_m[i]->Fill(Jet->mass()); 
	 }
	 /*if(counter60 >=6){
	   if(i<4 && adjJet->pt()>80.0 && fabs(adjJet->eta())<2.5){
	   v_PFjet_pt[i]->Fill(adjJet->pt()); 
	   v_PFjet_eta[i]->Fill(adjJet->eta()); 
	   v_PFjet_phi[i]->Fill(adjJet->phi()); 
	   v_PFjet_m[i]->Fill(adjJet->mass()); 
	   h_chargedHadronEnergyFraction->Fill(Jet->correctedJet("Uncorrected").chargedHadronEnergyFraction());
	   h_neutralHadronEnergyFraction->Fill(Jet->correctedJet("Uncorrected").neutralHadronEnergyFraction());
	   h_neutralEmEnergyFraction->Fill(Jet->correctedJet("Uncorrected").neutralEmEnergyFraction());
	   h_chargedEmEnergyFraction->Fill(Jet->correctedJet("Uncorrected").chargedEmEnergyFraction());
	   h_numberOfDaughters->Fill(Jet->correctedJet("Uncorrected").numberOfDaughters());
	   h_chargedMultiplicity->Fill(Jet->correctedJet("Uncorrected").chargedMultiplicity());


	   }
   if(i>=4 &&i<6 && adjJet->pt()>60.0  && fabs(adjJet->eta())<2.5){
	   v_PFjet_pt[i]->Fill(adjJet->pt()); 
	   v_PFjet_eta[i]->Fill(adjJet->eta()); 
	   v_PFjet_phi[i]->Fill(adjJet->phi()); 
	   v_PFjet_m[i]->Fill(adjJet->mass()); 
	   h_chargedHadronEnergyFraction->Fill(Jet->correctedJet("Uncorrected").chargedHadronEnergyFraction());
	   h_neutralHadronEnergyFraction->Fill(Jet->correctedJet("Uncorrected").neutralHadronEnergyFraction());
	   h_neutralEmEnergyFraction->Fill(Jet->correctedJet("Uncorrected").neutralEmEnergyFraction());
	   h_chargedEmEnergyFraction->Fill(Jet->correctedJet("Uncorrected").chargedEmEnergyFraction());
	   h_numberOfDaughters->Fill(Jet->correctedJet("Uncorrected").numberOfDaughters());
	   h_chargedMultiplicity->Fill(Jet->correctedJet("Uncorrected").chargedMultiplicity());
	 }

	 }*/


	 jet_PF_pt[i]=Jet->pt();

	 jet_PF_px[i]=Jet->px();
	 jet_PF_py[i]=Jet->py();
	 jet_PF_pz[i]=Jet->pz();
	 jet_PF_e[i]=Jet->energy();
	 jet_PF_mass[i]=Jet->mass();
	 jet_PF_et[i]=Jet->et();
	 jet_PF_jec_unc[i]=1;

	 jet_PF_eta[i]=Jet->eta();
	 jet_PF_phi[i]=Jet->phi();
	 jet_PF_theta[i]=Jet->theta();
	 jet_PF_NeutralHad[i]=Jet->correctedJet("Uncorrected").neutralHadronEnergyFraction();
	 jet_PF_area[i]=Jet->jetArea();
	 jet_PF_nJetDaughters[i]=Jet->numberOfDaughters();
	 jet_PF_PartonFlav[i]=Jet->partonFlavour();
     	 bdiscTCHE_PF[i]=Jet->bDiscriminator("trackCountingHighEffBJetTags");
	 bdiscTCHP_PF[i]=Jet->bDiscriminator("trackCountingHighPurBJetTags");
	 bdiscSSVHE_PF[i]=Jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	 bdiscSSSVHP_PF[i]=Jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	 bdiscCSV_PF[i]=Jet->bDiscriminator("combinedSecondaryVertexBJetTags");
	 bdiscJP_PF[i]=Jet->bDiscriminator("jetProbabilityBJetTags");
	 /*if (!_isData) {
	   if(1!=1){
	     int jetMom = -1; 
	     const reco::GenParticle * part = Jet->genParton();
	     if (part){
	       //cout<<"GenParton: "<<part->pdgId()<<" mother:"<<(*part).mother()->pdgId()<<endl;
	       const reco::GenParticle * mom;
	       //cout<<"Aftere "<<part->mother()<<endl;
	       if ((*part).mother()->pdgId() == 3000002){
		 const reco::GenParticle * teenmom = (const reco::GenParticle*) (*part).mother();
		 const reco::GenParticle * grammy = (const reco::GenParticle*) (*teenmom).mother();
		 mom = (const reco::GenParticle*) (*grammy).mother();
	       }
	       else{mom = (const reco::GenParticle*) (*part).mother();}
	       //cout<<part<<endl;
	       //cout<<mom->mass()<<endl;
	       for (int y = 0; y < int(parents->size()); ++y)
		 if (fabs((*parents)[y].p() - mom->p()) < 0.001) jetMom = y;
	       if (jetMom == -1){
		 jetMom = int(parents->size());
		 reco::GenParticle cand(*mom);
		 parents->push_back(cand);
		 std::cout << "Found Mom with number of daughters: " << parents->size() << std::endl;

	       }
	       
	     }
	     jet_PF_JetMom[i]=jetMom;	     
	     //	     cout<<"jetmomL: "<<jetMom<<endl;
	   }
	 }*/
	 i++;
	 
       }

  nCA8PFJets=nCleanCA8PFJets;

  i=0;
  for (std::vector<pat::Jet>::const_iterator Jet = fCleanCA8PFJets.begin(); Jet != fCleanCA8PFJets.end(); ++Jet) {       
    
    if(i<6){
      v_CA8PFjet_pt[i]->Fill(Jet->pt()); 
      v_CA8PFjet_eta[i]->Fill(Jet->eta()); 
      v_CA8PFjet_phi[i]->Fill(Jet->phi()); 
      v_CA8PFjet_m[i]->Fill(Jet->mass()); 
    }
    jet_CA8PF_pt[i]=Jet->pt();
    jet_CA8PF_et[i]=Jet->et();
    jet_CA8PF_eta[i]=Jet->eta();
    jet_CA8PF_phi[i]=Jet->phi();
    jet_CA8PF_theta[i]=Jet->theta();
    jet_CA8PF_mass[i]=Jet->mass();
    jet_CA8PF_area[i]=Jet->jetArea();
    jet_CA8PF_nJetDaughters[i]=Jet->numberOfDaughters();
    
    jet_CA8PF_px[i]=Jet->px();
    jet_CA8PF_py[i]=Jet->py();
    jet_CA8PF_pz[i]=Jet->pz();
    jet_CA8PF_e[i]=Jet->energy(); 
    bdiscTCHE_CA8PF[i]=Jet->bDiscriminator("trackCountingHighEffBJetTags");
    bdiscTCHP_CA8PF[i]=Jet->bDiscriminator("trackCountingHighPurBJetTags");
    bdiscSSVHE_CA8PF[i]=Jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
    bdiscSSSVHP_CA8PF[i]=Jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
    bdiscCSV_CA8PF[i]=Jet->bDiscriminator("combinedSecondaryVertexBJetTags");
    bdiscJP_CA8PF[i]=Jet->bDiscriminator("jetProbabilityBJetTags");
    //	 cout<<Jet->numberOfDaughters()<<endl;
    i++;
  }
  /*
       nCA8PrunedPFJets=nCleanCA8PrunedPFJets;
       i=0;
       for (std::vector<pat::Jet>::const_iterator Jet = fCleanCA8PrunedPFJets->begin(); Jet != fCleanCA8PrunedPFJets->end(); ++Jet) {       
	 if(i<6 && i < nCA8PrunedPFJets){
	   v_CA8PrunedPFjet_pt[i]->Fill(Jet->pt()); 
	   v_CA8PrunedPFjet_eta[i]->Fill(Jet->eta()); 
	   v_CA8PrunedPFjet_phi[i]->Fill(Jet->phi()); 
	   v_CA8PrunedPFjet_m[i]->Fill(Jet->mass()); 
	 }
	 jet_CA8PrunedPF_pt[i]=Jet->pt();
	 jet_CA8PrunedPF_px[i]=Jet->px();
	 jet_CA8PrunedPF_py[i]=Jet->py();
	 jet_CA8PrunedPF_pz[i]=Jet->pz();
	 jet_CA8PrunedPF_e[i]=Jet->energy(); 

	 jet_CA8PrunedPF_et[i]=Jet->et();
	 jet_CA8PrunedPF_eta[i]=Jet->eta();
	 jet_CA8PrunedPF_phi[i]=Jet->phi();
	 jet_CA8PrunedPF_theta[i]=Jet->theta();
	 jet_CA8PrunedPF_mass[i]=Jet->mass();
	 jet_CA8PrunedPF_area[i]=Jet->jetArea();
	 	 jet_CA8PrunedPF_nJetDaughters[i]=Jet->numberOfDaughters();
	 //	 cout<<Jet->numberOfDaughters()<<endl;
	 //if(Jet->numberOfDaughters() >=3) cout<<"FOUND ONE"<<endl;
	 bdiscTCHE_CA8PrunedPF[i]=Jet->bDiscriminator("trackCountingHighEffBJetTags");
	 bdiscTCHP_CA8PrunedPF[i]=Jet->bDiscriminator("trackCountingHighPurBJetTags");
	 bdiscSSVHE_CA8PrunedPF[i]=Jet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	 bdiscSSSVHP_CA8PrunedPF[i]=Jet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	 bdiscCSV_CA8PrunedPF[i]=Jet->bDiscriminator("combinedSecondaryVertexBJetTags");
	 bdiscJP_CA8PrunedPF[i]=Jet->bDiscriminator("jetProbabilityBJetTags");

	 if(Jet->numberOfDaughters() >=2){
	 reco::Jet const * subjet0 = dynamic_cast<reco::Jet const *>(Jet->daughter(0));
	 reco::Jet const * subjet1 = dynamic_cast<reco::Jet const *>(Jet->daughter(1));

	 jet_CA8PrunedPF_subJet1_mass[i]=subjet0->mass();
	 jet_CA8PrunedPF_subJet2_mass[i]=subjet1->mass();
	 jet_CA8PrunedPF_subJet3_mass[i]=subjet1->mass();

	 jet_CA8PrunedPF_subJet1_px[i]=subjet0->px();
	 jet_CA8PrunedPF_subJet2_px[i]=subjet1->px();
	 jet_CA8PrunedPF_subJet3_px[i]=subjet1->px();

	 jet_CA8PrunedPF_subJet1_py[i]=subjet0->py();
	 jet_CA8PrunedPF_subJet2_py[i]=subjet1->py();
	 jet_CA8PrunedPF_subJet3_py[i]=subjet1->py();

	 jet_CA8PrunedPF_subJet1_pz[i]=subjet0->pz();
	 jet_CA8PrunedPF_subJet2_pz[i]=subjet1->pz();
	 jet_CA8PrunedPF_subJet3_pz[i]=subjet1->pz();

	 jet_CA8PrunedPF_subJet1_e[i]=subjet0->energy();
	 jet_CA8PrunedPF_subJet2_e[i]=subjet1->energy();
	 jet_CA8PrunedPF_subJet3_e[i]=subjet1->energy();


	   
	 cout<<"Allthe masses "<<jet_CA8PrunedPF_nJetDaughters[i]<<" "<<subjet0->mass()<<endl;   
	 }
	 	 
	     for(int k =0; k<jet_CA8PrunedPF_nJetDaughters[i];k++){
	       reco::Jet const * subjetk = dynamic_cast<reco::Jet const *>(Jet->daughter(k));
	     cout<<subjetk->mass()<<endl;   
	     }
	 
	 i++;
	 }
       */
       nElectrons=nCleanElectrons;
       for(int i=0; i<nCleanElectrons; i++){
	 if(i<6){
	   v_e_pt[i]->Fill(fCleanElectrons[i].pt()); 
	   v_e_eta[i]->Fill(fCleanElectrons[i].eta()); 
	   v_e_phi[i]->Fill(fCleanElectrons[i].phi()); 
	 }
	 ept[i]=fCleanElectrons[i].pt();
	 epx[i]=fCleanElectrons[i].px();
	 epy[i]=fCleanElectrons[i].py();
	 epz[i]=fCleanElectrons[i].pz();
	 ee[i]=fCleanElectrons[i].energy(); 
	 echarge[i]=fCleanElectrons[i].charge();
	 ePFrelIso[i]=fCleanElectronsPFrelIso[i];
	 ePFabsIso[i]=fCleanElectronsPFabsIso[i];
     
       }
       nMuons=nCleanMuons;

       // GetPFIso(fCleanMuons);
       for(int i=0; i<nCleanMuons; i++){
	 if(i<6){
	   v_m_pt[i]->Fill(fCleanMuons[i].pt()); 
	   v_m_eta[i]->Fill(fCleanMuons[i].eta()); 
	   v_m_phi[i]->Fill(fCleanMuons[i].phi()); 
	 }
	 mpt[i]=fCleanMuons[i].pt();
	 mpx[i]=fCleanMuons[i].px();
	 mpy[i]=fCleanMuons[i].py();
	 mpz[i]=fCleanMuons[i].pz();
	 me[i]=fCleanMuons[i].energy();
	 mcharge[i]=fCleanMuons[i].charge();
	 mPFrelIso[i]=fCleanMuonsPFrelIso[i];
	 mPFabsIso[i]=fCleanMuonsPFabsIso[i];

     
       }
       nPhotons=nCleanPhotons;
       for(int i=0; i<nCleanPhotons; i++){
	 if(i<6){
	   v_ph_pt[i]->Fill(fCleanPhotons[i].pt()); 
	   v_ph_eta[i]->Fill(fCleanPhotons[i].eta()); 
	   v_ph_phi[i]->Fill(fCleanPhotons[i].phi()); 
	 }
	 phpt[i]=fCleanPhotons[i].pt();
	 phpx[i]=fCleanPhotons[i].px();
	 phpy[i]=fCleanPhotons[i].py();
	 phpz[i]=fCleanPhotons[i].pz();
	 phe[i]=fCleanPhotons[i].energy(); 
	 
       }


       GetMCTruth(iEvent);
       //       cout<<"--------------------"<<endl;

       MyTree->Fill();
       entry++;

     }//HasTrigger
   }//GoodRun
}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntupler::beginJob()
{
  if (_isData) {
    char c;
    int n;
    
    ifstream JSONFile(JSONFilename.data());
    if (JSONFile.fail()) {
			cout << "JSON file not found: " << JSONFilename.data() << endl;
			return;
		}
    nGoodRuns=0;
    int run, lb, le;
    bool startlb = false;
    while (! JSONFile.fail() && ! JSONFile.eof()) {
      c=JSONFile.peek();
      while(! JSONFile.fail() && !JSONFile.eof() && !( (c >= '0') && (c <= '9'))) {
				c = JSONFile.get();
				c = JSONFile.peek();
      }
      JSONFile>>n;
      if(n>100000){
	run = n;
	//cout<<run<<endl;
      }
      else{
	if(!startlb){
	  lb = n;
	  startlb = true;
	}
	else{
	  le = n;
	  GoodRuns.push_back(run);
	  GoodLumiStart.push_back(lb);
	  GoodLumiEnd.push_back(le);
	  ++nGoodRuns;
	  startlb = false;
	}
      }
    }
    cout << "Got: " << nGoodRuns << " specified as good." << endl;
  }
  char hTITLE[99];
  char hNAME[99];
  NtupleTree = new TFile(_ntupleTree.data(),"recreate");
  NtupleTree->cd();
  MyTree = new TTree("EvTree", "EvTree");
  SetBranches(MyTree);
  //create output file
  NtuplePlots = new TFile(_ntuplePlots.data(),"recreate");

  h_nGoodVtx = new TH1F("nVtx","Number of Vertices",30,0,30);
  h_zPosGoodVtx = new TH1F("zPosCleanVtx","Z position of the vertices",600,-30,30);
  //some plots to check on, before object cleaning
  sprintf(hTITLE, "Number of good Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nGoodJets = new TH1F("nJetsGood", hTITLE,20,0,20);

  sprintf(hTITLE, "Number of clean Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nCleanJets = new TH1F("nJetsClean", hTITLE,20,0,20);

  sprintf(hTITLE, "Number of good PF Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nGoodPFJets = new TH1F("nPFJetsGood", hTITLE,20,0,20);

  sprintf(hTITLE, "Number of good CA8PF Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nGoodCA8PFJets = new TH1F("nCA8PFJetsGood", hTITLE,20,0,20);

  sprintf(hTITLE, "Number of good CA8PrunedPF Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nGoodCA8PrunedPFJets = new TH1F("nCA8PrunedPFJetsGood", hTITLE,20,0,20);

  sprintf(hTITLE, "Number of good Electrons with Pt>%i and Eta<%i", (int) _ept,(int) _eeta);
  h_nGoodElectrons = new TH1F("nElectronsGood", hTITLE,10,0,10);

  sprintf(hTITLE, "Number of good Muons with Pt>%i and Eta<%i", (int) _mpt,(int) _meta);
  h_nGoodMuons = new TH1F("nMuonsGood", hTITLE,10,0,10);
  
  sprintf(hTITLE, "Number of good Photons with Pt>%i and Eta<%i", (int) _phpt,(int) _pheta);
  h_nGoodPhotons = new TH1F("nPhotonsGood", hTITLE,10,0,10);
  //after object cleaning

  sprintf(hTITLE, "Number of clean PF Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nCleanPFJets = new TH1F("nPFJetsClean", hTITLE,20,0,20);

  sprintf(hTITLE, "Number of clean CA8PF Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nCleanCA8PFJets = new TH1F("nCA8PFJetsClean", hTITLE,20,0,20);
  sprintf(hTITLE, "Number of clean CA8PrunedPF Jet with Pt>%i and Eta<%i", (int) _jetptcut,(int) _etacut);
  h_nCleanCA8PrunedPFJets = new TH1F("nCA8PrunedPFJetsClean", hTITLE,20,0,20);
  
  sprintf(hTITLE, "Number of clean Electrons with Pt>%i and Eta<%i", (int) _ept,(int) _eeta);
  h_nCleanElectrons = new TH1F("nElectronsClean", hTITLE,10,0,10);
  
  sprintf(hTITLE, "Number of clean Muons with Pt>%i and Eta<%i", (int) _mpt,(int) _meta);
  h_nCleanMuons = new TH1F("nMuonsClean", hTITLE,10,0,10);
  
  sprintf(hTITLE, "Number of clean Photons with Pt>%i and Eta<%i", (int) _phpt,(int) _pheta);
  h_nCleanPhotons = new TH1F("nPhotonsClean", hTITLE,10,0,10);

  h_chargedHadronEnergyFraction = new TH1F("chargedHadronEnergyFraction", "chargedHadronEnergyFraction",200,0,1);  
  h_neutralHadronEnergyFraction = new TH1F("neutralHadronEnergyFraction", "neutralHadronEnergyFraction",200,0,1);  
  h_neutralEmEnergyFraction = new TH1F("neutralEmEnergyFraction", "neutralEmEnergyFraction",200,0,1);  
  h_chargedEmEnergyFraction = new TH1F("chargedEmEnergyFraction", "chargedEmEnergyFraction",200,0,1);  
  h_numberOfDaughters = new TH1F("numberOfDaughters", "numberOfDaughters",200,0,200);  
  h_chargedMultiplicity = new TH1F("chargedMultiplicity", "chargedMultiplicity",200,0,200);  

  for(int i=0; i<6; i++)
    {
      if(i<=3){
	_jetptcut=80.0;
	_etacut=2.5;
      }
      else{
	_jetptcut=60.0;
	_etacut=2.5;
      }
      sprintf(hNAME, "PFjet_%i_pt", i);
      sprintf(hTITLE, "PFJetPt of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_PFjet_pt.push_back(new TH1F(hNAME,hTITLE,200,0,1000));
      
      sprintf(hNAME, "PFjet_%i_eta", i);
      sprintf(hTITLE, "PFJetEta of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_PFjet_eta.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "PFjet_%i_phi", i);
      sprintf(hTITLE, "PFJetPhi of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_PFjet_phi.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "PFjet_%i_m", i);
      sprintf(hTITLE, "PFJetMass of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_PFjet_m.push_back(new TH1F(hNAME,hTITLE,200,0,1000));

      sprintf(hNAME, "CA8PFjet_%i_pt", i);
      sprintf(hTITLE, "CA8PFJetPt of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PFjet_pt.push_back(new TH1F(hNAME,hTITLE,200,0,1000));
      
      sprintf(hNAME, "CA8PFjet_%i_eta", i);
      sprintf(hTITLE, "CA8PFJetEta of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PFjet_eta.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "CA8PFjet_%i_phi", i);
      sprintf(hTITLE, "CA8PFJetPhi of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PFjet_phi.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "CA8PFjet_%i_m", i);
      sprintf(hTITLE, "CA8PFJetMass of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PFjet_m.push_back(new TH1F(hNAME,hTITLE,200,0,1000));

      sprintf(hNAME, "CA8PrunedPFjet_%i_pt", i);
      sprintf(hTITLE, "CA8PrunedPFJetPt of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PrunedPFjet_pt.push_back(new TH1F(hNAME,hTITLE,200,0,1000));
      
      sprintf(hNAME, "CA8PrunedPFjet_%i_eta", i);
      sprintf(hTITLE, "CA8PrunedPFJetEta of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PrunedPFjet_eta.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "CA8PrunedPFjet_%i_phi", i);
      sprintf(hTITLE, "CA8PrunedPFJetPhi of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PrunedPFjet_phi.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "CA8PrunedPFjet_%i_m", i);
      sprintf(hTITLE, "CA8PrunedPFJetMass of the %i st Jet with Pt>%i and Eta<%i", i,(int) _jetptcut,(int) _etacut);
      v_CA8PrunedPFjet_m.push_back(new TH1F(hNAME,hTITLE,200,0,1000));

      sprintf(hNAME, "electron_%i_pt", i);
      sprintf(hTITLE, "ElectronPt of the %i st Electron with Pt>%i and Eta<%i", i,(int) _ept,(int) _eeta);
      v_e_pt.push_back(new TH1F(hNAME,hTITLE,200,0,1000));
      
      sprintf(hNAME, "electron_%i_eta", i);
      sprintf(hTITLE, "ElectronEta of the %i st Electron with Pt>%i and Eta<%i", i,(int) _ept,(int) _eeta);
      v_e_eta.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "electron_%i_phi", i);
      sprintf(hTITLE, "ElectronPhi of the %i st Electron with Pt>%i and Eta<%i", i,(int) _ept,(int) _eeta);
      v_e_phi.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "muon_%i_pt", i);
      sprintf(hTITLE, "MuonPt of the %i st Muon with Pt>%i and Eta<%i", i,(int) _mpt,(int) _meta);
      v_m_pt.push_back(new TH1F(hNAME,hTITLE,200,0,1000));

      sprintf(hNAME, "muon_%i_eta", i);
      sprintf(hTITLE, "MuonEta of the %i st Muon with Pt>%i and Eta<%i", i,(int) _mpt,(int) _meta);
      v_m_eta.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "muon_%i_phi", i);
      sprintf(hTITLE, "MuonPhi of the %i st Muon with Pt>%i and Eta<%i", i,(int) _mpt,(int) _meta);
      v_m_phi.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "photon_%i_pt", i);
      sprintf(hTITLE, "PhotonPt of the %i st Photon with Pt>%i and Eta<%i", i,(int) _phpt,(int) _pheta);
      v_ph_pt.push_back(new TH1F(hNAME,hTITLE,200,0,1000));
      
      sprintf(hNAME, "photon_%i_eta", i);
      sprintf(hTITLE, "PhotonEta of the %i st Photon with Pt>%i and Eta<%i", i,(int) _phpt,(int) _pheta);
      v_ph_eta.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
      
      sprintf(hNAME, "photon_%i_phi", i);
      sprintf(hTITLE, "PhotonPhi of the %i st Photon with Pt>%i and Eta<%i", i,(int) _phpt,(int) _pheta);
      v_ph_phi.push_back(new TH1F(hNAME,hTITLE,200,-5,5));
    }     
  h_MET  = new TH1F("PFMet","PFMET",1000,0.,1000.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntupler::endJob() 
{NtuplePlots->cd();
  h_nGoodVtx->Write();
  h_zPosGoodVtx->Write();
  NtuplePlots->cd();
  NtuplePlots->mkdir("Jets");
  NtuplePlots->cd("Jets");
  h_nGoodJets->Write();
  h_nCleanJets->Write();
  h_nGoodPFJets->Write();
  h_nCleanPFJets->Write();
  h_nGoodCA8PFJets->Write();
  h_nCleanCA8PFJets->Write();
  h_nGoodCA8PrunedPFJets->Write();
  h_nCleanCA8PrunedPFJets->Write();
  h_chargedHadronEnergyFraction->Write();
  h_neutralHadronEnergyFraction->Write();
  h_neutralEmEnergyFraction->Write();
  h_chargedEmEnergyFraction->Write();
  h_numberOfDaughters->Write();
  h_chargedMultiplicity->Write();
  for(int i=0; i<6; i++)
    { v_PFjet_m[i]->Write();
      v_PFjet_pt[i]->Write();
      v_PFjet_eta[i]->Write();
      v_PFjet_phi[i]->Write();

      v_CA8PFjet_m[i]->Write();
      v_CA8PFjet_pt[i]->Write();
      v_CA8PFjet_eta[i]->Write();
      v_CA8PFjet_phi[i]->Write();

      v_CA8PrunedPFjet_m[i]->Write();
      v_CA8PrunedPFjet_pt[i]->Write();
      v_CA8PrunedPFjet_eta[i]->Write();
      v_CA8PrunedPFjet_phi[i]->Write();
    }

  NtuplePlots->cd();
  NtuplePlots->mkdir("Electrons");
  NtuplePlots->cd("Electrons");
  h_nGoodElectrons->Write();
  h_nCleanElectrons->Write();
  for(int i=0; i<6; i++)
    { 
      v_e_pt[i]->Write();
      v_e_eta[i]->Write();
      v_e_phi[i]->Write();
    }  
  NtuplePlots->mkdir("Muons");
  NtuplePlots->cd("Muons");
  h_nGoodMuons->Write();
  h_nCleanMuons->Write();
  for(int i=0; i<6; i++)
    { 
      v_m_pt[i]->Write();
      v_m_eta[i]->Write();
      v_m_phi[i]->Write();
    }  

  NtuplePlots->mkdir("Photons");
  NtuplePlots->cd("Photons");
  h_nGoodPhotons->Write();
  h_nCleanPhotons->Write();
  for(int i=0; i<6; i++)
    { 
      v_ph_pt[i]->Write();
      v_ph_eta[i]->Write();
      v_ph_phi[i]->Write();
    }  
  h_MET->Write();

  //Write the tree out
  MyTree->GetCurrentFile()->Write();
  MyTree->GetCurrentFile()->Close();




}

// ------------ method called when starting to processes a run  ------------
void 
Ntupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Ntupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Ntupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Ntupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
Ntupler::UseJson( vector<int > GoodRuns, vector<int > GoodLumiStart,  vector<int > GoodLumiEnd,  Int_t nGoodRuns , int run, int lumis){
 
  if (_isData) {
    GoodRun=kFALSE;
    for (int ii=0;ii<nGoodRuns;++ii){
      if (run == GoodRuns[ii]){
	if (lumis >= GoodLumiStart[ii]
	    && lumis <=GoodLumiEnd[ii])
	  GoodRun = kTRUE;
      }
    }
  } else {  
    GoodRun = kTRUE;
  }
  return;
}

void 
Ntupler::getTriggerDecision(const edm::Event& iEvent, std::map<std::string, bool>& TriggerMap)
{
  edm::Handle<edm::TriggerResults> triggerResults;
  std::string menu = "HLT";
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", menu), triggerResults);
  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
	const std::vector<std::string> &nameList = triggerNames.triggerNames();
	int triggerIndex = -1;
  for (std::map<std::string, bool>::iterator It = TriggerMap.begin(); It != TriggerMap.end(); ++It) {
		size_t length = It->first.size();
		for (unsigned int index = 0; index < nameList.size(); ++index) {
			if (length <= nameList[index].size() && It->first == nameList[index].substr(0, length)) {
				 triggerIndex = index;
				 break;
			}
		}
		It->second = false;
    if (triggerIndex >= 0 && triggerIndex < (int) triggerResults->size()) {
			if (triggerResults->accept(triggerIndex)) {
				It->second = true;
				// cout << "Trigger " << It->first << " matches trigger " <<  nameList[triggerIndex];
				// cout << " and is true.\n";
			}
			return;
    }
  }
	cout << "Bad trigger names ";
  for (std::map<std::string, bool>::iterator It = TriggerMap.begin(); It != TriggerMap.end(); ++It)
		cout << It->first << ", ";
	cout << " with bad index " << triggerIndex << endl;
}


void
Ntupler::GetSUSYpoint(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

       edm::Handle<LHEEventProduct> product;
       bool commentIsThere = iEvent.getByLabel("source", product);
       comments_const_iterator comment;

       string tempString;
       vector<string> tempStrings;
       vector<string> parameters;
       for(comment = product->comments_begin(); comment != product->comments_end();
	   comment++)
	 {
	   if(commentIsThere)
	     {
	       tempString = comment->substr(0,comment->size());
	       tempStrings = split(tempString," ");
	       parameters = split(tempStrings[2], "_");
	       //cout<<parameters[1]<<" "<<parameters[2]<<endl;
	       MSquark= atof(parameters[1].c_str());
	       MLSP=atof(parameters[2].c_str());
	       
	     }
	 }

}


vector<string> Ntupler::split(string fstring, string splitter)
{
  vector<string> returnVector;
  size_t cursor;
  string beforeSplitter;
  string afterSplitter = fstring;
  if(fstring.find(splitter) == string::npos)
    {
      cout<<"No "<<splitter<<" found"<<endl;
      returnVector.push_back(fstring);
      return returnVector;
    }
  else
    {
      while(afterSplitter.find(splitter) != string::npos)
	{
	  cursor = afterSplitter.find(splitter);

	  beforeSplitter = afterSplitter.substr(0, cursor);
	  afterSplitter = afterSplitter.substr(cursor +1, afterSplitter.size());

	  returnVector.push_back(beforeSplitter);

	  if(afterSplitter.find(splitter) == string::npos)
            returnVector.push_back(afterSplitter);
	}
      return returnVector;
    }
}


void 
Ntupler::DoJetID(const edm::Event& iEvent,const edm::EventSetup& iSetup, std::string PatJetType, int type){
  // Jet Handle 
  //  for (unsigned int i =0; i < _patJetType.size(); i++){

  edm::Handle<edm::View<pat::Jet> > PatJets; 
  iEvent.getByLabel(PatJetType, PatJets); 
  
  //comment out jet id for tests
 for (unsigned int j=0; j<PatJets->size(); j++) {
    if (_debug && j==0){
      std::cout<<"--------------run: "<<run<<"--event: "<<event<<"--lumi: "<<lumis<<"--entry: "<<entry<<"-----------"<<endl;
      std::cout << "Number of PatJets: " << PatJets->size() << std::endl;
    }
     
    // Choose jets based on pt and abs(eta)
    if ((*PatJets)[j].pt()        >= 30.0  &&
	fabs((*PatJets)[j].eta()) <= 2.4)
      {
	if (_debug) std::cout << "Jet #" << j << " Passes pT cut of " << _jetptcut 
			      << " and eTa cut of " << _etacut << " with pT: " 
			      <<  (*PatJets)[j].pt() << ", and eTa: " 
			      << (*PatJets)[j].eta() << std::endl;
	bool jetID = false;

       	if (PatJetType == "slimmedJetsCA8"){
	  jetID=true;
	}
	else {
	  //ParticleFlow ID
	  jetID = 
	    ((*PatJets)[j].correctedJet("Uncorrected").neutralHadronEnergyFraction()   < 0.99 && 
	     (*PatJets)[j].correctedJet("Uncorrected").neutralEmEnergyFraction()       < 0.99 &&
	     (*PatJets)[j].correctedJet("Uncorrected").numberOfDaughters()             > 1   &&
	  
	    (fabs((*PatJets)[j].eta())                    > 2.4  ||
	    ((*PatJets)[j].correctedJet("Uncorrected").chargedHadronEnergyFraction() > 0.   &&
	    (*PatJets)[j].correctedJet("Uncorrected").chargedEmEnergyFraction()     < 0.99 &&
	    (*PatJets)[j].correctedJet("Uncorrected").chargedMultiplicity()         > 0.)));
	  
	}
	
	if(jetID){
	  if (_debug) std::cout<<"After JetID"
			       << "Jet #" << j << " Passes pT cut of " << _jetptcut 
			       << " and eTa cut of " << _etacut << " with pT: " 
			       <<  (*PatJets)[j].pt() << ", and eTa: " 
			       << (*PatJets)[j].eta() << std::endl;
	     
	     
	  if(type == 0){
	    fGoodPFJets.push_back((*PatJets)[j]);	   
	    nGoodPFJets++; 
	  }
	  
	  if(type == 1){
	    fGoodCA8PFJets.push_back((*PatJets)[j]);	   
	    nGoodCA8PFJets++; 
	  }
	  else{
	  fGoodJets.push_back((*PatJets)[j]);	   
	  nGoodJets++; 
	  }
	     
	}//JetID
      }//JetKinematics
  }//JetLoop

  if (_debug) std::cout << "Found "<< nGoodJets << " Jets" << std::endl;
  //JEC uncertainty doesn't work
  // const JetCorrector* corrector = JetCorrector::getJetCorrector("JetCorrectionService",iSetup);
  //JetCorrectionUncertainty *jecUnc(0);
  // JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("Jec11_V2_AK5PF_Uncertainty.txt");
  //  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  //(iSetup).get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
  //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);
  //     }
  return;
}


void 
Ntupler::GetTruePileUp(const edm::Event& iEvent){

  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

  std::vector<PileupSummaryInfo>::const_iterator PVI;

  nTruePileUp = -1;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

    int BX = PVI->getBunchCrossing();

    if(BX == 0) { 
      nTruePileUp = PVI->getTrueNumInteractions();
 
      continue;

    }

  }

  return;
}

void 
Ntupler::DoVertexID(const edm::Event& iEvent){

  edm::Handle<reco::VertexCollection>  recVtxs;
  iEvent.getByLabel(_primaryVertex, recVtxs);
  
 
  int CountVtx=0;
  for (size_t i=0; i<recVtxs->size(); ++i)
    if (!((*recVtxs)[i].isFake())) {
      if ( ((*recVtxs)[i].ndof() > 4) &&
           (fabs( (*recVtxs)[i].z()) <= 24) &&
           ((*recVtxs)[i].position().rho() <= 2) ){
        CountVtx++;
	h_zPosGoodVtx->Fill((*recVtxs)[i].z());
      }
    }
  
  if (CountVtx > 0) IsVtxGood = 1;
  h_nGoodVtx->Fill(CountVtx);
       nGoodVtx=CountVtx;

  // cout<< "Vertices " << CountVtx<<" "<< recVtxs->size()<<endl;
}



void
Ntupler::DoElectronID(const edm::Event& iEvent){

  // Get the PF electrons
  edm::Handle< std::vector<pat::Electron> > PatElectrons;
  iEvent.getByLabel("slimmedElectrons", PatElectrons);

  // Will need the vertices
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_primaryVertex, primaryVertices);

  // Loop over all electrons
  for (std::vector<pat::Electron>::const_iterator Electron = PatElectrons->begin(); Electron != PatElectrons->end(); ++Electron) {
    
    //calculate isolation

    
    if (Electron->pt() < 10.0) continue;
    if (fabs(Electron->superCluster()->eta()) > 2.5) continue;
    if (fabs(Electron->superCluster()->eta()) > 1.4442 && fabs(Electron->superCluster()->eta())< 1.566) continue;
    
    
    if (fabs(Electron->gsfTrack()->dxy((*primaryVertices)[0].position())) > 0.02) continue;
      

    
    // Barrel & endcap.. if in between say bye bye
    if (fabs(Electron->superCluster()->eta()) < 1.479) {
    if (Electron->scSigmaIEtaIEta() > 0.01) continue;
    if (fabs(Electron->deltaPhiSuperClusterTrackAtVtx()) > 0.06) continue;
    if (fabs(Electron->deltaEtaSuperClusterTrackAtVtx()) > 0.004) continue;
    if (Electron->hcalOverEcal() > 0.12) continue;
    } else if (fabs(Electron->superCluster()->eta()) > 1.479) {
    if (Electron->scSigmaIEtaIEta() > 0.03) continue;
    if (fabs(Electron->deltaPhiSuperClusterTrackAtVtx()) > 0.03) continue;
    if (fabs(Electron->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
    if (Electron->hcalOverEcal() > 0.1) continue;
    } else {
    // crackhead
    std::cerr << "ERROR: Electron selection. you should never see this message" << std::endl;
    continue;
    }
    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel("offlineBeamSpot", bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();
    
    edm::Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel("reducedEgamma","reducedConversions", hConversions);
    

    //reco::GsfElectron const *gsfele = dynamic_cast<reco::GsfElectron const *>((Electron->originalObjectRef().get()));
    bool isconversion = ConversionTools::hasMatchedConversion( *Electron , hConversions, beamspot.position());
    if(isconversion) continue;  
  
    
    fGoodElectrons.push_back(*Electron);
    nGoodElectrons++;


  } // pat electron loop

  return;
}










void
Ntupler::DoMuonID(const edm::Event& iEvent){
  
  // Get muon collection
  edm::Handle< std::vector<pat::Muon> > PatMuons; 
  iEvent.getByLabel("slimmedMuons", PatMuons);
  

  edm::Handle<reco::VertexCollection>  primaryVertices;
  iEvent.getByLabel( _primaryVertex, primaryVertices);
  //  cout<<primaryVertices->size()<<" "<<((*primaryVertices)[0].position())<<endl;

  
  for (std::vector<pat::Muon>::const_iterator Muon = PatMuons->begin(); Muon != PatMuons->end(); ++Muon) {
    
   
      // cout<<"is good muon: "<<muon::isTightMuon(*Muon,(*primaryVertices)[0])<<endl;
    if (!Muon->isGlobalMuon()) continue;
    if (!Muon->isPFMuon()) continue;
    if (!(Muon->globalTrack()->normalizedChi2() < 10.)) continue;
    if (!(Muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0)) continue;
    if (!(Muon->numberOfMatchedStations() > 1)) continue;
    //muonbestTrack doesn't work properly with miniAOD, segmentation violation seems to occur for certain muons
    //use innterTrack instead
    if (!(fabs(Muon->innerTrack()->dxy((*primaryVertices)[0].position())) < 0.2 )) continue;
    if (!(fabs(Muon->innerTrack()->dz((*primaryVertices)[0].position())) < 0.5)) continue;
    if (!(Muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0)) continue;
    if (!(Muon->track()->hitPattern().trackerLayersWithMeasurement() > 5)) continue;
  
      
    
    

    if ( Muon->pt() < 10.0 ) continue;

    if ( fabs(Muon->eta()) > 2.4 ) continue;



    
    if (!Muon->isGlobalMuon()) continue;
      fGoodMuons.push_back(*Muon);    
      nGoodMuons++;
    
  }
  
  return;
}

void
Ntupler::DoMETID(const edm::Event& iEvent){
  Handle< vector<MET> >      MetColl;
  iEvent.getByLabel(_METtype,  MetColl);
  fMET=(*MetColl)[0];
  
  
  return;
}
void
Ntupler::DoCleanUp(const edm::Event& iEvent,vector<pat::Muon >fGoodMuons,vector<pat::Electron >fGoodElectrons,vector<pat::Photon >fGoodPhotons,vector<pat::Jet >fGoodJets){
  for (size_t im = 0; im != fGoodMuons.size(); ++im) {
   float pfIsoCharged = fGoodMuons[im].pfIsolationR03().sumChargedHadronPt;
    float pfIsoNeutral = fGoodMuons[im].pfIsolationR03().sumNeutralHadronEt;
    float pfIsoPhoton  = fGoodMuons[im].pfIsolationR03().sumPhotonEt;
    float pfIsoPU      = fGoodMuons[im].pfIsolationR03().sumPUPt;      
    float  absoluteIsolation = pfIsoCharged + max((float) 0.0,(float) (pfIsoNeutral + pfIsoPhoton- 0.5*pfIsoPU ));
    float relIso =  absoluteIsolation/fGoodMuons[im].pt();
    fCleanMuonsPFrelIso.push_back(relIso);
    fCleanMuonsPFabsIso.push_back(relIso*fGoodMuons[im].pt());
 
    fCleanMuons.push_back(fGoodMuons[im] );
    nCleanMuons++;
  }
  // Keep non-overlapping electrons
  for (size_t ie = 0; ie != fGoodElectrons.size(); ++ie) {
    bool HasOverlap = false;
    TLorentzVector Electron(fGoodElectrons[ie].px(), fGoodElectrons[ie].py(), fGoodElectrons[ie].pz(), fGoodElectrons[ie].energy()); 
    for (size_t im = 0; im != fGoodMuons.size(); ++im) {
      TLorentzVector Muon(fGoodMuons[im].px(), fGoodMuons[im].py(), fGoodMuons[im].pz(), fGoodMuons[im].p()); 
      if (Muon.DeltaR( Electron ) < 0.4) {
        HasOverlap = true;
	if(_debug && HasOverlap) cout<<"Overlap Electron Muon with pt: "<<fGoodElectrons[ie].pt()<<" eta: "<<fGoodElectrons[ie].eta()<<endl;

      }
    }
    if (!HasOverlap) {
      
      //calculate isolation
    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(_rhoIsoInputTag, rhoHandle);
    double rhoIso = std::max(*(rhoHandle.product()), 0.0);
    double scEta = fGoodElectrons[ie].superCluster()->eta();
    double AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, scEta, ElectronEffectiveArea::kEleEAData2012);   
    double chIso = fGoodElectrons[ie].chargedHadronIso();
    double nhIso = fGoodElectrons[ie].neutralHadronIso();
    double phIso = fGoodElectrons[ie].photonIso();
    double relIso = ( chIso + max(0.0, nhIso + phIso - rhoIso*AEff) )/ fGoodElectrons[ie].ecalDrivenMomentum().pt(); 
    if(_debug)
      {    
	cout<<"chIso "<<chIso<<" nhIso "<<nhIso<<" phIso "<< phIso<<" rhoIso "<< rhoIso<< "AEff "<< AEff<<endl;
	cout<<"relIso "<<relIso<<endl;
      }
    fCleanElectronsPFrelIso.push_back(relIso);
    fCleanElectronsPFabsIso.push_back(relIso*fGoodElectrons[ie].ecalDrivenMomentum().pt());
    fCleanElectrons.push_back( fGoodElectrons[ie]);
    nCleanElectrons++;
    }
  }

  // Keep non-overlapping photons
  for (size_t ip = 0; ip != fGoodPhotons.size(); ++ip) {
    bool HasOverlap = false;
    TLorentzVector Photon(fGoodPhotons[ip].px(), fGoodPhotons[ip].py(), fGoodPhotons[ip].pz(), fGoodPhotons[ip].energy()); 
    for (size_t ie = 0; ie != fGoodElectrons.size(); ++ie) {
      TLorentzVector Electron(fGoodElectrons[ie].px(), fGoodElectrons[ie].py(), fGoodElectrons[ie].pz(), fGoodElectrons[ie].energy());
      if (Electron.DeltaR(Photon) < 0.4) {
        HasOverlap = true;
      }
    }
    if (!HasOverlap) {
      fCleanPhotons.push_back( fGoodPhotons[ip] );
      nCleanPhotons++;
    }
  }
  // Keep non-overlapping jets
  nBJets=0;
  for (size_t ij = 0; ij != fGoodJets.size(); ++ij) {
    bool HasOverlap = false;
    TLorentzVector Jet(fGoodJets[ij].px(), fGoodJets[ij].py(), fGoodJets[ij].pz(), fGoodJets[ij].energy()); 
    for (size_t ie = 0; ie != fCleanElectrons.size(); ++ie) {
      TLorentzVector Electron(fCleanElectrons[ie].px(), fCleanElectrons[ie].py(), fCleanElectrons[ie].pz(), fCleanElectrons[ie].energy());
      if (Electron.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t ip = 0; ip != fCleanPhotons.size(); ++ip) {
      TLorentzVector Photon(fCleanPhotons[ip].px(), fCleanPhotons[ip].py(), fCleanPhotons[ip].pz(), fCleanPhotons[ip].energy());
      if (Photon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t im = 0; im != fCleanMuons.size(); ++im) {
      TLorentzVector Muon(fCleanMuons[im].px(), fCleanMuons[im].py(), fCleanMuons[im].pz(), fCleanMuons[im].energy());
      if (Muon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }

    if (!HasOverlap) {
      fCleanJets.push_back( fGoodJets[ij] );
      nCleanJets++;
    }
  }


  for (size_t ij = 0; ij != fGoodPFJets.size(); ++ij){
    bool HasOverlap = false;
    TLorentzVector Jet(fGoodPFJets[ij].px(), fGoodPFJets[ij].py(), fGoodPFJets[ij].pz(), fGoodPFJets[ij].energy()); 
    for (size_t ie = 0; ie != fCleanElectrons.size(); ++ie) {
      TLorentzVector Electron(fCleanElectrons[ie].px(), fCleanElectrons[ie].py(), fCleanElectrons[ie].pz(), fCleanElectrons[ie].energy());
      if (Electron.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t ip = 0; ip != fCleanPhotons.size(); ++ip) {
      TLorentzVector Photon(fCleanPhotons[ip].px(), fCleanPhotons[ip].py(), fCleanPhotons[ip].pz(), fCleanPhotons[ip].energy());
      if (Photon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t im = 0; im != fCleanMuons.size(); ++im) {
      TLorentzVector Muon(fCleanMuons[im].px(), fCleanMuons[im].py(), fCleanMuons[im].pz(), fCleanMuons[im].energy());
      if (Muon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }

    if (!HasOverlap) {
      fCleanPFJets.push_back( fGoodPFJets[ij] );
      if (fGoodPFJets[ij].bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74) nBJets++;
      nCleanPFJets++;
    }
  }

  for (size_t ij = 0; ij != fGoodCA8PFJets.size(); ++ij){
    bool HasOverlap = false;
    TLorentzVector Jet(fGoodCA8PFJets[ij].px(), fGoodCA8PFJets[ij].py(), fGoodCA8PFJets[ij].pz(), fGoodCA8PFJets[ij].energy()); 
    for (size_t ie = 0; ie != fCleanElectrons.size(); ++ie) {
      TLorentzVector Electron(fCleanElectrons[ie].px(), fCleanElectrons[ie].py(), fCleanElectrons[ie].pz(), fCleanElectrons[ie].energy());
      if (Electron.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t ip = 0; ip != fCleanPhotons.size(); ++ip) {
      TLorentzVector Photon(fCleanPhotons[ip].px(), fCleanPhotons[ip].py(), fCleanPhotons[ip].pz(), fCleanPhotons[ip].energy());
      if (Photon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }
    for (size_t im = 0; im != fCleanMuons.size(); ++im) {
      TLorentzVector Muon(fCleanMuons[im].px(), fCleanMuons[im].py(), fCleanMuons[im].pz(), fCleanMuons[im].energy());
      if (Muon.DeltaR(Jet) < 0.4) {
        HasOverlap = true;
      }
    }

    if (!HasOverlap) {
      fCleanCA8PFJets.push_back( fGoodCA8PFJets[ij] );
      nCleanCA8PFJets++;
    }
  }

  return;
}
void 
Ntupler::GetMCTruth(const edm::Event& iEvent){
  if(!_isData){

    edm::Handle< std::vector<reco::GenParticle> > GenParticles; 
    iEvent.getByLabel("prunedGenParticles", GenParticles);  
       nGenPart=(*GenParticles).size();
    for (unsigned int p=0; p<(*GenParticles).size(); p++) { 
      //      if((*GenParticles)[p].status()==3){
      //cout<<p<<endl; 
      //use only that hard process
      if(p<200){
      pdgID[p]=(*GenParticles)[p].pdgId();
       
      MCpx[p]=(*GenParticles)[p].px();
      MCpy[p]=(*GenParticles)[p].py();
      MCpz[p]=(*GenParticles)[p].pz();
      MCe[p]=(*GenParticles)[p].energy();
      for (unsigned int ind = 0; ind < 200 && ind < (*GenParticles).size(); ind++) {
	if ((*GenParticles)[p].mother() == &((*GenParticles)[ind])) {
	  MCmotherind[p] = ind;
	  break;
	}
      }
      //}
      }

     
    }
  }
  return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(Ntupler);

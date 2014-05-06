#include <memory>
#include <iostream>
#include <TString.h> 
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH3D.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TBranch.h>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/HandleBase.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include <fstream>
#include <vector>
#include <functional>

//CMSSW includes
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"

#include "DesyAnalysis/Ntupler/interface/NtpReader.h"


using namespace std;
using namespace edm;
using namespace pat;

//                                                                                                                                                                                
// class declaration                                                                                                                                                             
//                                                                                                                                                                                  

class Ntupler : public edm::EDAnalyzer, public NtpReader {
 public:
  explicit Ntupler(const edm::ParameterSet&);
  ~Ntupler();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------                                                                                                                         
  void UseJson( vector<int > GoodRuns, vector<int > GoodLumiStart,  vector<int > GoodLumiEnd,  Int_t nGoodRuns, int run, int lumis );
  void DoJetID(const edm::Event&,const edm::EventSetup& iSetup, std::string PatJetType);
  void GetTruePileUp(const edm::Event& iEvent);
  void DoVertexID(const edm::Event&);
  void DoElectronID(const edm::Event&);
  void DoMuonID(const edm::Event&);
  void DoPhotonID(const edm::Event&);
  void DoMETID(const edm::Event&);
  void DoCleanUp(vector<Muon >fGoodMuons,vector<Electron >fGoodElectrons,vector<Photon >fGoodPhotons,vector<Jet >fGoodJets);
  void GetPFIso(vector<Muon >fCleanMuons);
  void GetMCTruth(const edm::Event&);
  void MakeTriplets(std::vector<pat::Jet >fCleanJets); 
  void getTriggerDecision(const edm::Event&, std::map<std::string, bool>&);
  void getTriggerDecision2(const edm::Event&, std::map<std::string, bool>&);


  void GetSUSYpoint(const edm::Event&, const edm::EventSetup&);
 bool SUSYfilter(const edm::Event&, const edm::EventSetup&);
 std::vector<std::string> split(std::string fstring, std::string splitter);
 typedef std::vector<std::string>::const_iterator comments_const_iterator;

  
  
  // ----------member data ---------------------------
  ofstream run_evt_ee, run_evt_mm, run_evt_em;
  
  double _sumPtMin, _etacut, _jetptcut, _eeta,_ept,_meta,_mpt,_pheta,_phpt;
  bool _isData;
  bool _isSUSY;
  bool _noTripletBtag;
  bool _debug;
  int _nbTagsMin; 
  int _nbTagsMax;      int _nTripletBtagsMin, _nTripletBtagsMax;
  bool _doBtagEff;
  std::string _ntuplePlots, _ntupleTree, corStep, corFlvr;   
  std::vector<std::string> _patJetType;
  std::string PatJetType;
  std::string _primaryVertex;
  std::string _jecAdj;
  std::string _jetCorrectionService;
  std::string _METtype;
  int _njetsMin, _njetsMax;
  std::string _htTrigger;
  std::vector<std::string>          fTriggerNamesSel;
  std::vector<std::string>          fTriggerNamesBase;
  std::vector<std::string>          fTriggerNamesBase2;
  std::map<std::string, bool>       fTriggerMap;
  std::map<std::string, bool>       fTriggerMap2;
  std::map<std::string, bool>       fTriggerMapBase2;
  
  int countE;  
  
  TFile* NtuplePlots; 
  TFile* NtupleTree; 
  //varibales to test for the JSON file                                                                                            
  bool isdata;
  string JSONFilename;
  TFile* JSONFile;
  vector<int > GoodRuns;
  vector<int > GoodLumiStart;
  vector<int >  GoodLumiEnd;
     std::vector<pat::Jet>  fCleanPFJets;  

  Int_t nGoodRuns;  
  Bool_t GoodRun;
  Int_t nGoodJets;
  Int_t nCleanJets;
  Int_t nGoodPFJets;
  Int_t nCleanPFJets;
  Int_t nGoodCA8PFJets;
  Int_t nCleanCA8PFJets;
  Int_t nGoodCA8PrunedPFJets;
  Int_t nCleanCA8PrunedPFJets;

  Int_t nBJets;
  Int_t nGoodElectrons;
  Int_t nCleanElectrons;
  
  Int_t nGoodMuons;
  Int_t nCleanMuons;
  
  Int_t nGoodPhotons;
  Int_t nCleanPhotons;
  
  Int_t IsVtxGood;
  TH1F* h_chargedHadronEnergyFraction;
  TH1F* h_neutralHadronEnergyFraction;
  TH1F* h_neutralEmEnergyFraction;
  TH1F* h_chargedEmEnergyFraction;
  TH1F* h_numberOfDaughters;
  TH1F* h_chargedMultiplicity;
  TH1F* h_nGoodVtx;
  TH1F* h_zPosGoodVtx;
  TH1F* h_NumEvtCut;

  TH1F* h_nCleanPFJets;
  TH1F* h_nGoodPFJets;
  TH1F* h_nCleanCA8PFJets;
  TH1F* h_nGoodCA8PFJets;
  TH1F* h_nCleanCA8PrunedPFJets;
  TH1F* h_nGoodCA8PrunedPFJets;

  TH1F* h_nGoodElectrons;
  TH1F* h_nCleanElectrons;
  TH1F* h_nGoodMuons;
  TH1F* h_nCleanMuons;
  TH1F* h_nGoodPhotons;
  TH1F* h_nCleanPhotons;
  TH1F* h_MET;
  //some event check
  TH1F* h_DiElectronMass;
  TH1F* h_DiMuonMass;
  TH1F* h_ElectronMuonMass;
  
  //redefining everything with vectors
  vector <string > VarList;
  vector <string > CutList;
  //define variables that we need for the triplets
  std::vector <std::vector<pat::Jet > > Triplet;
  std::vector<double> nMCBtagTriplet;
  std::vector<double> sumScalarPtTriplet;
  std::vector<double> sumVectorPtTriplet;
  std::vector<double> massTriplet;
  std::vector<double> massQuad;
  pat::Jet AntiTripletHighestJet;
  pat::Jet Jet1; pat::Jet Jet2; pat::Jet Jet3; pat::Jet Jet4; 

  Int_t q;
  //define our own jet collection (just in case if we want to clean up also)
  std::vector<pat::Jet >      fGoodJets;
  std::vector<pat::Jet >      fCleanJets;


  //  std::vector<pat::Jet >      fGoodPFJets;
  //std::vector<pat::Jet >      fCleanPFJets;
  //std::vector<pat::Jet >      fGoodCA8PFJets;
  //std::vector<pat::Jet >      fCleanCA8PFJets;
  //std::vector<pat::Jet >      fGoodCA8PrunedPFJets;
  //std::vector<pat::Jet >      fCleanCA8PrunedPFJets;
  std::vector<pat::Jet >      fBJets;
  std::vector<pat::Jet >      fNoBJets;
  //electron collection
  std::vector<pat::Electron> fGoodElectrons;
  std::vector<pat::Electron> fCleanElectrons;
  
  //muon collection
  std::vector<pat::Muon>     fGoodMuons;
  std::vector<pat::Muon>      fCleanMuons;
  std::vector<float >      fCleanMuonsPFIso;
  std::vector<pat::Muon>      fFakeMuons;
  
  //photon collection
  std::vector<pat::Photon>     fGoodPhotons;
  std::vector<pat::Photon>      fCleanPhotons;
  
  //vertices
  std::vector<reco::VertexCollection>     fGoodVtx;
  //MET
  //TLorentzVector  fMET;
  pat::MET  fMET;
  //some plots for the 6 highest jets in the event
  std::vector<TH1F*> v_PFjet_pt;
  std::vector<TH1F*> v_PFjet_eta;
  std::vector<TH1F*> v_PFjet_phi;
  std::vector<TH1F*> v_PFjet_m;

  std::vector<TH1F*> v_CA8PFjet_pt;
  std::vector<TH1F*> v_CA8PFjet_eta;
  std::vector<TH1F*> v_CA8PFjet_phi;
  std::vector<TH1F*> v_CA8PFjet_m;

  std::vector<TH1F*> v_CA8PrunedPFjet_pt;
  std::vector<TH1F*> v_CA8PrunedPFjet_eta;
  std::vector<TH1F*> v_CA8PrunedPFjet_phi;
  std::vector<TH1F*> v_CA8PrunedPFjet_m;

  std::vector<TH1F*> v_e_pt;
  std::vector<TH1F*> v_e_eta;
  std::vector<TH1F*> v_e_phi;
  
  std::vector<TH1F*> v_m_pt;
  std::vector<TH1F*> v_m_eta;
  std::vector<TH1F*> v_m_phi;
  
  std::vector<TH1F*> v_ph_pt;
  std::vector<TH1F*> v_ph_eta;
  std::vector<TH1F*> v_ph_phi;
  //some other event variables
  

  TH1F* h_LeptonPt;
  TH1F* h_ElectronPt;
  TH1F* h_MuonPt;
  
  //stuff for the tree
  TTree*  MyTree;
  
};





#ifndef GUARD_NtpReader_h
#define GUARD_NtpReader_h

#include <iostream>
#include <vector>

#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

class NtpReader
{
  public:
    NtpReader ();
    NtpReader (TTree*);
    NtpReader (std::vector<TString>&, bool const);
    ~NtpReader ();

    bool IsData ();

    void SetBranches (TTree*);
    void SetBranchAddresses ();
    int  GetEntry (int const);


  private:
    TChain fChain;
    bool   fIsData;



  // These members are in the ntuple
  public:
    int   run;
    int   event;
    int   lumis;
    int   entry;
    int   nTriplets;
    int   nGoodVtx;
    Bool_t HasBaseTrigger;
    Bool_t HasSelTrigger;
    Bool_t DataIs;
    double MSquark;
    double MLSP;
    float nTruePileUp;


    //TEST
    int NumJetsTest;
    int   mcsize;
    int   nGenPart;
    int   pdgID[200];
    float MCpx[200];
    float MCpy[200];
    float MCpz[200];
    float MCe[200];
    int   MCmotherind[200];

    // int MaxSize=500;
    int   nPFJets; 
    float jet_PF_px[200];
    float jet_PF_py[200];
    float jet_PF_pz[200];
    float jet_PF_e[200];
    float jet_PF_pt[200];
    float jet_PF_jec_unc[200];
    float  jet_PF_et[200];
    float  jet_PF_eta[200];
    float  jet_PF_phi[200];
    float  jet_PF_theta[200];
    float jet_PF_mass[200];
    float jet_PF_NeutralHad[200];
    float jet_PF_area[200];
    float jet_PF_nJetDaughters[200];
    int jet_PF_PartonFlav[200];
    int jet_PF_JetMom[200];
    float bdiscCSV_PF[200];

    float  bdiscJP_PF[200];
    float bdiscTCHE_PF[200];
    float bdiscTCHP_PF[200];
    float bdiscSSVHE_PF[200];
    float bdiscSSSVHP_PF[200];

    int   nCA8PFJets; 
    float jet_CA8PF_px[200];
    float jet_CA8PF_py[200];
    float jet_CA8PF_pz[200];
    float jet_CA8PF_e[200];
    float jet_CA8PF_pt[200];
    float  jet_CA8PF_et[200];
    float  jet_CA8PF_eta[200];
    float  jet_CA8PF_phi[200];
    float  jet_CA8PF_theta[200];
    float jet_CA8PF_mass[200];
    float jet_CA8PF_area[200];
    float jet_CA8PF_nJetDaughters[200];
    float bdiscCSV_CA8PF[200];
    float  bdiscJP_CA8PF[200];
    float bdiscTCHE_CA8PF[200];
    float bdiscTCHP_CA8PF[200];
    float bdiscSSVHE_CA8PF[200];
    float bdiscSSSVHP_CA8PF[200];

    int   nCA8PrunedPFJets; 
    float jet_CA8PrunedPF_px[200];
    float jet_CA8PrunedPF_py[200];
    float jet_CA8PrunedPF_pz[200];
    float jet_CA8PrunedPF_e[200];
    float jet_CA8PrunedPF_pt[200];
    float  jet_CA8PrunedPF_et[200];
    float  jet_CA8PrunedPF_eta[200];
    float  jet_CA8PrunedPF_phi[200];
    float  jet_CA8PrunedPF_theta[200];
    float jet_CA8PrunedPF_mass[200];
    float jet_CA8PrunedPF_area[200];
    float jet_CA8PrunedPF_nJetDaughters[200];
    float bdiscCSV_CA8PrunedPF[200];
    float  bdiscJP_CA8PrunedPF[200];
    float bdiscTCHE_CA8PrunedPF[200];
    float bdiscTCHP_CA8PrunedPF[200];
    float bdiscSSVHE_CA8PrunedPF[200];
    float bdiscSSSVHP_CA8PrunedPF[200];
    float jet_CA8PrunedPF_subJet1_mass[200];
    float jet_CA8PrunedPF_subJet2_mass[200];
    float jet_CA8PrunedPF_subJet3_mass[200];

    float jet_CA8PrunedPF_subJet1_px[200];
    float jet_CA8PrunedPF_subJet2_px[200];
    float jet_CA8PrunedPF_subJet3_px[200];

    float jet_CA8PrunedPF_subJet1_py[200];
    float jet_CA8PrunedPF_subJet2_py[200];
    float jet_CA8PrunedPF_subJet3_py[200];

    float jet_CA8PrunedPF_subJet1_pz[200];
    float jet_CA8PrunedPF_subJet2_pz[200];
    float jet_CA8PrunedPF_subJet3_pz[200];

    float jet_CA8PrunedPF_subJet1_e[200];
    float jet_CA8PrunedPF_subJet2_e[200];
    float jet_CA8PrunedPF_subJet3_e[200];



    int   nElectrons;
    float epx[200];
    float epy[200];
    float epz[200];
    float ee[200];
    float ept[200];
    int   echarge[200];

    int   nMuons;
    float mpx[200];
    float mpy[200];
    float mpz[200];
    float me[200];
    float mpt[200];
    float mcharge[200];
    float mPFIso[200];
    int   nPhotons;
    float phpx[200];
    float phpy[200];
    float phpz[200];
    float phe[200];
    float phpt[200];

    float pfMET;
    float pfMETphi;

    //the triplets here could handle 30 jets  with combinaotriz

    float triplet_jet1pt[5000];
    float triplet_jet1px[5000];
    float triplet_jet1py[5000];
    float triplet_jet1pz[5000];
    float triplet_jet1e[5000];

    float triplet_jet2pt[5000];
    float triplet_jet2px[5000];
    float triplet_jet2py[5000];
    float triplet_jet2pz[5000];
    float triplet_jet2e[5000];

    float triplet_jet3pt[5000];
    float triplet_jet3px[5000];
    float triplet_jet3py[5000];
    float triplet_jet3pz[5000];
    float triplet_jet3e[5000];

    float triplet_mass[5000];
    float triplet_sumScalarPt[5000];
    float triplet_sumVectorPt[5000];



};






#endif

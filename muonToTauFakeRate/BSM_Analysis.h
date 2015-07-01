////////////////////////////////////////////////////////////////////
// Developer: Andres Florez, Universidad de los Andes, Colombia. //
///////////////////////////////////////////////////////////////////


#ifndef BSM_Analysis_h
#define BSM_Analysis_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TDirectory.h>
using namespace std;

class BSM_Analysis {
public :
   BSM_Analysis(TFile*, TDirectory* dir[], int nDir, char*);
   ~BSM_Analysis();

   // create Histo maps
   void createHistoMaps (int);
   // Define maps for histograms
   std::map<unsigned int, TH1*> _hmap_events;
   // For muons
   std::map<unsigned int, TH1*> _hmap_probe_tau_pT;
   std::map<unsigned int, TH1*> _hmap_probe_tau_eta;
   std::map<unsigned int, TH1*> _hmap_probe_tau_phi;
   std::map<unsigned int, TH1*> _hmap_diLepton_mass;
   std::map<unsigned int, TH1*> _hmap_probe_tau_pT_fail;
   std::map<unsigned int, TH1*> _hmap_probe_tau_eta_fail;
   std::map<unsigned int, TH1*> _hmap_diLepton_mass_fail;
   std::map<unsigned int, TH1*> _hmap_probe_tau_phi_fail;
   
// Define Branches
   void setBranchAddress(TTree* BOOM);
   vector<string>  *Trigger_names;
   vector<int>     *Trigger_decision;
   vector<double>  *Muon_pt;
   vector<double>  *Muon_eta;
   vector<double>  *Muon_phi;
   vector<double>  *Muon_p;
   vector<double>  *Muon_energy;
   vector<double>  *Muon_charge;
   vector<bool>    *Muon_tight;
   vector<bool>    *Muon_soft;
   vector<bool>    *Muon_pf;
   vector<double>  *Muon_isoCharged;
   vector<double>  *Muon_isoSum;
   vector<double>  *Muon_isoCharParPt;
   vector<double>  *Muon_chi2;
   vector<double>  *Muon_validHits;
   vector<double>  *Muon_validHitsInner;
   vector<double>  *Muon_matchedStat;
   vector<double>  *Muon_dxy;
   vector<double>  *Muon_TLayers;
   vector<double>  *Muon_dz;
   vector<double>  *Muon_isoNeutralHadron;
   vector<double>  *Muon_isoPhoton;
   vector<double>  *Muon_isoPU;
   vector<double>  *patElectron_pt;
   vector<double>  *patElectron_eta;
   vector<double>  *patElectron_phi;
   vector<double>  *patElectron_energy;
   vector<double>  *patElectron_charge;
   vector<int>     *patElectron_isPassVeto;
   vector<int>     *patElectron_isPassLoose;
   vector<int>     *patElectron_isPassMedium;
   vector<int>     *patElectron_isPassTight;
   vector<int>     *patElectron_isPassHEEPId;
   vector<double>  *patElectron_d0;
   vector<double>  *patElectron_dz;
   vector<int>     *patElectron_expectedMissingInnerHits;
   vector<int>     *patElectron_passConversionVeto;
   vector<double>  *patElectron_isoChargedHadrons;
   vector<double>  *patElectron_isoNeutralHadrons;
   vector<double>  *patElectron_isoPhotons;
   vector<double>  *patElectron_isoPU;
   vector<double>  *Tau_eta;
   vector<double>  *Tau_phi;
   vector<double>  *Tau_pt;
   vector<double>  *Tau_energy;
   vector<double>  *Tau_charge;
   vector<int>     *Tau_decayModeFinding;
   vector<int>     *Tau_decayModeFindingNewDMs;
   vector<double>  *Tau_chargedIsoPtSum;
   vector<double>  *Tau_neutralIsoPtSum;
   vector<int>     *Tau_againstMuonTight3;
   vector<int>     *Tau_againstElectronMVATightMVA5;
   vector<double>  *Tau_nProngs;
   vector<double>  *Tau_puCorrPtSum;
   vector<int>     *Tau_byLooseCombinedIsolationDeltaBetaCorr;
   vector<int>     *Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byMediumCombinedIsolationDeltaBetaCorr;
   vector<int>     *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byTightCombinedIsolationDeltaBetaCorr;
   vector<int>     *Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byLooseIsolationMVA3newDMwLT;
   vector<int>     *Tau_byLooseIsolationMVA3newDMwoLT;
   vector<int>     *Tau_byLooseIsolationMva3oldDMwLT;
   vector<int>     *Tau_byLooseIsolationMVA3oldDMwoLT;
   vector<int>     *Tau_byMediumIsolationMVA3newDMwLT;
   vector<int>     *Tau_byMediumIsolationMVA3newDMwoLT;
   vector<int>     *Tau_byMediumIsolationMva3oldDMwLT;
   vector<int>     *Tau_byMediumIsolationMVA3oldDMwoLT;
   vector<int>     *Tau_byTightIsolationMVA3newDMwLT;
   vector<int>     *Tau_byTightIsolationMVA3newDMwoLT;
   vector<int>     *Tau_byTightIsolationMva3oldDMwLT;
   vector<int>     *Tau_byTightIsolationMVA3oldDMwoLT;
   vector<int>     *Tau_againstMuonLoose2;
   vector<int>     *Tau_againstMuonLoose3;
   vector<int>     *Tau_againstMuonTight2;
   vector<int>     *Tau_againstElectronMVALooseMVA5;
   vector<int>     *Tau_againstElectronMVAMediumMVA5;
   vector<int>     *Tau_byVLooseCombinedIsolationDeltaBetaCorr;
   vector<double>  *Tau_leadChargedCandPt;
   vector<double>  *Tau_leadChargedCandCharge;
   vector<double>  *Tau_leadChargedCandEta;
   vector<double>  *Tau_leadChargedCandPhi;
   vector<double>  *Jet_pt;
   vector<double>  *Jet_eta;
   vector<double>  *Jet_phi;
   vector<double>  *Jet_energy;
   vector<double>  *Jet_bDiscriminator;
   vector<double>  *Jet_mass;
   vector<double>  *Jet_neutralHadEnergyFraction;
   vector<double>  *Jet_neutralEmEmEnergyFraction;
   vector<double>  *Jet_chargedHadronEnergyFraction;
   vector<double>  *Jet_chargedEmEnergyFraction;
   vector<double>  *Jet_muonEnergyFraction;
   vector<double>  *Jet_electronEnergy;
   vector<double>  *Jet_photonEnergy;
   vector<double>  *UncorrJet_pt;
   Int_t           npuVertices;
   Float_t         trueInteractions;
   Int_t           ootnpuVertices;
   Int_t           npuVerticesp1;
   Int_t           bestVertices;
   Double_t        Met_pt;
   Double_t        Met_sumEt;
   Double_t        Met_phi;
   Double_t        Met_px;
   Double_t        Met_py;
   Double_t        Met_pz;
   Double_t        Gen_Met;
   Double_t        Met_shiftedPtUp;
   Double_t        Met_shiftedPtDown;

   // List of branches
   TBranch        *b_Trigger_names;
   TBranch        *b_Trigger_decision;
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_p;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_tight;   //!
   TBranch        *b_Muon_soft;   //!
   TBranch        *b_Muon_pf;   //!
   TBranch        *b_Muon_isoCharged;   //!
   TBranch        *b_Muon_isoSum;   //!
   TBranch        *b_Muon_isoCharParPt;   //!
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_validHits;   //!
   TBranch        *b_Muon_validHitsInner;   //!
   TBranch        *b_Muon_matchedStat;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_TLayers;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_isoNeutralHadron;   //!
   TBranch        *b_Muon_isoPhoton;   //!
   TBranch        *b_Muon_isoPU;   //!
   TBranch        *b_patElectron_pt;   //!
   TBranch        *b_patElectron_eta;   //!
   TBranch        *b_patElectron_phi;   //!
   TBranch        *b_patElectron_energy;   //!
   TBranch        *b_patElectron_charge;   //!
   TBranch        *b_patElectron_isPassVeto;   //!
   TBranch        *b_patElectron_isPassLoose;   //!
   TBranch        *b_patElectron_isPassMedium;   //!
   TBranch        *b_patElectron_isPassTight;   //!
   TBranch        *b_patElectron_isPassHEEPId;   //!
   TBranch        *b_patElectron_d0;   //!
   TBranch        *b_patElectron_dz;   //!
   TBranch        *b_patElectron_expectedMissingInnerHits;   //!
   TBranch        *b_patElectron_passConversionVeto;   //!
   TBranch        *b_patElectron_isoChargedHadrons;   //!
   TBranch        *b_patElectron_isoNeutralHadrons;   //!
   TBranch        *b_patElectron_isoPhotons;   //!
   TBranch        *b_patElectron_isoPU;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_energy;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayModeFinding;   //!
   TBranch        *b_Tau_decayModeFindingNewDMs;   //!
   TBranch        *b_Tau_chargedIsoPtSum;   //!
   TBranch        *b_Tau_neutralIsoPtSum;   //!
   TBranch        *b_Tau_againstMuonTight3;   //!
   TBranch        *b_Tau_againstElectronMVATightMVA5;   //!
   TBranch        *b_Tau_nProngs;   //!
   TBranch        *b_Tau_puCorrPtSum;   //!
   TBranch        *b_Tau_byLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byMediumCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byTightCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byLooseIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byMediumIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byTightIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_againstMuonLoose2;   //!
   TBranch        *b_Tau_againstMuonLoose3;   //!
   TBranch        *b_Tau_againstMuonTight2;   //!
   TBranch        *b_Tau_againstElectronMVALooseMVA5;   //!
   TBranch        *b_Tau_againstElectronMVAMediumMVA5;   //!
   TBranch        *b_Tau_byVLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_leadChargedCandPt;   //!
   TBranch        *b_Tau_leadChargedCandCharge;   //!
   TBranch        *b_Tau_leadChargedCandEta;   //!
   TBranch        *b_Tau_leadChargedCandPhi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_bDiscriminator;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_neutralHadEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_electronEnergy;   //!
   TBranch        *b_Jet_photonEnergy;   //!
   TBranch        *b_UncorrJet_pt;   //!
   TBranch        *b_npuVertices;   //!
   TBranch        *b_trueInteractions;   //!
   TBranch        *b_ootnpuVertices;   //!
   TBranch        *b_npuVerticesp1;   //!
   TBranch        *b_bestVertices;   //!
   TBranch        *b_Met_pt;   //!
   TBranch        *b_Met_sumEt;   //!
   TBranch        *b_Met_phi;   //!
   TBranch        *b_Met_px;   //!
   TBranch        *b_Met_py;   //!
   TBranch        *b_Met_pz;   //!
   TBranch        *b_Gen_Met;   //!
   TBranch        *b_Met_shiftedPtUp;   //!
   TBranch        *b_Met_shiftedPtDown;   //!

};
#endif

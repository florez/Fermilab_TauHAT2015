/////////////////////////////////////////////////////////////////////
/// Developer: Andres Florez, Universidad de los Andes, Colombia. //
///////////////////////////////////////////////////////////////////
#include "BSM_Analysis.h"

int main (int argc, char *argv[])
{

  // The line below creates a root file where the 
  // output histograms will be stored. 
  // The "argv[2]" paramter is the name of the output file
  // which is passed to the code from the console when 
  // running the code. 
  TFile * OutputHistos = new TFile(argv[2], "RECREATE");  
 
  // The four lines below are used to create a set of 
  // sub-directories insed the output root file.
  // each sub-directory will contain histograms after 
  // a set of selections criteria to be defined by the user.
  // The user can change the names of the strings used as 
  // example below: "AfterMuonChargeProduct" and "AfterTau_decayModeFinding".
  // The user can also include new directories. The "nDir" have to match 
  // the number of directories created by the user.
  // The array "theDirectory" should have unique IDs in order to 
  // point to the specific directory: e.g
  // theDirectory[0] = OutputHistos->mkdir("cut_0");
  // theDirectory[1] = OutputHistos->mkdir("cut_1");
  // theDirectory[2] = OutputHistos->mkdir("cut_2");
  // theDirectory[3] = OutputHistos->mkdir("cut_3");

  int nDir = 2;
  TDirectory *theDirectory[nDir];
  theDirectory[0] = OutputHistos->mkdir("AfterMuonChargeProduct");
  theDirectory[1] = OutputHistos->mkdir("AfterTau_decayModeFinding");
 
  // The "argv[1]" argument, is the location+name of the input file 
  // the user will run over. This parameter is also passed to the code 
  // from the console, and it is the first parameter to be passed to the code. 
  BSM_Analysis BSM_Analysis_(OutputHistos, theDirectory, nDir, argv[1]);
}

BSM_Analysis::BSM_Analysis(TFile* theFile, TDirectory *cdDir[], int nDir, char* fname) 
{
  // This function creates a set of histograms 
  // The histogrmas defined in this funtion will be 
  // created and stored inside each of the folders 
  // defined by the user in the "theDirectory" array.
  createHistoMaps(nDir);
  //load PU weights. We are using at the monent the old PU historgams
  // from 8 TeV. Nevertehelss, weight are not used at the moment. 
  // This peace of code is intended to show you how to implement PU 
  // weights in your analysis code in the future.
  TFile file_PUdata("PUdata.root","read");
  TH1F *PUweights = (TH1F*)file_PUdata.Get("analyzeHiMassTau/NVertices_0");
  PUweights->Scale(1/PUweights->Integral());
  TFile file_PUsim("PUsim.root","read");
  TH1F *PUsim = (TH1F*)file_PUsim.Get("analyzeHiMassTau/NVertices_0");
  PUsim->Scale(1/PUsim->Integral());

  PUweights->Divide(PUsim);

  //configure input file
  TFile *f = TFile::Open(fname);
  f->cd("TNT");
  TTree* BOOM = (TTree*)f->Get("TNT/BOOM");

  int nentries = (int) BOOM->GetEntries();
  setBranchAddress(BOOM);

  for (int i = 0; i < nentries; ++i)
    {
      BOOM->GetEntry(i);

       //define global event weight
       double weight =1.;
       weight=PUweights->GetBinContent(PUweights->FindBin(trueInteractions));
       // TLorentz vector to calculate the di-jet invariant mass
       TLorentzVector TagMuon_TL_vec(0., 0., 0., 0.); 
       TLorentzVector ProbeTagMuon_TL_vec(0., 0., 0., 0.);

       double charge_lead = 0.;
       double charge_slead = 0.;
       // This array is used to store the information 
       // of taus that pass or fail a given tau ID 
       // selection criterion.
       int pass_tau_id[nDir] = {0};

       bool pass_trigger = true;


       // For Trigger
       for (int t = 0 ; t < Trigger_decision->size(); t++){
          string theTriggers = Trigger_names->at(t);
          //cout <<theTriggers<<endl;
          string myTrigger   = "HLT_IsoMu24";
          std::size_t found = theTriggers.find(myTrigger);
          if (found!=std::string::npos){
            if (Trigger_decision->at(t) == 1){pass_trigger = true;}
          }
       }       

       int    matched_tau_index = 0;
       // Loop over muons:
       for (int j = 0; j < Muon_pt->size(); j++)
         {
            double DeltaR_matched = 99999.9;

            // select Z -> mumu events
            double lead_muon_pt = Muon_pt->at(0);
            if ( (Muon_pt->size() == 2) && (abs(Muon_eta->at(j)) < 2.4) && (Muon_pt->at(j) > 10.0) ){
              if ((Muon_tight->at(j) == 1) && (Muon_isoSum->at(j) < 5.0) && (pass_trigger)){
                  TagMuon_TL_vec.SetPtEtaPhiE(Muon_pt->at(j), Muon_eta->at(j), Muon_phi->at(j), Muon_energy->at(j)); 
                  charge_lead = Muon_charge->at(j);
              } else {
                  ProbeTagMuon_TL_vec.SetPtEtaPhiE(Muon_pt->at(j), Muon_eta->at(j), Muon_phi->at(j), Muon_energy->at(j));
                  charge_slead = Muon_charge->at(j);
                  // Now we loop over the taus in the event, if any, and find the 
                  // closes tau ot the probe muon. This information might be useful to 
                  // undertand what fraction of taus faking muons are near by the muon candidate.
                  for (int t = 0; t < Tau_pt->size(); t++)
                    {
                     double delta_eta = Muon_eta->at(j) - Tau_eta->at(t);
                     double delta_phi = Muon_phi->at(j) - Tau_phi->at(t);
                     double DeltaR_muon_tau = sqrt(pow(delta_eta, 2) + pow(delta_phi, 2));

                     if ( Tau_pt->size() != (t+1)){
                       if (DeltaR_muon_tau < DeltaR_matched) {
                          DeltaR_matched  = DeltaR_muon_tau;
                          matched_tau_index = t;
                       } 
                     }else {
                      if (DeltaR_muon_tau < DeltaR_matched) {
                          DeltaR_matched  = DeltaR_muon_tau;
                          matched_tau_index = t;
                       }
                        
                     } 
                      if( Tau_decayModeFinding->at(t) == 1){
                         pass_tau_id[1] = 1;
                      }
                   }
              }
                
            }
         }

        double charge_product = charge_lead*charge_slead;
        bool passMuonPt = false;
        _hmap_events[0]->Fill(0.0);
        if ((ProbeTagMuon_TL_vec.Pt() > 25.) || (TagMuon_TL_vec.Pt() > 25.)){
           passMuonPt = true;
           _hmap_events[0]->Fill(1.0);
        }
        if ((charge_product <= -1) && passMuonPt){
          _hmap_events[0]->Fill(2.0);
          _hmap_diMuon_mass[0]->Fill((TagMuon_TL_vec + ProbeTagMuon_TL_vec).M());
          _hmap_muon_pT[0]->Fill(ProbeTagMuon_TL_vec.Pt());
          _hmap_muon_eta[0]->Fill(ProbeTagMuon_TL_vec.Eta());
          _hmap_diMuon_mass_fail[0]->Fill((TagMuon_TL_vec + ProbeTagMuon_TL_vec).M());
          _hmap_muon_pT_fail[0]->Fill(ProbeTagMuon_TL_vec.Pt());
          _hmap_muon_eta_fail[0]->Fill(ProbeTagMuon_TL_vec.Eta());
          if (Tau_pt->size() != 0) fillHistos_matching(0, matched_tau_index, true);
        }
        for (int i = 1; i < nDir; i++){
          _hmap_events[i]->Fill(0.0);
          if ((charge_product == -1) && (passMuonPt) && (pass_tau_id[i] == 1)){
            _hmap_events[i]->Fill(1.0);
            _hmap_diMuon_mass[i]->Fill((TagMuon_TL_vec + ProbeTagMuon_TL_vec).M());
            _hmap_muon_pT[i]->Fill(ProbeTagMuon_TL_vec.Pt());
            _hmap_muon_eta[i]->Fill(ProbeTagMuon_TL_vec.Eta());
            if (Tau_pt->size() != 0) fillHistos_matching(i, matched_tau_index, true);
          } else if ((charge_product == -1) && (passMuonPt) && (pass_tau_id[i] == 0)) {
            _hmap_events[i]->Fill(2.0);
            _hmap_diMuon_mass_fail[i]->Fill((TagMuon_TL_vec + ProbeTagMuon_TL_vec).M());
            _hmap_muon_pT_fail[i]->Fill(ProbeTagMuon_TL_vec.Pt());
            _hmap_muon_eta_fail[i]->Fill(ProbeTagMuon_TL_vec.Eta());
          }
        }
     }
   theFile->cd();
   for (int d = 0; d < nDir; d++)
     {
       cdDir[d]->cd();
       _hmap_events[d]->Write();
       _hmap_diMuon_mass[d]->Write();
       _hmap_muon_pT[d]->Write();
       _hmap_muon_eta[d]->Write();
       _hmap_diMuon_mass_fail[d]->Write();
       _hmap_muon_pT_fail[d]->Write();
       _hmap_muon_eta_fail[d]->Write();
       _hmap_taus_probe_muon_w_matching_pT[d]->Write();
       _hmap_taus_probe_muon_w_matching_eta[d]->Write();
     } 
   theFile->Close();
}

BSM_Analysis::~BSM_Analysis()
{
  // do anything here that needs to be done at desctruction time
}

void BSM_Analysis::fillHistos_matching(int dir, int tau, bool match){

if(match) {
   _hmap_taus_probe_muon_w_matching_pT[dir]->Fill(Tau_pt->at(tau));
   _hmap_taus_probe_muon_w_matching_eta[dir]->Fill(Tau_eta->at(tau));
}

}


void BSM_Analysis::createHistoMaps (int directories)
{
   for (int i = 0; i < directories; i++)
     {
       // Events histogram
       _hmap_events[i]           = new TH1F("Events",          "Events", 3, 0., 3.);
       // Muon distributions
       _hmap_diMuon_mass[i]      = new TH1F("diMuonMass",      "m_{#mu, #mu}", 300., 0., 300.); 
       _hmap_muon_pT[i]          = new TH1F("muon_pT",         "#mu p_{T}",    300, 0., 300.);
       _hmap_muon_eta[i]         = new TH1F("muon_eta",        "#mu #eta",     50, -2.5, 2.5);
       _hmap_diMuon_mass_fail[i] = new TH1F("diMuonMass_fail", "m_{#mu, #mu}", 300., 0., 300.);
       _hmap_muon_pT_fail[i]     = new TH1F("muon_pT_fail",    "#mu p_{T}",    300, 0., 300.);
       _hmap_muon_eta_fail[i]    = new TH1F("muon_eta_fail",   "#mu #eta",     50, -2.5, 2.5);    
       // tau distributions
       _hmap_taus_probe_muon_w_matching_pT[i]      = new TH1F("taus_probe_muon_w_matching_pT",   "#tau p_{T}", 300, 0., 300.);
       _hmap_taus_probe_muon_w_matching_eta[i]     = new TH1F("taus_probe_muon_w_matching_eta",  "#tau #eta",  50, -2.5, 2.5);
     }
}

void BSM_Analysis::setBranchAddress(TTree* BOOM)
{

   // Set object pointer
   
   Muon_pt = 0;
   Muon_eta = 0;
   Muon_phi = 0;
   Muon_p = 0;
   Muon_energy = 0;
   Muon_charge = 0;
   Muon_tight = 0;
   Muon_soft = 0;
   Muon_pf = 0;
   Muon_isoCharged = 0;
   Muon_isoSum = 0;
   Muon_isoCharParPt = 0;
   Muon_chi2 = 0;
   Muon_validHits = 0;
   Muon_validHitsInner = 0;
   Muon_matchedStat = 0;
   Muon_dxy = 0;
   Muon_TLayers = 0;
   Muon_dz = 0;
   Muon_isoNeutralHadron = 0;
   Muon_isoPhoton = 0;
   Muon_isoPU = 0;
   Tau_eta = 0;
   Tau_phi = 0;
   Tau_pt = 0;
   Tau_energy = 0;
   Tau_charge = 0;
   Tau_decayModeFinding = 0;
   Tau_decayModeFindingNewDMs = 0;
   Tau_chargedIsoPtSum = 0;
   Tau_neutralIsoPtSum = 0;
   Tau_againstMuonTight3 = 0;
   Tau_againstElectronMVATightMVA5 = 0;
   Tau_nProngs = 0;
   Tau_puCorrPtSum = 0;
   Tau_byLooseCombinedIsolationDeltaBetaCorr = 0;
   Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byMediumCombinedIsolationDeltaBetaCorr = 0;
   Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byTightCombinedIsolationDeltaBetaCorr = 0;
   Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byLooseIsolationMVA3newDMwLT = 0;
   Tau_byLooseIsolationMVA3newDMwoLT = 0;
   Tau_byLooseIsolationMva3oldDMwLT = 0;
   Tau_byLooseIsolationMVA3oldDMwoLT = 0;
   Tau_byMediumIsolationMVA3newDMwLT = 0;
   Tau_byMediumIsolationMVA3newDMwoLT = 0;
   Tau_byMediumIsolationMva3oldDMwLT = 0;
   Tau_byMediumIsolationMVA3oldDMwoLT = 0;
   Tau_byTightIsolationMVA3newDMwLT = 0;
   Tau_byTightIsolationMVA3newDMwoLT = 0;
   Tau_byTightIsolationMva3oldDMwLT = 0;
   Tau_byTightIsolationMVA3oldDMwoLT = 0;
   Tau_againstMuonLoose2 = 0;
   Tau_againstMuonLoose3 = 0;
   Tau_againstMuonTight2 = 0;
   Tau_againstElectronMVALooseMVA5 = 0;
   Tau_againstElectronMVAMediumMVA5 = 0;
   Tau_byVLooseCombinedIsolationDeltaBetaCorr = 0;
   Tau_leadChargedCandPt = 0;
   Tau_leadChargedCandCharge = 0;
   Tau_leadChargedCandEta = 0;
   Tau_leadChargedCandPhi = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_energy = 0;
   Jet_bDiscriminator = 0;
   Jet_mass = 0;
   Jet_neutralHadEnergy = 0;
   Jet_neutralEmEmEnergy = 0;
   Jet_chargedHadronEnergy = 0;
   Jet_chargedEmEnergy = 0;
   Jet_muonEnergy = 0;
   Jet_electronEnergy = 0;
   Jet_photonEnergy = 0;
   UncorrJet_pt = 0;
   Trigger_names = 0;
   Trigger_decision = 0;
   // Set branch addresses and branch pointers
   if(!BOOM) return;
   BOOM->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
   BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision, &b_Trigger_decision);
   BOOM->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   BOOM->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   BOOM->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   BOOM->SetBranchAddress("Muon_p", &Muon_p, &b_Muon_p);
   BOOM->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
   BOOM->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   BOOM->SetBranchAddress("Muon_tight", &Muon_tight, &b_Muon_tight);
   BOOM->SetBranchAddress("Muon_soft", &Muon_soft, &b_Muon_soft);
   BOOM->SetBranchAddress("Muon_pf", &Muon_pf, &b_Muon_pf);
   BOOM->SetBranchAddress("Muon_isoCharged", &Muon_isoCharged, &b_Muon_isoCharged);
   BOOM->SetBranchAddress("Muon_isoSum", &Muon_isoSum, &b_Muon_isoSum);
   BOOM->SetBranchAddress("Muon_isoCharParPt", &Muon_isoCharParPt, &b_Muon_isoCharParPt);
   BOOM->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
   BOOM->SetBranchAddress("Muon_validHits", &Muon_validHits, &b_Muon_validHits);
   BOOM->SetBranchAddress("Muon_validHitsInner", &Muon_validHitsInner, &b_Muon_validHitsInner);
   BOOM->SetBranchAddress("Muon_matchedStat", &Muon_matchedStat, &b_Muon_matchedStat);
   BOOM->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy);
   BOOM->SetBranchAddress("Muon_TLayers", &Muon_TLayers, &b_Muon_TLayers);
   BOOM->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
   BOOM->SetBranchAddress("Muon_isoNeutralHadron", &Muon_isoNeutralHadron, &b_Muon_isoNeutralHadron);
   BOOM->SetBranchAddress("Muon_isoPhoton", &Muon_isoPhoton, &b_Muon_isoPhoton);
   BOOM->SetBranchAddress("Muon_isoPU", &Muon_isoPU, &b_Muon_isoPU);
   BOOM->SetBranchAddress("Tau_eta", &Tau_eta, &b_Tau_eta);
   BOOM->SetBranchAddress("Tau_phi", &Tau_phi, &b_Tau_phi);
   BOOM->SetBranchAddress("Tau_pt", &Tau_pt, &b_Tau_pt);
   BOOM->SetBranchAddress("Tau_energy", &Tau_energy, &b_Tau_energy);
   BOOM->SetBranchAddress("Tau_charge", &Tau_charge, &b_Tau_charge);
   BOOM->SetBranchAddress("Tau_decayModeFinding", &Tau_decayModeFinding, &b_Tau_decayModeFinding);
   BOOM->SetBranchAddress("Tau_decayModeFindingNewDMs", &Tau_decayModeFindingNewDMs, &b_Tau_decayModeFindingNewDMs);
   BOOM->SetBranchAddress("Tau_chargedIsoPtSum", &Tau_chargedIsoPtSum, &b_Tau_chargedIsoPtSum);
   BOOM->SetBranchAddress("Tau_neutralIsoPtSum", &Tau_neutralIsoPtSum, &b_Tau_neutralIsoPtSum);
   BOOM->SetBranchAddress("Tau_againstMuonTight3", &Tau_againstMuonTight3, &b_Tau_againstMuonTight3);
   BOOM->SetBranchAddress("Tau_againstElectronMVATightMVA5", &Tau_againstElectronMVATightMVA5, &b_Tau_againstElectronMVATightMVA5);
   BOOM->SetBranchAddress("Tau_nProngs", &Tau_nProngs, &b_Tau_nProngs);
   BOOM->SetBranchAddress("Tau_puCorrPtSum", &Tau_puCorrPtSum, &b_Tau_puCorrPtSum);
   BOOM->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr", &Tau_byLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   BOOM->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr", &Tau_byMediumCombinedIsolationDeltaBetaCorr, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   BOOM->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr", &Tau_byTightCombinedIsolationDeltaBetaCorr, &b_Tau_byTightCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &Tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwLT", &Tau_byLooseIsolationMVA3newDMwLT, &b_Tau_byLooseIsolationMVA3newDMwLT);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwoLT", &Tau_byLooseIsolationMVA3newDMwoLT, &b_Tau_byLooseIsolationMVA3newDMwoLT);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMva3oldDMwLT", &Tau_byLooseIsolationMva3oldDMwLT, &b_Tau_byLooseIsolationMva3oldDMwLT);
   BOOM->SetBranchAddress("Tau_byLooseIsolationMVA3oldDMwoLT", &Tau_byLooseIsolationMVA3oldDMwoLT, &b_Tau_byLooseIsolationMVA3oldDMwoLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwLT", &Tau_byMediumIsolationMVA3newDMwLT, &b_Tau_byMediumIsolationMVA3newDMwLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwoLT", &Tau_byMediumIsolationMVA3newDMwoLT, &b_Tau_byMediumIsolationMVA3newDMwoLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMva3oldDMwLT", &Tau_byMediumIsolationMva3oldDMwLT, &b_Tau_byMediumIsolationMva3oldDMwLT);
   BOOM->SetBranchAddress("Tau_byMediumIsolationMVA3oldDMwoLT", &Tau_byMediumIsolationMVA3oldDMwoLT, &b_Tau_byMediumIsolationMVA3oldDMwoLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMVA3newDMwLT", &Tau_byTightIsolationMVA3newDMwLT, &b_Tau_byTightIsolationMVA3newDMwLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMVA3newDMwoLT", &Tau_byTightIsolationMVA3newDMwoLT, &b_Tau_byTightIsolationMVA3newDMwoLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMva3oldDMwLT", &Tau_byTightIsolationMva3oldDMwLT, &b_Tau_byTightIsolationMva3oldDMwLT);
   BOOM->SetBranchAddress("Tau_byTightIsolationMVA3oldDMwoLT", &Tau_byTightIsolationMVA3oldDMwoLT, &b_Tau_byTightIsolationMVA3oldDMwoLT);
   BOOM->SetBranchAddress("Tau_againstMuonLoose2", &Tau_againstMuonLoose2, &b_Tau_againstMuonLoose2);
   BOOM->SetBranchAddress("Tau_againstMuonLoose3", &Tau_againstMuonLoose3, &b_Tau_againstMuonLoose3);
   BOOM->SetBranchAddress("Tau_againstMuonTight2", &Tau_againstMuonTight2, &b_Tau_againstMuonTight2);
   BOOM->SetBranchAddress("Tau_againstElectronMVALooseMVA5", &Tau_againstElectronMVALooseMVA5, &b_Tau_againstElectronMVALooseMVA5);
   BOOM->SetBranchAddress("Tau_againstElectronMVAMediumMVA5", &Tau_againstElectronMVAMediumMVA5, &b_Tau_againstElectronMVAMediumMVA5);
   BOOM->SetBranchAddress("Tau_byVLooseCombinedIsolationDeltaBetaCorr", &Tau_byVLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byVLooseCombinedIsolationDeltaBetaCorr);
   BOOM->SetBranchAddress("Tau_leadChargedCandPt", &Tau_leadChargedCandPt, &b_Tau_leadChargedCandPt);
   BOOM->SetBranchAddress("Tau_leadChargedCandCharge", &Tau_leadChargedCandCharge, &b_Tau_leadChargedCandCharge);
   BOOM->SetBranchAddress("Tau_leadChargedCandEta", &Tau_leadChargedCandEta, &b_Tau_leadChargedCandEta);
   BOOM->SetBranchAddress("Tau_leadChargedCandPhi", &Tau_leadChargedCandPhi, &b_Tau_leadChargedCandPhi);
   BOOM->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   BOOM->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   BOOM->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   BOOM->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   BOOM->SetBranchAddress("Jet_bDiscriminator", &Jet_bDiscriminator, &b_Jet_bDiscriminator);
   BOOM->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   BOOM->SetBranchAddress("Jet_neutralHadEnergy", &Jet_neutralHadEnergy, &b_Jet_neutralHadEnergy);
   BOOM->SetBranchAddress("Jet_neutralEmEmEnergy", &Jet_neutralEmEmEnergy, &b_Jet_neutralEmEmEnergy);
   BOOM->SetBranchAddress("Jet_chargedHadronEnergy", &Jet_chargedHadronEnergy, &b_Jet_chargedHadronEnergy);
   BOOM->SetBranchAddress("Jet_chargedEmEnergy", &Jet_chargedEmEnergy, &b_Jet_chargedEmEnergy);
   BOOM->SetBranchAddress("Jet_muonEnergy", &Jet_muonEnergy, &b_Jet_muonEnergy);
   BOOM->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
   BOOM->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
   BOOM->SetBranchAddress("UncorrJet_pt", &UncorrJet_pt, &b_UncorrJet_pt);
   BOOM->SetBranchAddress("npuVertices", &npuVertices, &b_npuVertices);
   BOOM->SetBranchAddress("trueInteractions", &trueInteractions, &b_trueInteractions);
   BOOM->SetBranchAddress("ootnpuVertices", &ootnpuVertices, &b_ootnpuVertices);
   BOOM->SetBranchAddress("npuVerticesp1", &npuVerticesp1, &b_npuVerticesp1);
   BOOM->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
   BOOM->SetBranchAddress("Met_pt", &Met_pt, &b_Met_pt);
   BOOM->SetBranchAddress("Met_sumEt", &Met_sumEt, &b_Met_sumEt);
   BOOM->SetBranchAddress("Met_phi", &Met_phi, &b_Met_phi);
   BOOM->SetBranchAddress("Met_px", &Met_px, &b_Met_px);
   BOOM->SetBranchAddress("Met_py", &Met_py, &b_Met_py);
   BOOM->SetBranchAddress("Met_pz", &Met_pz, &b_Met_pz);
   BOOM->SetBranchAddress("Gen_Met", &Gen_Met, &b_Gen_Met);
   BOOM->SetBranchAddress("Met_shiftedPtUp", &Met_shiftedPtUp, &b_Met_shiftedPtUp);
   BOOM->SetBranchAddress("Met_shiftedPtDown", &Met_shiftedPtDown, &b_Met_shiftedPtDown);
};

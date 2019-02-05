#define ZPeak_cxx
#include "ZPeak.h"
#include <TDirectory.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>

void ZPeak::Loop(TString dataset)
{
//   In a ROOT session, you can do:
//      root> .L ZPeak.C
//      root> ZPeak t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   bool isRealData = dataset.Contains("DoubleMuon") | dataset.Contains("DoubleEG");
   std::cout << "Processing " << dataset << ", real=" << isRealData << ", nEntries=" << nentries << std::endl;


   std::unique_ptr<TH2F> eleCorr;
   std::unique_ptr<TH2F> muCorr;
   std::unique_ptr<TH1D> puCorr;

   if ( !isRealData ) {
    TFile* fele = TFile::Open("correction_files/scalefactors_80x_egpog_37ifb.root");
    eleCorr.reset((TH2F*) fele->Get("scalefactors_Tight_Electron"));
    eleCorr->SetDirectory(0);
    delete fele;

    TFile* fmu = TFile::Open("correction_files/muon_scalefactors_37ifb.root");
    muCorr.reset((TH2F*) fmu->Get("scalefactors_Iso_MuonTightId"));
    muCorr->SetDirectory(0);
    delete fmu;

    TFile* fpu = TFile::Open("correction_files/puReweight0p38fb.root");
    puCorr.reset((TH1D*) fpu->Get(dataset));
    puCorr->SetDirectory(0);
    delete fpu;
   }


   TFile* fOut = new TFile(dataset+TString(".root"), "recreate");
   TDirectory* dOut = fOut->mkdir(dataset);
   auto mass_ee = new TH1D("mass_ee", "mass;M_{ll} [GeV];Counts", 120, 0., 120.);
   auto mass_em = new TH1D("mass_em", "mass;M_{ll} [GeV];Counts", 120, 0., 120.);
   auto mass_mm = new TH1D("mass_mm", "mass;M_{ll} [GeV];Counts", 120, 0., 120.);
   auto lepPt_ee = new TH2D("lepPt_ee", "LepPt;Leading p_{T}^{l} [GeV];Trailing p_{T}^{l} [GeV];Counts", 50, 0., 500., 50, 0., 500.);
   auto lepPt_em = new TH2D("lepPt_em", "LepPt;Leading p_{T}^{l} [GeV];Trailing p_{T}^{l} [GeV];Counts", 50, 0., 500., 50, 0., 500.);
   auto lepPt_mm = new TH2D("lepPt_mm", "LepPt;Leading p_{T}^{l} [GeV];Trailing p_{T}^{l} [GeV];Counts", 50, 0., 500., 50, 0., 500.);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      b_run->GetEntry(ientry);
      b_luminosityBlock->GetEntry(ientry);
      if( !isRealData ) {
        b_LHEWeight_originalXWGTUP->GetEntry(ientry);
        b_Pileup_nTrueInt->GetEntry(ientry);
      }
      b_nElectron->GetEntry(ientry);
      b_Electron_pt->GetEntry(ientry);
      b_Electron_eta->GetEntry(ientry);
      b_Electron_phi->GetEntry(ientry);
      b_Electron_mass->GetEntry(ientry);
      b_Electron_cutBased->GetEntry(ientry);
      b_Electron_pdgId->GetEntry(ientry);
      b_Electron_pfRelIso03_all->GetEntry(ientry);
      b_nMuon->GetEntry(ientry);
      b_Muon_pt->GetEntry(ientry);
      b_Muon_eta->GetEntry(ientry);
      b_Muon_phi->GetEntry(ientry);
      b_Muon_mass->GetEntry(ientry);
      b_Muon_tightId->GetEntry(ientry);
      b_Muon_pdgId->GetEntry(ientry);
      b_Muon_pfRelIso04_all->GetEntry(ientry);
      b_HLT_Ele32_WPTight_Gsf_L1DoubleEG->GetEntry(ientry);
      b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ->GetEntry(ientry);
      b_HLT_IsoMu24->GetEntry(ientry);
      b_HLT_IsoMu27->GetEntry(ientry);
      b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ->GetEntry(ientry);

      double weight = 1.;
      if ( !isRealData ) {
        weight *= LHEWeight_originalXWGTUP;
        weight *= puCorr->GetBinContent(puCorr->FindBin(Pileup_nTrueInt));
      }

      bool good_trigger{false};
      if ( (dataset.Contains("DoubleEG") or !isRealData) and (HLT_Ele32_WPTight_Gsf_L1DoubleEG or HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) ) {
        good_trigger = true;
      }
      if ( (dataset.Contains("DoubleMuon") or !isRealData) and (HLT_IsoMu24 or HLT_IsoMu27 or HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ) ) {
        good_trigger = true;
      }

      TLorentzVector lep1, lep2;
      int lep1_id{0}, lep2_id{0};

      bool tooManyLep{false};

      for(size_t iEl=0; iEl<nElectron; ++iEl) {
        if (Electron_pt[iEl] > 20. && std::abs(Electron_eta[iEl]) < 2.5 && Electron_cutBased[iEl] >= 4) {
          if (lep1_id==0) {
            lep1.SetPtEtaPhiM(Electron_pt[iEl], Electron_eta[iEl], Electron_phi[iEl], Electron_mass[iEl]);
            lep1_id = Electron_pdgId[iEl];
            if (!isRealData) weight *= eleCorr->GetBinContent(eleCorr->FindBin(Electron_eta[iEl], Electron_pt[iEl]));
          } else if (lep2_id==0) {
            lep2.SetPtEtaPhiM(Electron_pt[iEl], Electron_eta[iEl], Electron_phi[iEl], Electron_mass[iEl]);
            lep2_id = Electron_pdgId[iEl];
            if (!isRealData) weight *= eleCorr->GetBinContent(eleCorr->FindBin(Electron_eta[iEl], Electron_pt[iEl]));
          } else {
            tooManyLep = true;
          }
        }
      }

      for(size_t iMu=0; iMu<nMuon; ++iMu) {
        if (Muon_pt[iMu] > 20. && std::abs(Muon_eta[iMu]) < 2.4 && Muon_tightId[iMu] > 0) {
          if (lep1_id==0) {
            lep1.SetPtEtaPhiM(Muon_pt[iMu], Muon_eta[iMu], Muon_phi[iMu], Muon_mass[iMu]);
            lep1_id = Muon_pdgId[iMu];
            if (!isRealData) weight *= muCorr->GetBinContent(muCorr->FindBin(std::abs(Muon_eta[iMu]), Muon_pt[iMu]));
          } else if (lep2_id==0) {
            lep2.SetPtEtaPhiM(Muon_pt[iMu], Muon_eta[iMu], Muon_phi[iMu], Muon_mass[iMu]);
            lep2_id = Muon_pdgId[iMu];
            if (!isRealData) weight *= muCorr->GetBinContent(muCorr->FindBin(std::abs(Muon_eta[iMu]), Muon_pt[iMu]));
          } else {
            tooManyLep = true;
          }
        }
      }

      if (tooManyLep) continue;
      if ( lep2.Pt() > lep1.Pt() ) {
        std::swap(lep1, lep2);
        std::swap(lep1_id, lep2_id);
      }
      TLorentzVector zcand = lep1 + lep2;

      if ( lep1_id*lep2_id == -11*11 && lep1.Pt() > 25. ) {
        mass_ee->Fill(zcand.M(), weight);
        lepPt_ee->Fill(lep1.Pt(), lep2.Pt(), weight);
      }
      else if ( lep1_id*lep2_id == -11*13 ) {
        mass_em->Fill(zcand.M(), weight);
        lepPt_em->Fill(lep1.Pt(), lep2.Pt(), weight);
      }
      else if ( lep1_id*lep2_id == -13*13 ) {
        mass_mm->Fill(zcand.M(), weight);
        lepPt_mm->Fill(lep1.Pt(), lep2.Pt(), weight);
      }
   }
   std::cout << "Done " << dataset << std::endl;

   dOut->cd();
   mass_ee->Write();
   mass_em->Write();
   mass_mm->Write();
   lepPt_ee->Write();
   lepPt_em->Write();
   lepPt_mm->Write();
   fOut->Close();
}

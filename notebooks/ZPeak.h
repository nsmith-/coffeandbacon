//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  4 13:38:53 2019 by ROOT version 6.14/06
// from TChain Events/
//////////////////////////////////////////////////////////

#ifndef ZPeak_h
#define ZPeak_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

// Header file for the classes stored in the TTree if any.

class ZPeak {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   UInt_t          nElectron;
   Float_t         Electron_eta[8];   //[nElectron]
   Float_t         Electron_mass[8];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[8];   //[nElectron]
   Float_t         Electron_phi[8];   //[nElectron]
   Float_t         Electron_pt[8];   //[nElectron]
   Int_t           Electron_cutBased[8];   //[nElectron]
   Int_t           Electron_pdgId[8];   //[nElectron]
   Float_t         LHEWeight_originalXWGTUP;
   UInt_t          nMuon;
   Float_t         Muon_eta[7];   //[nMuon]
   Float_t         Muon_mass[7];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[7];   //[nMuon]
   Float_t         Muon_phi[7];   //[nMuon]
   Float_t         Muon_pt[7];   //[nMuon]
   Int_t           Muon_pdgId[7];   //[nMuon]
   Bool_t          Muon_tightId[7];   //[nMuon]
   Float_t         Pileup_nTrueInt;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_LHEWeight_originalXWGTUP;   //!
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_L1DoubleEG;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!

   ZPeak(TTree *tree=0);
   virtual ~ZPeak();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString dataset);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   bool isRealData_;
   TString dataset_;
   std::unique_ptr<TH2F> eleCorr_;
   std::unique_ptr<TH2F> muCorr_;
   std::unique_ptr<TH1D> puCorr_;

   TH1D* mass_ee;
   TH1D* mass_em;
   TH1D* mass_mm;
   TH2D* lepPt_ee;
   TH2D* lepPt_em;
   TH2D* lepPt_mm;

   void actuallyProcess();
};

#endif

#ifdef ZPeak_cxx
ZPeak::ZPeak(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("Events","");
      chain->Add("data/EA60F185-50A3-E811-8AC2-008CFA1113FC.root/Events");
      chain->Add("data/8EB617A3-3842-E811-BF83-001E67792558.root/Events");
      chain->Add("data/6CC2B88E-C444-E811-87B6-1866DAEA7E28.root/Events");
      chain->Add("data/142EA893-BC44-E811-95BB-FA163ECE9AA6.root/Events");
      chain->Add("data/B405DFEE-BC44-E811-A346-FA163EA1976B.root/Events");
      chain->Add("data/FC4BF610-3686-E811-A051-1866DA7F967F.root/Events");
      chain->Add("data/C2B508AE-4042-E811-B174-001E67792702.root/Events");
      chain->Add("data/8254B705-BB44-E811-A808-FA163E95DCD9.root/Events");
      chain->Add("data/38A88160-B6AB-E811-A835-D4AE526A0CEC.root/Events");
      chain->Add("data/382E4666-50A3-E811-855C-5065F382A241.root/Events");
      chain->Add("data/1AFC31AA-BB44-E811-B270-FA163E46AE2E.root/Events");
      chain->Add("data/1CEF9BDB-BC44-E811-996F-801844DF001C.root/Events");
      chain->Add("data/7C2B491C-BCAB-E811-88A2-002590DE6E5E.root/Events");
      chain->Add("data/0AF48D09-3686-E811-8187-00000086FE80.root/Events");
      tree = chain;
   }
   Init(tree);
}

ZPeak::~ZPeak()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZPeak::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZPeak::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZPeak::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);

   Notify();
}

Bool_t ZPeak::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZPeak::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef ZPeak_cxx

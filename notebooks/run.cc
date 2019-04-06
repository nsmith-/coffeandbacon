void run(void) {
  // Run using
  // $ time root -l -q -b -e ".L ZPeak.C+" -e ".x run.cc"

  TChain* chain;
  ZPeak* processor;

  chain = new TChain("Events");
  // chain->Add("data/nano_dy_lzma_combined.root");
  chain->Add("data/nano_dy_lz4_combined.root");
  processor = new ZPeak(chain);
  processor->Loop("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8");
  delete processor;
  delete chain;

}

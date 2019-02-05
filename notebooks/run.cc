void run(void) {

  TChain* chain;
  ZPeak* processor;

  chain = new TChain("Events");
  chain->Add("data/7C2B491C-BCAB-E811-88A2-002590DE6E5E.root");
  chain->Add("data/38A88160-B6AB-E811-A835-D4AE526A0CEC.root");
  processor = new ZPeak(chain);
  processor->Loop("DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8");
  delete processor;
  delete chain;

  chain = new TChain("Events");
  chain->Add("data/8EB617A3-3842-E811-BF83-001E67792558.root");
  chain->Add("data/C2B508AE-4042-E811-B174-001E67792702.root");
  processor = new ZPeak(chain);
  processor->Loop("DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8");
  delete processor;
  delete chain;

  chain = new TChain("Events");
  chain->Add("data/0AF48D09-3686-E811-8187-00000086FE80.root");
  chain->Add("data/FC4BF610-3686-E811-A051-1866DA7F967F.root");
  chain->Add("data/382E4666-50A3-E811-855C-5065F382A241.root");
  chain->Add("data/EA60F185-50A3-E811-8AC2-008CFA1113FC.root");
  processor = new ZPeak(chain);
  processor->Loop("TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8");
  delete processor;
  delete chain;

  chain = new TChain("Events");
  chain->Add("data/8254B705-BB44-E811-A808-FA163E95DCD9.root");
  chain->Add("data/1AFC31AA-BB44-E811-B270-FA163E46AE2E.root");
  chain->Add("data/142EA893-BC44-E811-95BB-FA163ECE9AA6.root");
  processor = new ZPeak(chain);
  processor->Loop("DoubleEG");
  delete processor;
  delete chain;

  chain = new TChain("Events");
  chain->Add("data/1CEF9BDB-BC44-E811-996F-801844DF001C.root");
  chain->Add("data/6CC2B88E-C444-E811-87B6-1866DAEA7E28.root");
  chain->Add("data/B405DFEE-BC44-E811-A346-FA163EA1976B.root");
  processor = new ZPeak(chain);
  processor->Loop("DoubleMuon");
  delete processor;
  delete chain;

}

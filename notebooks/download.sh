mkdir -p data

xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/7C2B491C-BCAB-E811-88A2-002590DE6E5E.root ./data/
xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/38A88160-B6AB-E811-A835-D4AE526A0CEC.root ./data/
xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/40000/8EB617A3-3842-E811-BF83-001E67792558.root ./data/
xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/40000/C2B508AE-4042-E811-B174-001E67792702.root ./data/
xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0AF48D09-3686-E811-8187-00000086FE80.root ./data/
xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/FC4BF610-3686-E811-A051-1866DA7F967F.root ./data/
xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/382E4666-50A3-E811-855C-5065F382A241.root ./data/
xrdcp root://cmseos.fnal.gov//store/mc/RunIIFall17NanoAOD/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/EA60F185-50A3-E811-8AC2-008CFA1113FC.root ./data/
xrdcp root://cmseos.fnal.gov//store/data/Run2017B/DoubleEG/NANOAOD/31Mar2018-v1/70000/8254B705-BB44-E811-A808-FA163E95DCD9.root ./data/
xrdcp root://cmseos.fnal.gov//store/data/Run2017B/DoubleEG/NANOAOD/31Mar2018-v1/70000/1AFC31AA-BB44-E811-B270-FA163E46AE2E.root ./data/
xrdcp root://cmseos.fnal.gov//store/data/Run2017B/DoubleEG/NANOAOD/31Mar2018-v1/70000/142EA893-BC44-E811-95BB-FA163ECE9AA6.root ./data/
xrdcp root://cmseos.fnal.gov//store/data/Run2017B/DoubleMuon/NANOAOD/31Mar2018-v1/70000/1CEF9BDB-BC44-E811-996F-801844DF001C.root ./data/
xrdcp root://cmseos.fnal.gov//store/data/Run2017B/DoubleMuon/NANOAOD/31Mar2018-v1/70000/6CC2B88E-C444-E811-87B6-1866DAEA7E28.root ./data/
xrdcp root://cmseos.fnal.gov//store/data/Run2017B/DoubleMuon/NANOAOD/31Mar2018-v1/70000/B405DFEE-BC44-E811-A346-FA163EA1976B.root ./data/

# {run: [[lumi#, int. lumi, sigma inst. lumi / bx, inst. lumi / bx], ...], ...}
curl https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/PileUp/pileup_latest.txt > metadata/pileup_2017.json

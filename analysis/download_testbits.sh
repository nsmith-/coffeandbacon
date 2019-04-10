#!/usr/bin/env bash

mkdir -p data
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/skim/TTToHadronic_TuneCP5_13TeV_powheg_pythia8_0.root data/
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.02/skim/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_0.root data/
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/skim/JetHTRun2017F_17Nov2017_v1_24.rootnodupl.root data/
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/skim/SingleMuonRun2017B_17Nov2017_v1_2.root data/

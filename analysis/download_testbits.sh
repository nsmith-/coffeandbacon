#!/usr/bin/env bash

mkdir -p data

# skimmed bits
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/skim/TTToHadronic_TuneCP5_13TeV_powheg_pythia8_0.root data/
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.02/skim/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_0.root data/
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/skim/JetHTRun2017F_17Nov2017_v1_24.rootnodupl.root data/
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/skim/SingleMuonRun2017B_17Nov2017_v1_2.root data/

# unskimmed bits
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/GluGluHToBB_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8/Output_job0_job0_file0to59.root data/GluGluHToBB_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8.Output_job0_job0_file0to59.root
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/GluGluHToCC_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8/Output_job15_job15_file900to959.root data/GluGluHToCC_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8.Output_job15_job15_file900to959.root
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/Output_job116_job116_file1160to1169.root data/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8.Output_job116_job116_file1160to1169.root
xrdcp root://cmseos.fnal.gov//eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.01/TTToHadronic_TuneCP5_13TeV_powheg_pythia8/Output_job101_job101_file4040to4079.root data/TTToHadronic_TuneCP5_13TeV_powheg_pythia8.Output_job101_job101_file4040to4079.root

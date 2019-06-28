import ROOT
ROOT.gROOT.SetBatch(True)
import json
from tqdm import tqdm

with open('files.json') as fin:
        files = json.load(fin)

redir = 'root://cmsxrootd.fnal.gov/'
samples = {k: [redir+f for f in v] for k,v in files.items()}
# `hadd -O -fk404 -experimental-io-features kGenerateOffsetMap nano_dy_lzma_combined.root`

for sample in samples:
    if sample == "GluGluHToBB_M125_LHEHpT_250-Inf_13TeV_amcatnloFXFX_pythia8":
        continue
    for i, file in enumerate(tqdm(samples[sample], desc=sample)):
        fin = ROOT.TFile.Open(file)
        events = fin.Get("Events")
        runs = fin.Get("Runs")
        fout = ROOT.TFile.Open("data/%s.%d.root" % (sample, i), "recreate", "", 404)
        fout.cd()
        tev = events.CopyTree("(GenPart_pdgId==25) && (GenPart_statusFlags&128) && (GenPart_pt>250.)")
        trun = runs.CopyTree("")
        tev.Write()
        trun.Write()
        fout.Close()


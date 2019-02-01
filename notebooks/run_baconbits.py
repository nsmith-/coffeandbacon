#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time
import numexpr

import uproot
import numpy as np
from fnal_column_analysis_tools import hist, lookup_tools

with open("metadata/datadef.json") as fin:
    datadef = json.load(fin)

extractor = lookup_tools.extractor()
extractor.add_weight_sets(["* * correction_files/n2ddt_transform_2017MC.root"])
extractor.add_weight_sets(["* * correction_files/TriggerEfficiencies_Run2017_noPS.root"])
extractor.finalize()
evaluator = extractor.make_evaluator()
n2ddt_rho_pt = evaluator[b"Rho2D"]

gpar = np.array([1.00626, -1.06161, 0.0799900, 1.20454])
cpar = np.array([1.09302, -0.000150068, 3.44866e-07, -2.68100e-10, 8.67440e-14, -1.00114e-17])
fpar = np.array([1.27212, -0.000571640, 8.37289e-07, -5.20433e-10, 1.45375e-13, -1.50389e-17])

def msd_weight(pt, eta):
    genw = gpar[0] + gpar[1]*np.power(pt*gpar[2], -gpar[3])
    ptpow = np.power.outer(pt, np.arange(cpar.size))
    cenweight = np.dot(ptpow, cpar)
    forweight = np.dot(ptpow, fpar)
    weight = np.where(np.abs(eta)<1.3, cenweight, forweight)
    return weight


with gzip.open("pileup_mc.pkl.gz", "rb") as fin:
    pileup_corr = pickle.load(fin)

with uproot.open("correction_files/pileup_Cert_294927-306462_13TeV_PromptReco_Collisions17_withVar.root") as fin:
    data_pu = fin["pileup"].values
    for k in pileup_corr.keys():
        mc_pu = pileup_corr[k]
        corr = (data_pu / np.maximum(mc_pu, 1)) / (data_pu.sum() / mc_pu.sum())
        corr[mc_pu==0.] = 1.
        pileup_corr[k] = lookup_tools.dense_lookup.dense_lookup(corr, fin["pileup"].edges)


with uproot.open("correction_files/TriggerEfficiencies_Run2017_noPS.root") as fin:
    denom = fin["h_runBtoF_pass_Mu50"]
    num = fin["h_runBtoF_pass_Main"]
    eff = num.values/np.maximum(denom.values, 1)
    msd_bins, pt_bins = num.edges
    # Cut pt < 200
    pt_bins = pt_bins[8:]
    # ROOT is insane.. the array shape is [y,x]
    eff = eff[8:,:]
    jetTriggerEff = lookup_tools.dense_lookup.dense_lookup(eff, (msd_bins, pt_bins))


# [pb]
dataset_xs = {k: v['xs'] for k,v in datadef.items()}
lumi = 1000.  # [1/pb]

dataset = hist.Cat("dataset", "Primary dataset")

gencat = hist.Bin("AK8Puppijet0_isHadronicV", "Matched", [0,1,2,3,9,10,11])
# one can relabel intervals, although process mapping obviates this
titles = ["QCD", "V(light) matched", "V(c) matched", "V(b) matched", "Top W(ud)+b", "Top W(cs)+b"]
for i,v in enumerate(gencat.identifiers()):
    setattr(v, 'label', titles[i])

jetpt = hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", [450, 500, 550, 600, 675, 800, 1000])
jetpt_coarse = hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", [450, 800])
jetmass = hist.Bin("AK8Puppijet0_msd", "Jet $m_{sd}$", 23, 40, 201)
jetmass_coarse = hist.Bin("AK8Puppijet0_msd", "Jet $m_{sd}$", [40, 100, 140, 200])
jetrho = hist.Bin("jetrho", r"Jet $\rho$", 13, -6, -2.1)
doubleb = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", 20, 0., 1)
doublec = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", 20, 0., 1.)
doublecvb = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", 20, 0., 1.)
doubleb_coarse = [1., 0.93, 0.92, 0.89, 0.85, 0.7]
doubleb_coarse = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", doubleb_coarse[::-1])
doublec_coarse = [0.87, 0.84, 0.83, 0.79, 0.69, 0.58]
doublec_coarse = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", doublec_coarse[::-1])
doublecvb_coarse = [0.93, 0.91, 0.86, 0.76, 0.6, 0.17, 0.12]
doublecvb_coarse = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", doublecvb_coarse[::-1])
n2ddt_coarse = hist.Bin("AK8Puppijet0_N2sdb1_ddt", "N2 DDT", [0.])


hists = {}
hists['sumw'] = hist.Hist("sumw", dataset, hist.Bin("sumw", "Weight value", [0.]))
hists['hjetpt'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
hists['hjetpt_SR'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
#hists['htagtensor'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, n2ddt_coarse, jetmass_coarse, doubleb, doublec, doublecvb, dtype='f')
hists['hsculpt'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
hists['hsculpt_SR'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')

hists['pfmet_nminus1_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("pfmet", r"PF $p_{T}^{miss}$", 40, 0, 200))
hists['opposite_ak8_n3sdb1_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0, 4))
hists['opposite_ak8_tau32_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 4))
hists['opposite_ak4_leadingDeepCSV_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak4_leadingDeepCSV", "Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1))

branches = [
    "AK8Puppijet0_pt",
    "AK8Puppijet0_eta",
    "AK8Puppijet0_msd",
    "AK8Puppijet0_isHadronicV",
    "AK8Puppijet0_deepdoubleb",
    "AK8Puppijet0_deepdoublec",
    "AK8Puppijet0_deepdoublecvb",
    "AK8Puppijet0_N2sdb1",
    "AK8Puppijet0_isTightVJet",
    "nAK4PuppijetsPt30dR08_0",
    "npu",
    "scale1fb",
    "kfactorEWK",
    "kfactorQCD",
    "neleLoose",
    "nmuLoose",
    "ntau",
    "pfmet",
    "AK8Puppijet0_phi",
    "AK8Puppijet1_phi",
    "AK8Puppijet1_e4_v2_sdb1",
    "AK8Puppijet1_e3_v1_sdb1",
    "AK8Puppijet1_tau32",
    "AK4Puppijet0_dPhi08",
    "AK4Puppijet1_dPhi08",
    "AK4Puppijet2_dPhi08",
    "AK4Puppijet3_dPhi08",
    "AK4Puppijet0_deepcsvb",
    "AK4Puppijet1_deepcsvb",
    "AK4Puppijet2_deepcsvb",
    "AK4Puppijet3_deepcsvb",
]

tstart = time.time()


for h in hists.values(): h.clear()
nevents = defaultdict(lambda: 0.)

def clean(val, default):
    val[np.isnan(val)|(val==-999.)] = default
    return val

def processfile(dataset, file):
    # Many 'invalid value encountered in ...' due to pt and msd sometimes being zero
    # This will just fill some NaN bins in the histogram, which is fine
    tree = uproot.open(file)["Events"]
    arrays = tree.arrays(branches, namedecode='ascii')

    # jet |eta|<2.5 sometimes gives no events
    # or other cuts in: https://github.com/DAZSLE/BaconAnalyzer/blob/102x/Analyzer/src/VJetLoader.cc#L270-L272
    arrays["AK8Puppijet0_pt"] = clean(arrays["AK8Puppijet0_pt"], 0.001)
    arrays["AK8Puppijet0_N2sdb1"] = clean(arrays["AK8Puppijet0_N2sdb1"], np.inf)

    # we'll take care of cross section later, just check if +/-1
    genW = np.sign(arrays["scale1fb"])
    weight = genW
    if dataset in pileup_corr:
        weight *= pileup_corr[dataset](arrays["npu"])
    if 'ZJetsToQQ_HT' in dataset or 'WJetsToQQ_HT' in dataset:
        weight *= arrays["kfactorEWK"] * arrays["kfactorQCD"]
    weight *= jetTriggerEff(arrays["AK8Puppijet0_msd"], arrays["AK8Puppijet0_pt"])
    weight *= (arrays["AK8Puppijet0_pt"] > 200)

    weight_SR = weight * ((arrays["neleLoose"]==0) & (arrays["nmuLoose"]==0) & (arrays["ntau"]==0))

    arrays["AK8Puppijet0_msd"] *= msd_weight(arrays["AK8Puppijet0_pt"], arrays["AK8Puppijet0_eta"])
    arrays["jetrho"] = 2*np.log(np.maximum(1e-4, arrays["AK8Puppijet0_msd"]/arrays["AK8Puppijet0_pt"]))
    arrays["AK8Puppijet0_N2sdb1_ddt"] = arrays["AK8Puppijet0_N2sdb1"] - n2ddt_rho_pt(arrays["jetrho"], arrays["AK8Puppijet0_pt"])
    weight_SR *= (arrays["AK8Puppijet0_N2sdb1_ddt"] < 0) & (arrays["AK8Puppijet0_isTightVJet"]!=0)
    weight_metCut = (arrays["pfmet"] < 140.)

    e4_v2_jet1 = clean(arrays['AK8Puppijet1_e4_v2_sdb1'], 1.)
    e3_v1_jet1 = clean(arrays['AK8Puppijet1_e3_v1_sdb1'], -1.)
    dphi = np.unwrap(arrays['AK8Puppijet1_phi'] - arrays['AK8Puppijet0_phi'])
    arrays['opposite_ak8_n3sdb1'] = np.where(np.abs(dphi) > np.pi/2., e4_v2_jet1/np.maximum(1e-4, e3_v1_jet1)**2, np.inf)
    arrays['opposite_ak8_tau32'] = np.where(np.abs(dphi) > np.pi/2., arrays['AK8Puppijet1_tau32'], np.inf)
    dphi04 = np.column_stack(arrays['AK4Puppijet%d_dPhi08' % i] for i in range(4))
    btag04 = np.column_stack(arrays['AK4Puppijet%d_deepcsvb' % i] for i in range(4))
    btag04[np.abs(dphi04)<np.pi/2] = -np.inf
    arrays['opposite_ak4_leadingDeepCSV'] = np.max(btag04, axis=1)

    hout = {}
    for k in hists.keys():
        h = hists[k].copy(content=False)
        if k == 'sumw':
            h.fill(dataset=dataset, sumw=genW)
        elif k == 'pfmet_nminus1_SR':
            h.fill(dataset=dataset, **arrays, weight=weight_SR)
        elif '_SR' in k:
            h.fill(dataset=dataset, **arrays, weight=weight_SR*weight_metCut)
        else:
            h.fill(dataset=dataset, **arrays, weight=weight)
        hout[k] = h
    return dataset, tree.numentries, hout


nworkers = 10
#fileslice = slice(None, 5)
fileslice = slice(None)
#with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = set()
    for dataset, info in datadef.items():
        futures.update(executor.submit(processfile, dataset, file) for file in info['files'][fileslice])
    try:
        total = len(futures)
        processed = 0
        while len(futures) > 0:
            finished = set(job for job in futures if job.done())
            for job in finished:
                dataset, nentries, hout = job.result()
                nevents[dataset] += nentries
                for k in hout.keys():
                    hists[k] += hout[k]
                processed += 1
                print("Processing: done with % 4d / % 4d files" % (processed, total))
            futures -= finished
        del finished
    except KeyboardInterrupt:
        print("Ok quitter")
        for job in futures: job.cancel()
    except:
        for job in futures: job.cancel()
        raise


sumw = hists.pop('sumw')
scale = {}
for ds in nevents.keys():
    ds_sumw = sumw.values(overflow='all')[(ds,)]
    print(ds, nevents[ds], ds_sumw)
    scale[ds] = lumi*dataset_xs[ds] / (ds_sumw[1]-ds_sumw[0])

for h in hists.values(): h.scale(scale, axis="dataset")

dt = time.time() - tstart
print("%.2f us*cpu/event" % (1e6*dt*nworkers/sum(nevents.values()), ))
nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in hists.values())
nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in hists.values())
print("Processed %.1fM events" % (sum(nevents.values())/1e6, ))
print("Filled %.1fM bins" % (nbins/1e6, ))
print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

# Pickle is not very fast or memory efficient, will be replaced by something better soon
with gzip.open("hists.pkl.gz", "wb") as fout:
    pickle.dump(hists, fout)


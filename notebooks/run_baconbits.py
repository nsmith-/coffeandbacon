from __future__ import print_function, division
from collections import defaultdict, OrderedDict, deque
import concurrent.futures
import gzip
import pickle
import json
import time
import numexpr

import uproot
import numpy as np
from fnal_column_analysis_tools import hist, lookup_tools

with open("data/datadef.json") as fin:
    datadef = json.load(fin)

extractor = lookup_tools.extractor()
extractor.add_weight_sets(["* * data/n2ddt_transform_2017MC.root"])
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

# [pb]
dataset_xs = {k: v['xs'] for k,v in datadef.items()}

lumi = 1000.  # [1/pb]

dataset = hist.Cat("dataset", "Primary dataset")

gencat = hist.Bin("AK8Puppijet0_isHadronicV", "Matched", 4, 0., 4)
# one can relabel intervals, although process mapping obviates this
titles = ["QCD", "V(light) matched", "V(c) matched", "V(b) matched"]
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
n2ddt = hist.Bin("AK8Puppijet0_N2sdb1_ddt", "N2 DDT", 20, -0.25, 0.25)
n2ddt_coarse = hist.Bin("AK8Puppijet0_N2sdb1_ddt", "N2 DDT", [-0.1, 0.])

hists = {}
hists['hjetpt'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
hists['htagtensor'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, n2ddt_coarse, jetmass_coarse, doubleb, doublec, doublecvb, dtype='f')
hists['hsculpt'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, n2ddt, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')

branches = [
    "AK8Puppijet0_pt",
    "AK8Puppijet0_eta",
    "AK8Puppijet0_msd",
    "AK8Puppijet0_isHadronicV",
    "AK8Puppijet0_deepdoubleb",
    "AK8Puppijet0_deepdoublec",
    "AK8Puppijet0_deepdoublecvb",
    "AK8Puppijet0_N2sdb1",
]

tstart = time.time()


for h in hists.values(): h.clear()
nevents = defaultdict(lambda: 0.)


def processfile(dataset, file):
    tree = uproot.open(file)["Events"]
    arrays = tree.arrays(branches, namedecode='ascii')
    arrays["AK8Puppijet0_msd"] *= msd_weight(arrays["AK8Puppijet0_pt"], arrays["AK8Puppijet0_eta"])
    arrays["jetrho"] = 2*np.log(np.maximum(arrays["AK8Puppijet0_msd"], 0.01)/arrays["AK8Puppijet0_pt"])
    arrays["AK8Puppijet0_N2sdb1_ddt"] = arrays["AK8Puppijet0_N2sdb1"] - n2ddt_rho_pt(arrays["jetrho"], arrays["AK8Puppijet0_pt"])
    hout = {}
    for k in hists.keys():
        h = hists[k].copy(content=False)
        h.fill(dataset=dataset, **arrays)
        hout[k] = h
    return dataset, tree.numentries, hout


nworkers = 10
#fileslice = slice(None, 5)
fileslice = slice(None)
#with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = deque()
    for dataset, info in datadef.items():
        futures.extend(executor.submit(processfile, dataset, file) for file in info['files'][fileslice])
    try:
        nfiles = len(futures)
        for i in range(nfiles):
            fut = futures.popleft()
            dataset, nentries, hout = fut.result()
            nevents[dataset] += nentries
            for k in hout.keys():
                hists[k] += hout[k]
            print("Processing: done with % 4d / % 4d files" % (i+1, nfiles))
            del fut
    except KeyboardInterrupt:
        print("Ok quitter")
        for fut in futures: fut.cancel()
    except:
        for fut in futures: fut.cancel()
        raise

scale = dict((ds, lumi * dataset_xs[ds] / nevents[ds]) for ds in nevents.keys())
for h in hists.values(): h.scale(scale, axis="dataset")

dt = time.time() - tstart
print("%.2f us*cpu/event" % (1e6*dt*nworkers/sum(nevents.values()), ))
nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in hists.values())
nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in hists.values())
print("Processed %.1fM events" % (sum(nevents.values())/1e6, ))
print("Filled %.1fM bins" % (nbins/1e6, ))
print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

with gzip.open("hists.pkl.gz", "wb") as fout:
    pickle.dump(hists, fout)


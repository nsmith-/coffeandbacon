#!/usr/bin/env python
from __future__ import print_function, division
import collections
import warnings
import concurrent.futures
import gzip
import pickle
import cloudpickle
import json
import time
import numexpr

import uproot
import awkward
import numpy as np
from fnal_column_analysis_tools import hist, lookup_tools

with open("metadata/samplefiles.json") as fin:
    samplefiles = json.load(fin)


# instrument xrootd source
def _read(self, chunkindex):
    self.bytesread = getattr(self, 'bytesread', 0) + self._chunkbytes
    return self._read_real(chunkindex)

uproot.source.xrootd.XRootDSource._read_real = uproot.source.xrootd.XRootDSource._read
uproot.source.xrootd.XRootDSource._read = _read

with gzip.open("corrections.cpkl.gz", "rb") as fin:
    corrections = cloudpickle.load(fin)


# axis definitions
dataset = hist.Cat("dataset", "Primary dataset")
gencat = hist.Bin("AK8Puppijet0_isHadronicV", "Matched", [0,1,2,3,9,10,11])
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
doublecvb_coarse = [0.93, 0.91, 0.6, 0.2, 0.17]
doublecvb_coarse = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", doublecvb_coarse[::-1])
n2ddt_coarse = hist.Bin("AK8Puppijet0_N2sdb1_ddt", "N2 DDT", [0.])


hists = {}
hists['sumw'] = hist.Hist("sumw", dataset, hist.Bin("sumw", "Weight value", [0.]))
hists['hjetpt'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
hists['hjetpt_SR'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
#hists['hsculpt'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
hists['hsculpt_SR'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
hists['htagtensor_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, doubleb, doublec, doublecvb, dtype='f')

hists['pfmet_nminus1_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("pfmet", r"PF $p_{T}^{miss}$", 40, 0, 200))
hists['opposite_ak8_n3sdb1_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0.5, 3))
hists['opposite_ak8_tau32_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 1))
hists['opposite_ak8_msd_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_msd", r"Jet $\m_{sd}$", 40, 50, 200))
hists['opposite_ak4_leadingDeepCSV_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak4_leadingDeepCSV", "Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1))
hists['njets_ak4_SR'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("nAK4PuppijetsPt30", "Number AK4 Jets", 8, 0, 8))


def clean(val, default):
    val[np.isnan(val)|(val==-999.)] = default
    return val


def process(df):
    isData = df['dataset'] == "data_obs"
    # jet |eta|<2.5 sometimes gives no events
    # or other cuts in: https://github.com/DAZSLE/BaconAnalyzer/blob/102x/Analyzer/src/VJetLoader.cc#L270-L272
    df["AK8Puppijet0_pt"] = clean(df["AK8Puppijet0_pt"], 0.001)
    df["AK8Puppijet0_N2sdb1"] = clean(df["AK8Puppijet0_N2sdb1"], np.inf)

    # we'll take care of cross section later, just check if +/-1
    genW = np.sign(df["scale1fb"])
    weight = genW
    if dataset in corrections['2017_pileupweight_dataset']:
        weight *= corrections['2017_pileupweight_dataset'][dataset](df["npu"])
    if 'ZJetsToQQ_HT' in dataset or 'WJetsToQQ_HT' in dataset:
        weight *= df["kfactorEWK"] * df["kfactorQCD"]
    weight *= corrections['2017_trigweight_msd_pt'](df["AK8Puppijet0_msd"], df["AK8Puppijet0_pt"])
    weight *= (df["AK8Puppijet0_pt"] > 200)

    weight_SR = weight * ((df["neleLoose"]==0) & (df["nmuLoose"]==0) & (df["ntau"]==0))

    df["AK8Puppijet0_msd"] *= corrections['msdweight'](df["AK8Puppijet0_pt"], df["AK8Puppijet0_eta"])
    df["jetrho"] = 2*np.log(np.maximum(1e-4, df["AK8Puppijet0_msd"]/df["AK8Puppijet0_pt"]))
    df["AK8Puppijet0_N2sdb1_ddt"] = df["AK8Puppijet0_N2sdb1"] - corrections['2017_n2ddt_rho_pt'](df["jetrho"], df["AK8Puppijet0_pt"])
    weight_SR *= (df["AK8Puppijet0_N2sdb1_ddt"] < 0) & (df["AK8Puppijet0_isTightVJet"]!=0)
    weight_metCut = (df["pfmet"] < 140.)

    e4_v2_jet1 = clean(df['AK8Puppijet1_e4_v2_sdb1'], 1.)
    e3_v1_jet1 = clean(df['AK8Puppijet1_e3_v1_sdb1'], -1.)
    dphi = np.unwrap(df['AK8Puppijet1_phi'] - df['AK8Puppijet0_phi'])
    df['opposite_ak8_n3sdb1'] = np.where(np.abs(dphi) > np.pi/2., e4_v2_jet1/np.maximum(1e-4, e3_v1_jet1)**2, np.inf)
    df['opposite_ak8_tau32'] = np.where(np.abs(dphi) > np.pi/2., df['AK8Puppijet1_tau32'], np.inf)
    df['opposite_ak8_msd'] = np.where(np.abs(dphi) > np.pi/2., df['AK8Puppijet1_msd'], np.inf)
    dphi04 = np.column_stack([df['AK4Puppijet%d_dPhi08' % i] for i in range(4)])
    btag04 = np.column_stack([df['AK4Puppijet%d_deepcsvb' % i] for i in range(4)])
    btag04[np.abs(dphi04)<np.pi/2] = -np.inf
    df['opposite_ak4_leadingDeepCSV'] = np.max(btag04, axis=1)

    hout = {}
    for k in hists.keys():
        h = hists[k].copy(content=False)
        fields = {k: df[k] for k in h.fields if k in df}
        if k == 'sumw':
            if 'skim_sumw' in df:
                h.fill(dataset=dataset, sumw=1, weight=df['skim_sumw'])
            else:
                h.fill(dataset=dataset, sumw=genW)
        elif k == 'pfmet_nminus1_SR':
            h.fill(**fields, weight=weight_SR)
        elif '_SR' in k:
            h.fill(**fields, weight=weight_SR*weight_metCut)
        else:
            h.fill(**fields, weight=weight)
        hout[k] = h

    return hout


class LazyArrays(collections.abc.MutableMapping):
    def __init__(self, tree):
        self._tree = tree
        self._dict = self._tree.lazyarrays(namedecode='ascii')
        self._materialized = set()

    def __delitem__(self, key):
        del self._dict[key]

    def __getitem__(self, key):
        if key in self._dict:
            value = self._dict[key]
            if isinstance(value, uproot.tree.LazyArray):
                self._materialized.add(key)
                return value[:]
            return value
        else:
            raise KeyError(key)

    def __iter__(self):
        print("uh oh")
        for item in self._dict:
            self._materialized.add(item[0])
            yield item

    def __len__(self):
        return len(self._dict)

    def __setitem__(self, key, value):
        self._dict[key] = value

    @property
    def materialized(self):
        return self._materialized


def processfile(dataset, file):
    fin = uproot.open(file)
    skim_sumw = None
    if 'otree' in fin:
        tree = fin['otree']
        if 'SumWeights' in fin:
            skim_sumw = fin['SumWeights'].values[0]
    else:
        tree = fin['Events']

    df = LazyArrays(tree)
    df['dataset'] = dataset
    df['skim_sumw'] = skim_sumw
    tic = time.time()

    hout = process(df)

    toc = time.time()
    bytesread = fin.source.bytesread if isinstance(fin.source, uproot.source.xrootd.XRootDSource) else 0
    output = {
        'dataset': dataset,
        'nentries': tree.numentries,
        'histograms': hout,
        'bytesread': bytesread,
        'elapsedtime': toc-tic,
        'columns': df.materialized,
    }
    return output


test = False

tstart = time.time()
for h in hists.values(): h.clear()
nevents = collections.defaultdict(lambda: 0.)
nbytes = collections.defaultdict(lambda: 0.)
sumworktime = 0.
columns_accessed = set()

def collect(output):
    global sumworktime, columns_accessed
    dataset = output['dataset']
    nevents[dataset] += output['nentries']
    nbytes[dataset] += output['bytesread']
    sumworktime += output['elapsedtime']
    for k, v in output['histograms'].items():
        hists[k] += v
    columns_accessed.update(output['columns'])


if test:
    nworkers = 1
    testfiles = [
        ("TTToHadronic_TuneCP5_13TeV_powheg_pythia8", "TTToHadronic_TuneCP5_13TeV_powheg_pythia8_0.root"),
        ("TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_0.root"),
    ]
    for i, (dataset, file) in enumerate(testfiles):
        collect(processfile(dataset, file))
        print("Done processing test file %d" % i)
else:
    nworkers = 10
    #fileslice = slice(None, 5)
    fileslice = slice(None)
    #with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
    with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
        futures = set()
        samples = samplefiles['Hbb_create_2017']
        for datasets in samples.values():
            if isinstance(datasets, list):
                # raw data, no norm
                dataset = "data_obs"
                futures.update(executor.submit(processfile, dataset, file) for file in datasets)
            elif isinstance(datasets, dict):
                for dataset, files in datasets.items():
                    futures.update(executor.submit(processfile, dataset, file) for file in files)
        try:
            total = len(futures)
            processed = 0
            while len(futures) > 0:
                finished = set(job for job in futures if job.done())
                for job in finished:
                    collect(job.result())
                    processed += 1
                    print("Processing: done with % 4d / % 4d files" % (processed, total))
                futures -= finished
                del finished
                time.sleep(1)
        except KeyboardInterrupt:
            print("Ok quitter")
            for job in futures:
                job.cancel()
            print("Killed pending jobs")
            print("Running jobs:", sum(1 for j in futures if j.running()))
        except:
            for job in futures: job.cancel()
            raise


def read_xsections(filename):
    out = {}
    with open(filename) as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            dataset, xsexpr, *_ = line.split()
            try:
                xs = float(numexpr.evaluate(xsexpr))
            except:
                print("numexpr evaluation failed for line: %s" % line)
                raise
            if xs <= 0:
                warnings.warn("Cross section is <= 0 in line: %s" % line, RuntimeWarning)
            out[dataset] = xs
    return out


# curl -O https://raw.githubusercontent.com/kakwok/ZPrimePlusJet/DDB/analysis/ggH/xSections.dat
xsections = read_xsections("metadata/xSections.dat")
lumi = 41100  # [1/pb]

sumw = hists.pop('sumw')
scale = {}
for ds in nevents.keys():
    ds_sumw = sumw.values(overflow='all')[(ds,)]
    print(ds, nevents[ds], ds_sumw)
    if ds != "data_obs":
        scale[ds] = lumi*xsections[ds] / (ds_sumw[1]-ds_sumw[0])

for h in hists.values(): h.scale(scale, axis="dataset")

dt = time.time() - tstart
print("Columns accessed:")
for col in sorted(list(columns_accessed)):
    print("\t", col)
print("%.2f us*cpu/event" % (1e6*dt*nworkers/sum(nevents.values()), ))
print("%.2f us*cpu/event work time" % (1e6*sumworktime/sum(nevents.values()), ))
nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in hists.values())
nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in hists.values())
print("Processed %.1fM events" % (sum(nevents.values())/1e6, ))
print("Read %.1fM bytes" % (sum(nbytes.values())/1e6, ))
print("Filled %.1fM bins" % (nbins/1e6, ))
print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

# Pickle is not very fast or memory efficient, will be replaced by something better soon
with gzip.open("hists.pkl.gz", "wb") as fout:
    pickle.dump(hists, fout)


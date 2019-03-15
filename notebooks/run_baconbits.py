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
jetrho = hist.Bin("ak8jet_rho", r"Jet $\rho$", 13, -6, -2.1)
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
hists['hjetpt_signalregion'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
hists['hsculpt_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
hists['htagtensor_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, doubleb, doublec, doublecvb, dtype='f')
hists['opposite_ak8_n3sdb1_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0.5, 3))
hists['opposite_ak8_tau32_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 1))
hists['opposite_ak8_msd_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_msd", r"Jet $\m_{sd}$", 40, 50, 200))
hists['opposite_ak4_leadingDeepCSV_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak4_leadingDeepCSV", "Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1))
hists['njets_ak4_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("nAK4PuppijetsPt30", "Number AK4 Jets", 8, 0, 8))

hists['nminus1_pfmet140_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("pfmet", r"PF $p_{T}^{miss}$", 40, 0, 200))


class PackedSelection(object):
    def __init__(self, dtype='uint64'):
        self._dtype = np.dtype(dtype)
        self._names = []
        self._mask = None

    def add(self, name, selection):
        if isinstance(selection, np.ndarray) and selection.dtype == np.dtype('bool'):
            if len(self._names) == 0:
                self._mask = np.zeros(shape=selection.shape, dtype=self._dtype)
            elif len(self._names) == 64:
                raise RuntimeError("Exhausted all slots for %r, consider a larger dtype or fewer selections" % self._dtype)
            elif self._mask.shape != selection.shape:
                raise ValueError("New selection '%s' has different shape than existing ones (%r vs. %r)" % (name, selection.shape, self._mask.shape))
            self._mask |= selection.astype(self._dtype) << len(self._names)
            self._names.append(name)
        else:
            raise ValueError("PackedSelection only understands numpy boolean arrays, got %r" % selection)

    def require(self, **names):
        mask = 0
        require = 0
        for name, val in names.items():
            if not isinstance(val, bool):
                raise ValueError("Please use only booleans in PackedSelection.require(), received %r for %s" % (val, name))
            idx = self._names.index(name)
            mask |= 1<<idx
            require |= int(val)<<idx
        return (self._mask & mask) == require

    def all(self, *names):
        if len(names) == 1 and isinstance(names[0], collections.abc.Iterable):
            names = names[0]
        return self.require(**{name: True for name in names})


class Weights(object):
    def __init__(self, size):
        self._weight = np.ones(size)
        self._modifiers = {}
        self._weightStats = {}

    def add(self, name, weight, weightUp=None, weightDown=None, shift=False):
        """
            name: name of correction weight
            weight: nominal weight
            weightUp: weight with correction uncertainty shifted up
            weightDown: weight with correction uncertainty shifted down (leave None if symmetric)
            shift: if True, interpret weightUp and weightDown as a difference relative to the nominal value
        """
        self._weight *= weight
        if weightUp is not None:
            if shift:
                weightUp += weight
            weightUp[weight != 0.] /= weight[weight != 0.]
            self._modifiers[name+'Up'] = weightUp
        if weightDown is not None:
            if shift:
                weightDown = weight - weightDown
            weightDown[weight != 0.] /= weight[weight != 0.]
            self._modifiers[name+'Down'] = weightDown
        self._weightStats[name] = {
            'sumw': weight.sum(),
            'sumw2': (weight**2).sum(),
            'min': weight.min(),
            'max': weight.max(),
            'n': weight.size,
        }

    def weight(self, modifier=None):
        if modifier is None:
            return self._weight
        elif 'Down' in modifier and modifier not in self._modifiers:
            return self._weight / self._modifiers[modifier.replace('Down', 'Up')]
        return self._weight * self._modifiers[modifier]


def clean(df, val, default, positive=False):
    if positive:
        df[val][np.isnan(df[val])|(df[val]<=0.)] = default
    else:
        df[val][np.isnan(df[val])|(df[val]==-999.)] = default


def build_leading_ak8_variables(df):
    # jet |eta|<2.5 sometimes gives no events
    # or other cuts in: https://github.com/DAZSLE/BaconAnalyzer/blob/102x/Analyzer/src/VJetLoader.cc#L270-L272
    # set safe dummy values to avoid domain errors (FPU exceptions slow things down!)
    clean(df, 'AK8Puppijet0_pt', 100., positive=True) # msdweight goes negative for pt < 13
    clean(df, 'AK8Puppijet0_msd', 1e-7, positive=True)
    clean(df, 'AK8Puppijet0_N2sdb1', np.inf)
    clean(df, 'AK8Puppijet0_pt_JESUp', 1e-3)
    clean(df, 'AK8Puppijet0_pt_JESDown', 1e-3)
    clean(df, 'AK8Puppijet0_pt_JERUp', 1e-3)
    clean(df, 'AK8Puppijet0_pt_JERDown', 1e-3)
    clean(df, 'AK8Puppijet0_deepdoubleb', -1.)
    df['AK8Puppijet0_msd'] *= corrections['msdweight'](df['AK8Puppijet0_pt'], df['AK8Puppijet0_eta'])
    df['ak8jet_rho'] = 2*np.log(df['AK8Puppijet0_msd']/df['AK8Puppijet0_pt'])
    df['ak8jet_n2ddt'] = df['AK8Puppijet0_N2sdb1'] - corrections['2017_n2ddt_rho_pt'](df['ak8jet_rho'], df['AK8Puppijet0_pt'])


def subleading_n3(df):
    e4_v2_jet1 = clean(df, 'AK8Puppijet1_e4_v2_sdb1', 1.)
    e3_v1_jet1 = clean(df, 'AK8Puppijet1_e3_v1_sdb1', 1e-4, positive=True)
    return df['AK8Puppijet1_e4_v2_sdb1']/df['AK8Puppijet1_e3_v1_sdb1']**2


def build_subleading_ak8_variables(df):
    dphi = np.abs(np.unwrap(df['AK8Puppijet1_phi'] - df['AK8Puppijet0_phi']))
    df['opposite_ak8_n3sdb1'] = np.where(dphi > np.pi/2., subleading_n3(df), np.inf)
    df['opposite_ak8_tau32'] = np.where(dphi > np.pi/2., df['AK8Puppijet1_tau32'], np.inf)
    df['opposite_ak8_msd'] = np.where(dphi > np.pi/2., df['AK8Puppijet1_msd'], np.inf)


def build_ak4_variables(df):
    # dR08, dPhi08 with respect to leading ak8 jet: https://github.com/DAZSLE/BaconAnalyzer/blob/102x/Analyzer/src/JetLoader.cc#L478-L479
    n_ak4 = 4
    stack = lambda var: np.column_stack([df['AK4Puppijet%d_%s' % (i, var)] for i in range(n_ak4)])
    dR = stack('dR08')
    dphi = stack('dPhi08')
    btag = stack('deepcsvb')
    pt = stack('pt')
    # seems |eta|<2.5 already in tuple
    require = (np.abs(dphi) > np.pi/2) & (pt > 30.)
    btag_ttrej = np.where(require, btag, -np.inf)
    df['opposite_ak4_leadingDeepCSV'] = np.max(btag_ttrej, axis=1)
    require = (dR > 0.8) & (pt > 50.)
    btag_muCR = np.where(require, btag, -np.inf)
    df['ak4_leadingDeepCSV_dR08'] = np.max(btag_muCR, axis=1)


def build_met_systematics(df):
    metx = df['pfmet']*np.sin(df['pfmetphi'])
    mety = df['pfmet']*np.cos(df['pfmetphi'])
    df['pfmet_JESUp'] = np.hypot(metx + df['MetXCorrjesUp'], mety + df['MetYCorrjesUp'])
    df['pfmet_JESDown'] = np.hypot(metx + df['MetXCorrjesDown'], mety + df['MetYCorrjesDown'])
    df['pfmet_JERUp'] = np.hypot(metx + df['MetXCorrjerUp'], mety + df['MetYCorrjerUp'])
    df['pfmet_JERDown'] = np.hypot(metx + df['MetXCorrjerDown'], mety + df['MetYCorrjerDown'])


def process(df):
    isData = df['dataset'] == 'data_obs'

    weights = Weights(df.size)

    # we'll take care of cross section later, just check if +/-1
    if not isData:
        weights.add('genweight', np.sign(df['scale1fb']))

    if dataset in corrections['2017_pileupweight_dataset']:
        weights.add('pileupweight',
                    corrections['2017_pileupweight_dataset'][dataset](df['npu']),
                    corrections['2017_pileupweight_dataset_puUp'][dataset](df['npu']),
                    corrections['2017_pileupweight_dataset_puUp'][dataset](df['npu']),
                    )

    if 'ZJetsToQQ_HT' in dataset or 'WJetsToQQ_HT' in dataset:
        weights.add('kfactor', df['kfactorEWK'] * df['kfactorQCD'])
        # TODO unc.

    # trigger weight uses uncorrected jet mass
    if not isData:
        weights.add('trigweight',
                    corrections['2017_trigweight_msd_pt'](df['AK8Puppijet0_msd'], df['AK8Puppijet0_pt']),
                    corrections['2017_trigweight_msd_pt_trigweightUp'](df['AK8Puppijet0_msd'], df['AK8Puppijet0_pt']),
                    corrections['2017_trigweight_msd_pt_trigweightDown'](df['AK8Puppijet0_msd'], df['AK8Puppijet0_pt']),
                    )

    # muon CR weights
    if not isData:
        mu_abseta = np.abs(df['vmuoLoose0_eta'])
        weights.add('mutrigweight',
                    corrections['2017_mutrigweight_pt_abseta'](df['vmuoLoose0_pt'], mu_abseta),
                    corrections['2017_mutrigweight_pt_abseta_mutrigweightShift'](df['vmuoLoose0_pt'], mu_abseta),
                    shift=True
                    )
        weights.add('muidweight',
                    corrections['2017_muidweight_abseta_pt'](mu_abseta, df['vmuoLoose0_pt']),
                    corrections['2017_muidweight_abseta_pt_muidweightShift'](mu_abseta, df['vmuoLoose0_pt']),
                    shift=True
                    )
        weights.add('muisoweight',
                    corrections['2017_muisoweight_abseta_pt'](mu_abseta, df['vmuoLoose0_pt']),
                    corrections['2017_muisoweight_abseta_pt_muisoweightShift'](mu_abseta, df['vmuoLoose0_pt']),
                    shift=True
                    )


    build_leading_ak8_variables(df)
    build_subleading_ak8_variables(df)
    build_ak4_variables(df)
    build_met_systematics(df)

    selection = PackedSelection()
    if isData:
        selection.add('trigger', df['triggerBits'] & corrections['2017_triggerMask'])
    else:
        selection.add('trigger', np.ones(df.size, dtype='bool'))

    selection.add('noLeptons', (df['neleLoose']==0) & (df['nmuLoose']==0) & (df['ntau']==0))
    selection.add('oneMuon', (df['neleLoose']==0) & (df['nmuLoose']==1) & (df['ntau']==0))
    selection.add('muonAcceptance', (df['vmuoLoose0_pt'] > 55.) & (np.abs(df['vmuoLoose0_eta']) < 2.1))
    selection.add('ak4btagMediumDR08', df['ak4_leadingDeepCSV_dR08'] > 0.4941)  # at least one passes medium cut
    selection.add('muonDphiAK8', np.abs(np.unwrap(df['vmuoLoose0_phi'] - df['AK8Puppijet0_phi'])) > 2*np.pi/3)
    selection.add('antiak4btagMediumOppHem', df['opposite_ak4_leadingDeepCSV'] < 0.4941)  # none pass
    selection.add('tightVjet', df['AK8Puppijet0_isTightVJet'] != 0)
    selection.add('n2ddtPass', df['ak8jet_n2ddt'] < 0)
    selection.add('doublebtagPass', df['AK8Puppijet0_deepdoubleb'] > 0.9)
    selection.add('jetMass', df['AK8Puppijet0_msd'] > 40.)

    selection.add('jetKinematics', df['AK8Puppijet0_pt'] > 450.)
    selection.add('jetKinematicsMuonCR', df['AK8Puppijet0_pt'] > 400.)
    selection.add('pfmet140', df['pfmet'] < 140.)

    for syst in ['JESUp', 'JESDown', 'JERUp', 'JERDown']:
        selection.add('jetKinematics'+syst, df['AK8Puppijet0_pt_'+syst] > 450)
        selection.add('jetKinematicsMuonCR'+syst, df['AK8Puppijet0_pt_'+syst] > 400.)
        selection.add('pfmet140'+syst, df['pfmet_'+syst] < 140.)

    regions = {}
    regions['signalregion'] = {'trigger', 'n2ddtPass', 'noLeptons', 'jetKinematics', 'tightVjet', 'doublebtagPass', 'antiak4btagMediumOppHem'}
    # TODO: mutrigger
    regions['muoncontrol'] = {'n2ddtPass', 'oneMuon', 'jetKinematicsMuonCR', 'tightVjet', 'doublebtagPass', 'ak4btagMediumDR08', 'muonDphiAK8'}

    print(weights._weightStats)
    hout = {}
    for histname in hists.keys():
        h = hists[histname].copy(content=False)
        fields = {k: df[k] for k in h.fields if k in df}
        if histname == 'sumw':
            if 'skim_sumw' in df:
                h.fill(dataset=dataset, sumw=1, weight=df['skim_sumw'])
            else:
                h.fill(dataset=dataset, sumw=df['scale1fb'])
        elif 'nminus1' in histname:
            _, sel, region = histname.split('_')
            cut = regions[region] - {sel}
            weight = weights.weight() * selection.all(cut)
            h.fill(**fields, weight=weight)
        # TODO: nested hists? search region mames in histname?
        elif 'signalregion' in histname:
            region = 'signalregion'
            cut = regions[region]
            weight = weights.weight() * selection.all(cut)
            h.fill(**fields, weight=weight)
        else:
            weight = weights.weight()
            h.fill(**fields, weight=weight)
        hout[histname] = h

    return hout


class DataFrame(collections.abc.MutableMapping):
    def __init__(self, tree):
        self._tree = tree
        self._dict = {}
        self._materialized = set()

    def __delitem__(self, key):
        del self._dict[key]

    def __getitem__(self, key):
        if key in self._dict:
            return self._dict[key]
        elif key in self._tree:
            self._materialized.add(key)
            self._dict[key] = self._tree[key].array()
            return self._dict[key]
        else:
            raise KeyError(key)

    def __iter__(self):
        print('uh oh')
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

    @property
    def size(self):
        return self._tree.numentries


def processfile(dataset, file):
    fin = uproot.open(file)
    skim_sumw = None
    if 'otree' in fin:
        tree = fin['otree']
        if 'SumWeights' in fin:
            skim_sumw = fin['SumWeights'].values[0]
    else:
        tree = fin['Events']

    df = DataFrame(tree)
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


test = True

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


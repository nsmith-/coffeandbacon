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
import pprint
import copy

import uproot
import awkward
import numpy as np
from fnal_column_analysis_tools import hist, lookup_tools, processor

test = False

if test:
    from pyinstrument import Profiler
    profiler = Profiler()
    profiler.start()

tstart = time.time()

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
gencat = hist.Bin("AK8Puppijet0_isHadronicV", "V matching index", [0,1,2,3,9,10,11])
jetpt = hist.Bin("AK8Puppijet0_pt", r"Jet $p_T$", [450, 500, 550, 600, 675, 800, 1200])
jetmass = hist.Bin("AK8Puppijet0_msd", r"Jet $m_{sd}$", 23, 40, 201)
jetpt_coarse = hist.Bin("AK8Puppijet0_pt", r"Jet $p_T$", [450, 1200])
jetmass_coarse = hist.Bin("AK8Puppijet0_msd", r"Jet $m_{sd}$", [40, 103, 152, 201])
jetrho = hist.Bin("ak8jet_rho", r"Jet $\rho$", 13, -6, -2.1)
doubleb = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", 20, 0., 1)
doublec = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", 20, 0., 1.)
doublecvb = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", 20, 0., 1.)
doubleb_coarse = [1., 0.9, 0.89, 0.85, 0.7]
doubleb_coarse = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", doubleb_coarse[::-1])
doublec_coarse = [0.87, 0.84, 0.83, 0.79, 0.69]
doublec_coarse = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", doublec_coarse[::-1])
doublecvb_coarse = [0.93, 0.91, 0.6, 0.2, 0.17]
doublecvb_coarse = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", doublecvb_coarse[::-1])


accumulator_def = {}
accumulator_def['sumw'] = hist.Hist("sumw", dataset, hist.Bin("sumw", "Weight value", [0.]))
accumulator_def['hjetpt_signalregion'] = hist.Hist("Events", dataset, gencat, hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300), dtype='f')
accumulator_def['hsculpt_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt, jetmass, doubleb_coarse, doublec_coarse, doublecvb_coarse, dtype='f')
accumulator_def['htagtensor_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, doubleb, doublec, doublecvb, dtype='f')
accumulator_def['opposite_ak8_n3sdb1_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0.5, 3))
accumulator_def['opposite_ak8_tau32_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 1))
accumulator_def['opposite_ak8_msd_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak8_msd", r"Jet $\m_{sd}$", 40, 50, 200))
accumulator_def['opposite_ak4_leadingDeepCSV_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("opposite_ak4_leadingDeepCSV", "Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1))
accumulator_def['njets_ak4_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, hist.Bin("nAK4PuppijetsPt30", "Number AK4 Jets", 8, 0, 8))

accumulator_def['nminus1_pfmet_signalregion'] = hist.Hist("Events", dataset, gencat, jetpt_coarse, jetmass_coarse, doubleb_coarse, hist.Bin("pfmet", r"PF $p_{T}^{miss}$", 40, 0, 200))
accumulator_def['nminus1_n2ddtPass_signalregion'] = hist.Hist("Events", dataset, gencat, jetmass_coarse, doubleb_coarse, hist.Bin("ak8jet_n2ddt", r"Jet $N_{2,DDT}^{\beta=1}$", 40, -1, 1))
accumulator_def['templates_signalregion'] = hist.Hist("Events", dataset, gencat, hist.Cat("systematic", "Systematic"), jetpt, jetmass, doubleb_coarse)
accumulator_def['templates_muoncontrol'] = hist.Hist("Events", dataset, gencat, hist.Cat("systematic", "Systematic"), jetpt, jetmass, doubleb_coarse)

accumulator_def['nentries'] = 0
accumulator_def['bytesread'] = 0
accumulator_def['sumworktime'] = 0.
accumulator_def['columns_accessed'] = set()


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
    dataset = df['dataset']
    isRealData = dataset in ["JetHT", "SingleMuon"]

    weights = processor.Weights(df.size)

    # SumWeights is sum(scale1fb), so we need to use full value here
    if not isRealData:
        weights.add('genweight', df['scale1fb'])

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
    if not isRealData:
        weights.add('trigweight',
                    corrections['2017_trigweight_msd_pt'](df['AK8Puppijet0_msd'], df['AK8Puppijet0_pt']),
                    corrections['2017_trigweight_msd_pt_trigweightUp'](df['AK8Puppijet0_msd'], df['AK8Puppijet0_pt']),
                    corrections['2017_trigweight_msd_pt_trigweightDown'](df['AK8Puppijet0_msd'], df['AK8Puppijet0_pt']),
                    )

    # muon CR weights
    if not isRealData:
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

    selection = processor.PackedSelection()
    if isRealData:
        # Only take jet triggers from JetHT, single muon triggers from SingleMuon dataset
        # necessary but not sufficient condition to prevent double-counting
        # (this plus mutually exclusive offline selections are sufficient)
        selection.add('trigger', (df['triggerBits'] & corrections['2017_triggerMask']).astype('bool') & (dataset=="JetHT"))
        selection.add('mutrigger', ((df['triggerBits']&1) & df['passJson']).astype('bool') & (dataset=="SingleMuon"))
    else:
        selection.add('trigger', np.ones(df.size, dtype='bool'))
        selection.add('mutrigger', np.ones(df.size, dtype='bool'))

    selection.add('noLeptons', (df['neleLoose']==0) & (df['nmuLoose']==0) & (df['ntau']==0))
    selection.add('oneMuon', (df['neleLoose']==0) & (df['nmuLoose']==1) & (df['ntau']==0))
    selection.add('muonAcceptance', (df['vmuoLoose0_pt'] > 55.) & (np.abs(df['vmuoLoose0_eta']) < 2.1))
    selection.add('ak4btagMediumDR08', df['ak4_leadingDeepCSV_dR08'] > 0.4941)  # at least one passes medium cut
    selection.add('muonDphiAK8', np.abs(np.unwrap(df['vmuoLoose0_phi'] - df['AK8Puppijet0_phi'])) > 2*np.pi/3)
    selection.add('antiak4btagMediumOppHem', df['opposite_ak4_leadingDeepCSV'] < 0.4941)  # none pass
    selection.add('tightVjet', df['AK8Puppijet0_isTightVJet'] != 0)
    selection.add('n2ddtPass', df['ak8jet_n2ddt'] < 0)
    selection.add('jetMass', df['AK8Puppijet0_msd'] > 40.)

    selection.add('jetKinematics', df['AK8Puppijet0_pt'] > 450.)
    selection.add('jetKinematicsMuonCR', df['AK8Puppijet0_pt'] > 400.)
    selection.add('pfmet', df['pfmet'] < 140.)

    shiftSystematics = ['JESUp', 'JESDown', 'JERUp', 'JERDown']
    shiftedQuantities = {'AK8Puppijet0_pt', 'pfmet'}
    shiftedSelections = {'jetKinematics', 'jetKinematicsMuonCR', 'pfmet'}
    for syst in shiftSystematics:
        selection.add('jetKinematics'+syst, df['AK8Puppijet0_pt_'+syst] > 450)
        selection.add('jetKinematicsMuonCR'+syst, df['AK8Puppijet0_pt_'+syst] > 400.)
        selection.add('pfmet'+syst, df['pfmet_'+syst] < 140.)

    regions = {}
    regions['signalregion'] = {'trigger', 'n2ddtPass', 'noLeptons', 'jetKinematics', 'tightVjet', 'antiak4btagMediumOppHem'}
    regions['muoncontrol'] = {'mutrigger', 'n2ddtPass', 'oneMuon', 'jetKinematicsMuonCR', 'tightVjet', 'ak4btagMediumDR08', 'muonDphiAK8'}

    if test:
        print("Weight statistics:")
        pprint.pprint(weights._weightStats, indent=4)

    hout = {}
    for histname in accumulator_def.keys():
        h = accumulator_def[histname]
        if not isinstance(h, hist.Hist):
            continue
        h = h.copy(content=False)
        fields = {k: df[k] for k in h.fields if k in df}
        region = [r for r in regions.keys() if r in histname]

        if histname == 'sumw':
            if isRealData:
                pass
            elif 'skim_sumw' in df:
                # hacky way to only accumulate file-level information once
                if df['skim_sumw'] is not None:
                    h.fill(dataset=dataset, sumw=1, weight=df['skim_sumw'])
            else:
                h.fill(dataset=dataset, sumw=np.sign(df['scale1fb']))
        elif 'nminus1' in histname:
            _, sel, region = histname.split('_')
            cut = regions[region] - {sel}
            weight = weights.weight() * selection.all(cut)
            h.fill(**fields, weight=weight)
        elif len(region) == 1:
            region = region[0]
            weight = weights.weight()
            cut = selection.all(regions[region])
            h.fill(systematic="", **fields, weight=weight*cut)
            if 'systematic' in h.fields:
                for syst in weights.variations:
                    h.fill(systematic=syst, **fields, weight=weights.weight(syst)*cut)
                for syst in shiftSystematics:
                    cut = {s for s in regions[region] if s not in shiftedSelections}
                    cut.update({s+syst for s in regions[region] if s in shiftedSelections})
                    cut = selection.all(cut)
                    for val in shiftedQuantities:
                        fields[val] = df[val+'_'+syst]
                    h.fill(systematic=syst, **fields, weight=weight*cut)
        elif len(region) > 1:
            raise ValueError("Histogram '%s' has a name matching multiple region definitions: %r" % (histname, region))
        else:
            raise ValueError("Histogram '%s' does not fall into any region definitions." % (histname, ))

        hout[histname] = h

    return hout


def collect(output, accumulators):
    for k, v in output.items():
        if k in accumulators:
            try:
                accumulators[k] += v
            except TypeError:
                try:
                    accumulators[k] |= v
                except TypeError:
                    raise
        else:
            raise ValueError("output of processor contains unknown accumulator '%s'" % k)


def processfile(dataset, file):
    fin = uproot.open(file)
    skim_sumw = None
    if 'otree' in fin:
        tree = fin['otree']
        if 'SumWeights' in fin:
            skim_sumw = fin['SumWeights'].values[0]
    else:
        tree = fin['Events']

    tic = time.time()

    output = copy.deepcopy(accumulator_def)
    # would be cool to use columns_accessed and work time to dynamically optimize this
    stride = 500000
    for index in range(tree.numentries//stride + 1):
        df = processor.DataFrame(tree, stride, index)
        df['dataset'] = dataset
        # hacky way to only accumulate file-level information once
        df['skim_sumw'] = skim_sumw if index == 0 else None
        collect(process(df), output)

    toc = time.time()

    output['nentries'] = tree.numentries
    output['bytesread'] = fin.source.bytesread if isinstance(fin.source, uproot.source.xrootd.XRootDSource) else 0
    output['sumworktime'] = toc-tic
    output['columns_accessed'] = df.materialized
    return output


final_accumulators = copy.deepcopy(accumulator_def)

# TODO: rip this out and into fnal_column_analysis_tools.processor
if test:
    nworkers = 1
    testfiles = [
        ("TTToHadronic_TuneCP5_13TeV_powheg_pythia8", "TTToHadronic_TuneCP5_13TeV_powheg_pythia8_0.root"),
        ("TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_0.root"),
        ("JetHT", "JetHTRun2017F_17Nov2017_v1_24.rootnodupl.root"),
        ("SingleMuon", "SingleMuonRun2017B_17Nov2017_v1_2.root"),
    ]
    for i, (dataset, file) in enumerate(testfiles):
        collect(processfile(dataset, file), final_accumulators)
        print("Done processing test file %d" % i)
else:
    nworkers = 12
    with open("metadata/samplefiles.json") as fin:
        samplefiles = json.load(fin)
    #fileslice = slice(None, 5)
    fileslice = slice(None)
    #with concurrent.futures.ThreadPoolExecutor(max_workers=nworkers) as executor:
    with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
        futures = set()
        samples = samplefiles['Hbb_create_2017']
        for group, datasets in samples.items():
            for dataset, files in datasets.items():
                futures.update(executor.submit(processfile, dataset, file) for file in files)
        try:
            total = len(futures)
            processed = 0
            while len(futures) > 0:
                finished = set(job for job in futures if job.done())
                for job in finished:
                    collect(job.result(), final_accumulators)
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

scale = {}
sumw = final_accumulators.pop('sumw').values(overflow='all')
for (ds,), ds_sumw in sumw.items():
    print(ds, ds_sumw)
    scale[ds] = lumi*xsections[ds] / (ds_sumw[1]-ds_sumw[0])

for h in final_accumulators.values():
    if isinstance(h, hist.Hist):
        h.scale(scale, axis="dataset")

print("Columns accessed:")
for col in sorted(list(final_accumulators['columns_accessed'])):
    print("\t", col)
print("%.2f us*cpu/event work time" % (1e6*final_accumulators['sumworktime']/final_accumulators['nentries'], ))
nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in final_accumulators.values() if isinstance(h, hist.Hist))
nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in final_accumulators.values() if isinstance(h, hist.Hist))
print("Processed %.1fM events" % (final_accumulators['nentries']/1e6, ))
print("Read %.1fM bytes" % (final_accumulators['bytesread']/1e6, ))
print("Filled %.1fM bins" % (nbins/1e6, ))
print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

# Pickle is not very fast or memory efficient, will be replaced by something better soon
with gzip.open("hists.pkl.gz", "wb") as fout:
    pickle.dump(final_accumulators, fout)

dt = time.time() - tstart
print("%.2f us*cpu/event overall" % (1e6*dt*nworkers/final_accumulators['nentries'], ))

if test:
    profiler.stop()
    with open("run_baconbits.html", "w") as fout:
        fout.write(profiler.output_html())


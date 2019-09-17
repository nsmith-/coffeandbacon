import uproot
import numpy as np
import json
import subprocess
import awkward
import lz4.frame as lz4f
import hashlib
import os

def make_samples():
    samples = {
        "ttHTobb_M125_TuneCP5_13TeV_powheg_pythia8": {
            'dataset': "/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
        "WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8": {
            'dataset': "/WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
        "WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8": {
            'dataset': "/WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
        "ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8": {
            'dataset': "/ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
        "ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8": {
            'dataset': "/ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
        "ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8": {
            'dataset': "/ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
        "VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix": {
            'dataset': "/VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
        "GluGluHToBB_M125_13TeV_powheg_pythia8": {
            'dataset': "/GluGluHToBB_M125_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM",
        },
    }

    for sample, info in samples.items():
        sfiles = subprocess.getoutput('dasgoclient --query="file dataset=%s"' % info['dataset'])
        info['nano_files'] = sfiles

    with open('metadata/samples_allyears.json') as fin:
        loc = json.load(fin)

    for sample, info in samples.items():
        samplefiles = loc['Hbb_create_2018'][sample]
        info['bits_files'] = samplefiles

    with open('metadata/join_nano.json', 'w') as fout:
        json.dump(samples, fout)


#make_samples()
with open('metadata/join_nano.json') as fin:
    samples = json.load(fin)

newcolumns = {'HTXS_Higgs_pt', 'HTXS_stage1_1_cat_pTjet30GeV', 'HTXS_stage1_1_fine_cat_pTjet30GeV', 'HTXS_njets30', 'FatJet_btagDDBvL', 'event'}

def makeindex(lumi, event):
    lumibits = 18
    if np.any(lumi > (1 << lumibits)):
        raise RuntimeError("lumi too big:", lumi.max())
    if np.any(event > (1 << (64 - lumibits))):
        raise RuntimeError("event too big:", event.max())
    return (lumi.astype(np.uint64) << (64 - lumibits)) + event.astype(np.uint64)

for sample, info in samples.items():
    print(sample)
    samplefiles = info['bits_files']
    nanofiles = info['nano_files'].split('\n')
    for file in samplefiles:
        print("    ", file)
        hashed = hashlib.sha256(file.encode('ascii')).hexdigest()
        ofile = 'adata/' + hashed + '.awkd'
        if os.path.exists(ofile):
            continue
        fbits = uproot.open(file)
        tbits = fbits['otree']
        run = tbits['runNum'].array()
        if not np.all(run == 1):
            print("MC not all run 1??")
            print(np.unique(run, return_counts=True))
        lumi = tbits['lumiSec'].array()
        event = tbits['evtNum'].array()
        index = makeindex(lumi, event)
        idxmap = np.arange(len(index))
        allfound = np.zeros(index.shape, dtype=bool)
        put = np.array([], dtype=idxmap.dtype)
        cols = {}
        for nfile in nanofiles:
            fnano = uproot.open('root://cmsxrootd.fnal.gov/' + nfile)
            tnano = fnano['Events']
            nlumi = tnano['luminosityBlock'].array()
            nevent = tnano['event'].array()
            nindex = makeindex(nlumi, nevent)
            nidxmap = np.argsort(nindex)
            pos = np.searchsorted(nindex[nidxmap], index)
            found = pos < len(nindex)
            found[found] = nindex[nidxmap[pos[found]]] == index[found]
            allfound |= found
            take = nidxmap[pos[found]]
            put = np.append(put, idxmap[found])
            for cname, cval in tnano.arrays(newcolumns, namedecode='ascii').items():
                if cname in cols:
                    cols[cname] = awkward.concatenate([cols[cname], cval[take]])
                else:
                    cols[cname] = cval[take]
            print("        %d %s" % (found.sum(), nfile))
        cols = awkward.Table(cols)
        cols = cols[np.argsort(put)]
        print("allfound: ", allfound.sum() == len(allfound))
        print("event matches:", (cols['event'] == event).all())
        print("saving:", ofile)
        awkward.save(ofile, cols, mode='w')



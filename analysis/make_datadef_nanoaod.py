#!/usr/bin/env python
import uproot
import json
import pyxrootd.client
import fnmatch
import numpy as np
import numexpr
import subprocess
import concurrent.futures
import warnings
import os
import difflib

uproot_xrootd_opts = dict(chunkbytes=30*1024, limitbytes=20*(1024**2))
fnaleos = "root://cmsxrootd.fnal.gov/"

datasets = [
    '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM',
    '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/NANOAODSIM',
    '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM',
    '/DoubleEG/Run2017B-31Mar2018-v1/NANOAOD',
    '/DoubleEG/Run2017C-31Mar2018-v1/NANOAOD',
    '/DoubleEG/Run2017D-31Mar2018-v1/NANOAOD',
    '/DoubleEG/Run2017E-31Mar2018-v1/NANOAOD',
    '/DoubleEG/Run2017F-31Mar2018-v1/NANOAOD',
    '/DoubleMuon/Run2017B-31Mar2018-v1/NANOAOD',
    '/DoubleMuon/Run2017C-31Mar2018-v1/NANOAOD',
    '/DoubleMuon/Run2017D-31Mar2018-v1/NANOAOD',
    '/DoubleMuon/Run2017E-31Mar2018-v1/NANOAOD',
    '/DoubleMuon/Run2017F-31Mar2018-v1/NANOAOD',
]
getentries = False

xsections = {
    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8": 18610.,
    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8": 1921.8*3,
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8": 87.3,
}

xrdfs = pyxrootd.client.FileSystem(fnaleos)

def xrdls(directory, fullpath=True):
    status, listing = xrdfs.dirlist(directory)
    if status['status'] != 0:
        raise Exception("XRootD failed to stat %s" % directory)
    prefix = directory+"/" if fullpath else ""
    return ["%s%s" % (prefix, d['name']) for d in listing['dirlist']]


def xrdfstat(path):
    status, stat = xrdfs.stat(path)
    if status['status'] != 0:
        raise Exception("XRootD failed to stat %s, status=%d" % (path, status))
    return stat


datadef = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    for dataset in datasets:
        flist = subprocess.getoutput('dasgoclient -query="file dataset=%s"' % dataset).split()

        nbytes = executor.map(lambda path: xrdfstat(path)['size'], flist)
        nbytes = np.array(list(nbytes))

        urllist = [fnaleos+path for path in flist]
        nentries = np.array(0)
        if getentries:
            nentries = uproot.numentries(urllist, "Events", total=False, executor=executor)
            nentries = np.array(list(nentries.values()))

        pd = dataset.split('/')[1]

        print(dataset)
        print("    # Files:", len(flist))
        print("    Total bytes: %d" % nbytes.sum())
        print("    Avg. bytes: %.0f" % (nbytes.sum()/nbytes.size))
        if getentries:
            print("    Total entries: %d" % nentries.sum())
            print("    Avg. entries: %.0f" % (nentries.sum()/nentries.size))
            if pd in xsections:
                print("    Effective lumi (assuming weight=1): %.0f /pb" % (nentries.sum()/xsections[pd]))
        if pd in xsections:
            xs = xsections[pd]
        if pd not in xsections and 'Run201' not in dataset:
            nearest = list(xsections.keys())
            nearest.sort(key=lambda s: difflib.SequenceMatcher(None, s, pd).ratio())
            print(pd, " missing xsection, taking closest name:", nearest[-1])
            xs = xsections[nearest[-1]]
        elif 'Run201' in dataset:
            xs = 0.

        if pd in datadef:
            datadef[pd]['files'].extend(urllist)
            datadef[pd]['bytes'] += int(nbytes.sum())
            datadef[pd]['entries'] += int(nentries.sum())
        else:
            datadef[pd] = {
                'files': urllist,
                'xs': xs,
                'bytes': int(nbytes.sum()),
                'entries': int(nentries.sum()),
            }

with open("data/datadef_nano.json", "w") as fout:
    json.dump(datadef, fout, indent=4)


#!/usr/bin/env python
import uproot
import json
import pyxrootd.client
import fnmatch
import numpy as np
import numexpr
import concurrent.futures
import warnings
import os
import difflib

uproot_xrootd_opts = dict(chunkbytes=30*1024, limitbytes=20*(1024**2))
fnaleos = "root://cmseos.fnal.gov/"
dazsle_root = "/eos/uscms/store/user/lpcbacon/dazsle/zprimebits-v15.04"
#patterns = ["DYJetsToQQ*", "WJetsToQQ*"]
patterns = ["ZJetsToQQ*", "WJetsToQQ*"]
#patterns = ["WJetsToLNu*"]
getentries = True

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
        raise Exception("XRootD failed to stat %s" % path)
    return stat


datadef = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    for directory in xrdls(dazsle_root):
        dataset = os.path.basename(directory)
        if not any(fnmatch.fnmatch(dataset, pattern) for pattern in patterns):
            continue
        print(directory)
        flist = fnmatch.filter(xrdls(directory), "*.root")
        if len(flist) == 0:
            print("    NO FILES")
            continue

        nbytes = executor.map(lambda path: xrdfstat(path)['size'], flist)
        nbytes = np.array(list(nbytes))

        urllist = [fnaleos+path for path in flist]
        if getentries:
            nentries = uproot.numentries(urllist, "Events", total=False, executor=executor)
            nentries = np.array(list(nentries.values()))

        print("    # Files:", len(flist))
        print("    Total bytes: %d" % nbytes.sum())
        print("    Avg. bytes: %.0f" % (nbytes.sum()/nbytes.size))
        if dataset in xsections:
            xs = xsections[dataset]
        elif dataset not in xsections and 'Run201' not in dataset:
            nearest = list(xsections.keys())
            nearest.sort(key=lambda s: difflib.SequenceMatcher(None, s, dataset).ratio())
            print("    ", dataset, " missing xsection, taking closest name:", nearest[-1])
            xs = xsections[nearest[-1]]
        elif 'Run201' in dataset:
            xs = 0.
        if getentries:
            print("    Total entries: %d" % nentries.sum())
            print("    Avg. entries: %.0f" % (nentries.sum()/nentries.size))
            print("    Effective lumi (assuming weight=1): %.0f /pb" % (nentries.sum()/xs))
        datadef[dataset] = urllist

samples = {
    'datadef': datadef
}
with open("metadata/datadef.json", "w") as fout:
    json.dump(samples, fout, indent=4)


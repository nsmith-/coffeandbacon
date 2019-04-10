#!/usr/bin/env python
from __future__ import print_function, division
import collections
import warnings
import gzip
import pickle
import json
import time
import cloudpickle

import uproot
import numpy as np
from fnal_column_analysis_tools import hist, processor

# instrument xrootd source
def _read(self, chunkindex):
    self.bytesread = getattr(self, 'bytesread', 0) + self._chunkbytes
    return self._read_real(chunkindex)

uproot.source.xrootd.XRootDSource._read_real = uproot.source.xrootd.XRootDSource._read
uproot.source.xrootd.XRootDSource._read = _read


def process_file(dataset, file, processor_instance, stats_accumulator, preload_items=None):
    fin = uproot.open(file)
    skim_sumw = None
    if 'otree' in fin:
        tree = fin['otree']
        if 'SumWeights' in fin:
            skim_sumw = fin['SumWeights'].values[0]
    else:
        tree = fin['Events']

    tic = time.time()

    output = processor_instance.accumulator.identity()
    # would be cool to use columns_accessed and work time to dynamically optimize this
    stride = 500000
    for index in range(tree.numentries//stride + 1):
        df = processor.DataFrame(tree, stride, index, preload_items=preload_items)
        df['dataset'] = dataset
        # hacky way to only accumulate file-level information once
        df['skim_sumw'] = skim_sumw if index == 0 else None
        output += processor_instance.process(df)

    toc = time.time()

    stats = stats_accumulator.identity()
    stats['nentries'] += tree.numentries
    stats['bytesread'] += fin.source.bytesread if isinstance(fin.source, uproot.source.xrootd.XRootDSource) else 0
    stats['sumworktime'] += toc-tic
    stats['columns_accessed'] += df.materialized
    return output, stats


if __name__ == '__main__':
    test = True
    executor = 'futures'
    # executor = 'iterative'
    nworkers = 1 if executor == 'iterative' else 2
    preload_items = {}
    output_filename = "hists.pkl.gz"
    processor_instance_filename = "boostedHbbProcessor.cpkl.gz"
    profile_html_filename = 'run_baconbits.html'

    if test:
        filelist = [
            ("TTToHadronic_TuneCP5_13TeV_powheg_pythia8", "TTToHadronic_TuneCP5_13TeV_powheg_pythia8_0.root"),
            ("TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_0.root"),
            ("JetHT", "JetHTRun2017F_17Nov2017_v1_24.rootnodupl.root"),
            ("SingleMuon", "SingleMuonRun2017B_17Nov2017_v1_2.root"),
        ]
    else:
        with open("metadata/samplefiles.json") as fin:
            samplefiles = json.load(fin)
        samples = samplefiles['Hbb_create_2017']
        filelist = []
        for group, datasets in samples.items():
            for dataset, files in datasets.items():
                for file in files[:]:
                    filelist.append((dataset, file))


    with gzip.open(processor_instance_filename, "rb") as fin:
        processor_instance = cloudpickle.load(fin)

    combined_accumulator = processor.dict_accumulator({
        'stats': processor.dict_accumulator({
            'nentries': processor.accumulator(0),
            'bytesread': processor.accumulator(0),
            'sumworktime': processor.accumulator(0.),
            'columns_accessed': processor.set_accumulator(),
        }),
        'job': processor_instance.accumulator.identity(),
    })

    def work_function(item):
        dataset, file = item
        out, stats = process_file(dataset, file, processor_instance, combined_accumulator['stats'], preload_items)
        return processor.dict_accumulator({'stats': stats, 'job': out})

    tstart = time.time()
    if executor == 'iterative':
        from pyinstrument import Profiler
        profiler = Profiler()
        profiler.start()
        processor.iterative_executor(filelist, work_function, combined_accumulator)
        profiler.stop()
        with open(profile_html_filename, "w") as fout:
            fout.write(profiler.output_html())
    elif executor == 'futures':
        processor.futures_executor(filelist, work_function, combined_accumulator, workers=nworkers)

    final_accumulator = combined_accumulator['job']
    stats = combined_accumulator['stats']
    processor_instance.postprocess(final_accumulator)

    print("Columns accessed:", set(stats['columns_accessed']))
    print("%.2f us*cpu/event work time" % (1e6*stats['sumworktime'].value/stats['nentries'].value, ))
    print("Processed %.1fM events" % (stats['nentries'].value/1e6, ))
    print("Read %.1fM bytes" % (stats['bytesread'].value/1e6, ))

    nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in final_accumulator.values() if isinstance(h, hist.Hist))
    nfilled = sum(sum(np.sum(arr>0) for arr in h._sumw.values()) for h in final_accumulator.values() if isinstance(h, hist.Hist))
    print("Filled %.1fM bins" % (nbins/1e6, ))
    print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

    # Pickle is not very fast or memory efficient, will be replaced by something better soon
    with gzip.open(output_filename, "wb") as fout:
        pickle.dump(final_accumulator, fout)

    dt = time.time() - tstart
    print("%.2f us*cpu/event overall" % (1e6*dt*nworkers/stats['nentries'].value, ))

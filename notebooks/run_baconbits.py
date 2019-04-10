#!/usr/bin/env python
import gzip
import pickle
import json
import time
import cloudpickle
import argparse

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
    parser = argparse.ArgumentParser(description='Run analysis on baconbits files using processor cloudpickle files')
    parser.add_argument('--processor', default='boostedHbbProcessor.cpkl.gz', help='The name of the compiled processor file')
    parser.add_argument('--output', default='hists.pkl.gz', help='Output histogram filename')
    parser.add_argument('--samplejson', default='metadata/samplefiles.json', help='JSON file containing dataset and file locations')
    parser.add_argument('--sample', default='Hbb_create_2017', help='The sample to use in the sample JSON')
    parser.add_argument('--limit', type=int, default=None, metavar='N', help='Limit to the first N files of each dataset in sample JSON')
    parser.add_argument('--test', action='store_true', help='Process the test files defined in this execuatble')
    parser.add_argument('--executor', choices=['iterative', 'futures'], default='iterative', help='The type of executor to use')
    parser.add_argument('-j', '--workers', type=int, default=12, help='Number of workers to use for multi-worker executors (e.g. futures or condor)')
    parser.add_argument('--profile-out', dest='profilehtml', default=None, help='Filename for the pyinstrument HTML profile output')
    args = parser.parse_args()

    # Set a list of preloaded columns, to profile the execution separately from the uproot deserialization
    preload_items = {}

    if args.test:
        filelist = [
            ("TTToHadronic_TuneCP5_13TeV_powheg_pythia8", "TTToHadronic_TuneCP5_13TeV_powheg_pythia8_0.root"),
            ("TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8", "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_0.root"),
            ("JetHT", "JetHTRun2017F_17Nov2017_v1_24.rootnodupl.root"),
            ("SingleMuon", "SingleMuonRun2017B_17Nov2017_v1_2.root"),
        ]
    else:
        with open(args.samplejson) as fin:
            samplefiles = json.load(fin)
        samples = samplefiles[args.sample]
        filelist = []
        for group, datasets in samples.items():
            for dataset, files in datasets.items():
                for file in files[:args.limit]:
                    filelist.append((dataset, file))

    with gzip.open(args.processor, "rb") as fin:
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
    if args.executor == 'iterative':
        if args.profilehtml is not None:
            from pyinstrument import Profiler
            profiler = Profiler()
            profiler.start()
        processor.iterative_executor(filelist, work_function, combined_accumulator)
        if args.profilehtml is not None:
            profiler.stop()
            with open(args.profilehtml, "w") as fout:
                fout.write(profiler.output_html())
    elif args.executor == 'futures':
        processor.futures_executor(filelist, work_function, combined_accumulator, workers=args.workers)

    final_accumulator = combined_accumulator['job']
    stats = combined_accumulator['stats']
    processor_instance.postprocess(final_accumulator)

    print("Columns accessed:", set(stats['columns_accessed']))
    print("%.2f us*cpu/event work time" % (1e6*stats['sumworktime'].value/stats['nentries'].value, ))
    print("Processed %.1fM events" % (stats['nentries'].value/1e6, ))
    print("Read %.1fM bytes" % (stats['bytesread'].value/1e6, ))

    nbins = sum(sum(arr.size for arr in h._sumw.values()) for h in final_accumulator.values() if isinstance(h, hist.Hist))
    nfilled = sum(sum(np.sum(arr > 0) for arr in h._sumw.values()) for h in final_accumulator.values() if isinstance(h, hist.Hist))
    print("Filled %.1fM bins" % (nbins/1e6, ))
    print("Nonzero bins: %.1f%%" % (100*nfilled/nbins, ))

    # Pickle is not very fast or memory efficient, will be replaced by something better soon
    with gzip.open(args.output, "wb") as fout:
        pickle.dump(final_accumulator, fout)

    dt = time.time() - tstart
    nworkers = 1 if args.executor == 'iterative' else args.workers
    print("%.2f us*cpu/event overall" % (1e6*dt*nworkers/stats['nentries'].value, ))

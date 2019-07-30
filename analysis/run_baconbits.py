#!/usr/bin/env python
import json
import time
import argparse
import glob
from tqdm import tqdm

import uproot
import numpy as np
from coffea import hist, processor
from coffea.util import load, save


# instrument xrootd source
if not hasattr(uproot.source.xrootd.XRootDSource, '_read_real'):
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
        df = processor.LazyDataFrame(tree, stride, index, preload_items=preload_items)
        df['dataset'] = dataset
        # hacky way to only accumulate file-level information once
        if 'otree' in fin:
            df['skim_sumw'] = skim_sumw if index == 0 else None
        output += processor_instance.process(df)

    toc = time.time()

    stats = stats_accumulator.identity()
    stats['nentries'] += tree.numentries
    stats['bytesread'] += fin.source.bytesread if isinstance(fin.source, uproot.source.xrootd.XRootDSource) else 0
    stats['sumworktime'] += toc-tic
    stats['columns_accessed'] += df.materialized
    return output, stats


def validate(dataset, file):
    fin = uproot.open(file)
    if 'otree' in fin:
        return fin['otree'].numentries
    else:
        return fin['Events'].numentries


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run analysis on baconbits files using processor coffea files')
    parser.add_argument('--processor', default='boostedHbbProcessor.coffea', help='The name of the compiled processor file (default: %(default)s)')
    parser.add_argument('--output', default=r'hists_{sample}.coffea', help='Output histogram filename (default: %(default)s)')
    parser.add_argument('--samplejson', default='metadata/samplefiles.json', help='JSON file containing dataset and file locations (default: %(default)s)')
    parser.add_argument('--sample', default='test_skim', help='The sample to use in the sample JSON (default: %(default)s)')
    parser.add_argument('--limit', type=int, default=None, metavar='N', help='Limit to the first N files of each dataset in sample JSON')
    parser.add_argument('--validate', action='store_true', help='Do not process, just check all files are accessible')
    parser.add_argument('--executor', choices=['iterative', 'futures'], default='iterative', help='The type of executor to use (default: %(default)s)')
    parser.add_argument('-j', '--workers', type=int, default=12, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    parser.add_argument('--profile-out', dest='profilehtml', default=None, help='Filename for the pyinstrument HTML profile output')
    args = parser.parse_args()
    if args.output == parser.get_default('output'):
        args.output = 'hists_%s.coffea' % args.sample

    # Set a list of preloaded columns, to profile the execution separately from the uproot deserialization
    preload_items = {}

    with open(args.samplejson) as fin:
        samplefiles = json.load(fin)
    if args.sample not in samplefiles:
        print("Sample '%s' not available in %s, available:" % (args.sample, args.samplejson), list(samplefiles.keys()))
    sample = samplefiles[args.sample]
    filelist = []
    for dataset, files in sample.items():
        if type(files) is str:
            if files.startswith('root://'):
                raise ValueError("glob not supported over xrootd (yet)")
            prefix = 'root://cmseos.fnal.gov/' if files.startswith('/eos/uscms/') else ''
            files = [prefix + f for f in glob.glob(files)]
        for file in files[:args.limit]:
            filelist.append((dataset, file))

    if args.validate:
        for ds, fn in tqdm(filelist, desc='Validating files'):
            try:
                validate(ds, fn)
            except OSError:
                print("File open error for %s, %s" % (ds, fn))
        exit(0)

    processor_instance = load(args.processor)

    combined_accumulator = processor.dict_accumulator({
        'stats': processor.dict_accumulator({
            'nentries': processor.value_accumulator(int),
            'bytesread': processor.value_accumulator(int),
            'sumworktime': processor.value_accumulator(float),
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

    save(final_accumulator, args.output)

    dt = time.time() - tstart
    nworkers = 1 if args.executor == 'iterative' else args.workers
    print("%.2f us*cpu/event overall" % (1e6*dt*nworkers/stats['nentries'].value, ))

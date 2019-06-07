#!/usr/bin/env python
from __future__ import print_function, division
import lz4.frame as lz4f
import cloudpickle
import json
import argparse

import uproot
import numpy as np
from coffea import processor


def get_pileup(item):
    dataset, filename = item
    file = uproot.open(filename)
    puhist = file["Pu"]
    pileup = processor.accumulator(np.zeros_like(puhist.values))
    pileup += puhist.values
    sumwhist = file["SumWeights"]
    sumw = processor.accumulator(np.zeros(1))
    sumw += sumwhist.values[0]
    return processor.dict_accumulator({
        'pileup': processor.dict_accumulator({dataset: pileup}),
        'sumw': processor.dict_accumulator({dataset: sumw}),
    })


def make_pileup(args):
    with open(args.samplejson) as fin:
        samplefiles = json.load(fin)
    sample = samplefiles[args.sample]

    filelist = []
    for dataset, files in sample.items():
        if dataset == 'JetHT' or dataset == 'SingleMuon':
            continue
        for file in files:
            filelist.append((dataset, file))

    final_accumulator = processor.dict_accumulator({
        'pileup': processor.dict_accumulator(),
        'sumw': processor.dict_accumulator(),
    })
    processor.futures_executor(filelist, get_pileup, final_accumulator, workers=args.workers)

    with lz4f.open("correction_files/pileup_mc.cpkl.lz4", "wb") as fout:
        cloudpickle.dump(final_accumulator['pileup'], fout)

    with lz4f.open("correction_files/sumw_mc.cpkl.lz4", "wb") as fout:
        cloudpickle.dump(final_accumulator['sumw'], fout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Precompute pileup correction for a given set of samples')
    parser.add_argument('--samplejson', default='metadata/samplefiles.json', help='JSON file containing dataset and file locations')
    parser.add_argument('--sample', default='Hbb_2017', help='The sample to use in the sample JSON')
    parser.add_argument('-j', '--workers', type=int, default=8, help='Number of workers to use for multi-worker executors (e.g. futures or condor)')
    args = parser.parse_args()

    make_pileup(args)

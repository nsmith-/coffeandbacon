#!/usr/bin/env python
from __future__ import print_function, division
import json
import argparse
from functools import partial

import uproot
import numpy as np
from coffea import processor
from coffea.util import load, save


def get_pileup(item):
    dataset, filename = item
    file = uproot.open(filename)
    puhist = file["Pu"]
    pileup = processor.value_accumulator(partial(np.zeros, puhist.values.size))
    pileup += puhist.values
    sumwhist = file["SumWeights"]
    sumw = processor.value_accumulator(int)
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

    save(final_accumulator['pileup'], 'correction_files/pileup_mc.coffea')
    save(final_accumulator['sumw'], 'correction_files/sumw_mc.coffea')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Precompute MC pileup distribution for a given set of samples (mainly useful for 2017)')
    parser.add_argument('--samplejson', default='metadata/samplefiles.json', help='JSON file containing dataset and file locations (default: %(default)s)')
    parser.add_argument('--sample', default='Hbb_2017', help='The sample to use in the sample JSON (default: %(default)s)')
    parser.add_argument('-j', '--workers', type=int, default=8, help='Number of workers to use for multi-worker executors (e.g. futures or condor) (default: %(default)s)')
    args = parser.parse_args()

    make_pileup(args)

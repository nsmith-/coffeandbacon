#!/usr/bin/env python
from __future__ import print_function, division
import lz4.frame as lz4f
import cloudpickle
import json
import sys

import uproot
import numpy as np
from fnal_column_analysis_tools import processor

with open("metadata/samplefiles.json") as fin:
    samplefiles = json.load(fin)
sample = samplefiles[sys.argv[-1]]

filelist = []
for dataset, files in sample.items():
    if dataset == 'JetHT' or dataset == 'SingleMuon':
        continue
    for file in files:
        filelist.append((dataset, file))


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


final_accumulator = processor.dict_accumulator({
    'pileup': processor.dict_accumulator(),
    'sumw': processor.dict_accumulator(),
})
processor.futures_executor(filelist, get_pileup, final_accumulator, workers=8)

with lz4f.open("correction_files/pileup_mc.cpkl.lz4", "wb") as fout:
    cloudpickle.dump(final_accumulator['pileup'], fout)

with lz4f.open("correction_files/sumw_mc.cpkl.lz4", "wb") as fout:
    cloudpickle.dump(final_accumulator['sumw'], fout)


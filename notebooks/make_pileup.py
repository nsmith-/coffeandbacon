#!/usr/bin/env python
from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import warnings
import concurrent.futures
import gzip
import pickle
import json
import time

import uproot
import numpy as np

with open("metadata/datadef.json") as fin:
    datadef = json.load(fin)

def get_pileup(dataset, file):
    pu = uproot.open(file)["Pu"]
    return dataset, pu.values

pileup_mc = {}
nworkers = 10
fileslice = slice(None)
with concurrent.futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
    futures = set()
    for dataset, info in datadef.items():
        futures.update(executor.submit(get_pileup, dataset, file) for file in info['files'][fileslice])
    try:
        total = len(futures)
        processed = 0
        while len(futures) > 0:
            finished = set(job for job in futures if job.done())
            for job in finished:
                dataset, hpu = job.result()
                if dataset in pileup_mc:
                    pileup_mc[dataset] += hpu
                else:
                    pileup_mc[dataset] = hpu
                processed += 1
                if processed % 10 == 0:
                    print("Processing: done with % 4d / % 4d files" % (processed, total))
            futures -= finished
        del finished
    except KeyboardInterrupt:
        print("Ok quitter")
        for job in futures: job.cancel()
    except:
        for job in futures: job.cancel()
        raise


with gzip.open("pileup_mc.pkl.gz", "wb") as fout:
    pickle.dump(pileup_mc, fout)


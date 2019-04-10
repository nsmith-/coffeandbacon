# coding: utf-8
from __future__ import print_function, division
from collections import defaultdict
import gzip
import pickle
import json
import re
import os

import uproot
import numpy as np

from fnal_column_analysis_tools import hist
from fnal_column_analysis_tools.hist import export
import processmap

with gzip.open("hists.pkl.gz") as fin:
    hists_unmapped = pickle.load(fin)


hists = {}
for key, val in hists_unmapped.items():
    if isinstance(val, hist.Hist):
        hists[key] = processmap.apply(val)


if os.path.exists("templates.root"):
    os.remove("templates.root")
fout = uproot.create("templates.root")

nodata = re.compile("(?!data_obs)")
h = hists['templates_signalregion'][nodata]
lumi = 41.1
h.scale({p: lumi for p in h[nodata].identifiers('process')}, axis="process")

for proc in h.identifiers('process'):
    for i, ptbin in enumerate(h.identifiers('AK8Puppijet0_pt')):
        for syst in h.identifiers('systematic'):
            mproj = (slice(None), 'all')
            systreal = syst
            fail_template = (h.project('process', proc)
                              .project('AK8Puppijet0_isHadronicV', *mproj)
                              .project('systematic', systreal)
                              .project('AK8Puppijet0_pt', ptbin)
                              .project('AK8Puppijet0_deepdoubleb', slice(None,0.89), overflow='under')
                            )
            pass_template = (h.project('process', proc)
                              .project('AK8Puppijet0_isHadronicV', *mproj)
                              .project('systematic', systreal)
                              .project('AK8Puppijet0_pt', ptbin)
                              .project('AK8Puppijet0_deepdoubleb', slice(0.89,None))
                            )
            content = fail_template.sum('AK8Puppijet0_msd').values()
            if content == {} or content[()] == 0.:
                print(proc, ptbin, syst)
                continue
            sname = "_%s" % syst if syst != '' else ''
            name = "%s_pass%s_bin%d" % (proc, sname, i)
            fout[name] = export.export1d(pass_template)
            name = "%s_fail%s_bin%d" % (proc, sname, i)
            fout[name] = export.export1d(fail_template)

fout.close()

if os.path.exists("hist_1DZbb_muonCR.root"):
    os.remove("hist_1DZbb_muonCR.root")
fout = uproot.create("hist_1DZbb_muonCR.root")

h = hists['templates_muoncontrol']
lumi = 41.1
h.scale({p: lumi for p in h[nodata].identifiers('process')}, axis="process")

rename = {
    'trigweight': 'trigger',
    'pileupweight': 'Pu',
    'mutrigweight': 'mutrigger',
    'muidweight': 'muid',
    'muisoweight': 'muiso',
    'matchedUp': 'matched',
    'matchedDown': 'unmatched',
}

for proc in h.identifiers('process'):
    for syst in h.identifiers('systematic'):
        mproj = (slice(None), 'all')
        systreal = syst
        fail_template = (h.project('process', proc)
                            .project('AK8Puppijet0_isHadronicV', *mproj)
                            .project('systematic', systreal)
                            .project('AK8Puppijet0_pt', overflow='all')
                            .project('AK8Puppijet0_deepdoubleb', slice(None,0.89), overflow='under')
                        )
        pass_template = (h.project('process', proc)
                            .project('AK8Puppijet0_isHadronicV', *mproj)
                            .project('systematic', systreal)
                            .project('AK8Puppijet0_pt', overflow='all')
                            .project('AK8Puppijet0_deepdoubleb', slice(0.89,None))
                        )
        content = fail_template.sum('AK8Puppijet0_msd').values()
        if content == {} or content[()] == 0.:
            print(proc, syst)
            continue
        sname = "_%s" % syst if syst != '' else ''
        for k,v in rename.items():
            sname = sname.replace(k, v)
        name = "%s_pass%s" % (proc, sname)
        fout[name] = export.export1d(pass_template)
        name = "%s_fail%s" % (proc, sname)
        fout[name] = export.export1d(fail_template)

fout.close()

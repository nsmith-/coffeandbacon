    # coding: utf-8
from __future__ import print_function, division
from collections import defaultdict
import gzip
import json
import re
import os

import uproot
import numpy as np

from coffea import hist
from coffea.util import load, save
import processmap

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cc", default=False, action='store_true', help="Make templates for Hcc")
parser.add_argument("--split", default=False, action='store_true', help="Split W/Z by flavor")
parser.add_argument("--three-regions", dest='threeregions', default=False, action='store_true', help="Split by max prob pbb/pcc/pqq")
args = parser.parse_args()

hists_unmapped = load('hists_Hbb_create_2017.coffea')


hists = {}
for key, val in hists_unmapped.items():
    if isinstance(val, hist.Hist):
        hists[key] = processmap.apply(val)

if args.cc: template_file = "templatesCC.root"
elif args.threeregions: template_file = "templates3.root"
else: template_file = "templates.root"
if os.path.exists(template_file):
    os.remove(template_file)
fout = uproot.create(template_file)


if args.cc: h = hists['templates_hCCsignalregion']
elif args.threeregions: h = hists['templates_signalregion']
else: h = hists['templates_signalregion']

# Scale MC
nodata = re.compile("(?!data_obs)")
lumi = 41.1
h.scale({p: lumi for p in h[nodata].identifiers('process')}, axis="process")

proc_names = h.identifiers('process')
if args.split: proc_names += ['wcq', 'zbb', 'zcc']

for proc in proc_names:
    print(proc)
    for i, ptbin in enumerate(h.identifiers('AK8Puppijet0_pt')):
        for syst in h.identifiers('systematic'):
            source_proc = proc
            if args.split:
                if proc == 'zbb':
                    mproj = (3,)
                elif proc == 'wcq' or proc == 'zcc':
                    mproj = (2,)
                elif proc == 'wqq' or proc == 'zqq':
                    mproj = (1,)
                else:
                    mproj = (slice(None), 'all')
                if proc in ['wcq', 'wqq']: source_proc = 'wqq'
                if proc in ['zbb', 'zcc', 'zqq']: source_proc = 'zqq'
            else:
                mproj = (slice(None), 'all')
            systreal = syst
            if args.threeregions:
                pqq_template = (h.integrate('process', source_proc)
                                .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                .integrate('systematic', systreal)
                                .integrate('AK8Puppijet0_pt', ptbin)
                                .integrate('pxx', 1)
                                .integrate('AK8Puppijet0_deepdoubleb', overflow='all')
                                )
                pcc_template = (h.integrate('process', source_proc)
                                .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                .integrate('systematic', systreal)
                                .integrate('AK8Puppijet0_pt', ptbin)
                                .integrate('pxx', 2)
                                .integrate('AK8Puppijet0_deepdoubleb', overflow='all')
                                )
                pbb_template = (h.integrate('process', source_proc)
                                .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                .integrate('systematic', systreal)
                                .integrate('AK8Puppijet0_pt', ptbin)
                                .integrate('pxx', 3)
                                .integrate('AK8Puppijet0_deepdoubleb', overflow='all')
                                )
            elif args.cc:
                fail_template = (h.integrate('process', source_proc)
                                  .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                  .integrate('systematic', systreal)
                                  .integrate('AK8Puppijet0_pt', ptbin)
                                  .integrate('pxx')
                                  .integrate('AK8Puppijet0_deepdoublec', slice(None,0.83), overflow='under')
                                 )
                pass_template = (h.integrate('process', source_proc)
                                  .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                  .integrate('systematic', systreal)
                                  .integrate('AK8Puppijet0_pt', ptbin)
                                  .integrate('pxx')
                                  .integrate('AK8Puppijet0_deepdoublec', slice(0.83,None))
                                 )

            else:
                fail_template = (h.integrate('process', source_proc)
                                  .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                  .integrate('systematic', systreal)
                                  .integrate('AK8Puppijet0_pt', ptbin)
                                  .integrate('pxx')
                                  .integrate('AK8Puppijet0_deepdoubleb', slice(None,0.89), overflow='under')
                                 )
                pass_template = (h.integrate('process', source_proc)
                                  .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                  .integrate('systematic', systreal)
                                  .integrate('AK8Puppijet0_pt', ptbin)
                                  .integrate('pxx')
                                  .integrate('AK8Puppijet0_deepdoubleb', slice(0.89,None))
                                 )
            try:
                content = fail_template.sum('AK8Puppijet0_msd').values()
            except:
                content = pqq_template.sum('AK8Puppijet0_msd').values()
            if content == {} or content[()] == 0.:
                print("Missing", proc, ptbin, syst)
                continue

            sname = "_%s" % syst if syst.name != '' else ''
            if args.threeregions:
                name = "%s_pqq%s_bin%d" % (proc, sname, i)
                fout[name] = hist.export1d(pqq_template)
                name = "%s_pcc%s_bin%d" % (proc, sname, i)
                fout[name] = hist.export1d(pcc_template)
                name = "%s_pbb%s_bin%d" % (proc, sname, i)
                fout[name] = hist.export1d(pbb_template)
            else:
                name = "%s_pass%s_bin%d" % (proc, sname, i)
                fout[name] = hist.export1d(pass_template)
                name = "%s_fail%s_bin%d" % (proc, sname, i)
                fout[name] = hist.export1d(fail_template)

fout.close()

if args.cc: outmu_file = "hist_1DZcc_muonCR.root"
else: outmu_file = "hist_1DZbb_muonCR.root"

if os.path.exists(outmu_file):
    os.remove(outmu_file)
fout = uproot.create(outmu_file)

if args.cc: h = hists['templates_hCCmuoncontrol']
else: h = hists['templates_muoncontrol']
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
        if args.cc:
            fail_template = (h.integrate('process', proc)
                                .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                .integrate('systematic', systreal)
                                .integrate('AK8Puppijet0_pt', overflow='all')
                                .integrate('pxx')
                                .integrate('AK8Puppijet0_deepdoublec', slice(None,0.83), overflow='under')
                            )
            pass_template = (h.integrate('process', proc)
                                .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                .integrate('systematic', systreal)
                                .integrate('AK8Puppijet0_pt', overflow='all')
                                .integrate('pxx')
                                .integrate('AK8Puppijet0_deepdoublec', slice(0.83,None))
                            )
        else:
            fail_template = (h.integrate('process', proc)
                                .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                .integrate('systematic', systreal)
                                .integrate('AK8Puppijet0_pt', overflow='all')
                                .integrate('pxx')
                                .integrate('AK8Puppijet0_deepdoubleb', slice(None,0.89), overflow='under')
                            )
            pass_template = (h.integrate('process', proc)
                                .integrate('AK8Puppijet0_isHadronicV', *mproj)
                                .integrate('systematic', systreal)
                                .integrate('AK8Puppijet0_pt', overflow='all')
                                .integrate('pxx')
                                .integrate('AK8Puppijet0_deepdoubleb', slice(0.89,None))
                            )
        content = fail_template.sum('AK8Puppijet0_msd').values()
        if content == {} or content[()] == 0.:
            print(proc, syst)
            continue
        sname = "_%s" % syst if syst != '' else ''
        for k,v in rename.items():
            sname = sname.replace(k, v)
        name = "%s_pass%s" % (proc, sname)
        fout[name] = hist.export1d(pass_template)
        name = "%s_fail%s" % (proc, sname)
        fout[name] = hist.export1d(fail_template)

fout.close()

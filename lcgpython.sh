#!/usr/bin/env bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc62-opt/setup.sh

xrdcp -s root://cmseos.fnal.gov//store/user/$(whoami)/pylocal.tgz .
tar -zxf pylocal.tgz
export PYTHONPATH=$(find site-packages/ -name *.egg |tr '\n' ':')$PYTHONPATH
echo $@
python $@ || exit $?

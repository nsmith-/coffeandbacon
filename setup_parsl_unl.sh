!/usr/bin/env bash

echo 'source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-slc6-gcc62-opt/setup.sh' >> ${HOME}/.bash_profile
echo 'export PATH=${HOME}/.local/bin:$PATH' >> ${HOME}/.bash_profile
echo 'export PYTHONPATH=${HOME}/.local/lib/python3.6/site-packages:$PYTHONPATH' >> ${HOME}/.bash_profile

source ${HOME}/.bash_profile

pip install --upgrade pip --user
pip install entrypoints --upgrade --user
pip install https://github.com/Parsl/parsl/zipball/master --user
pip install numpy scipy numba pycairo matplotlib pandas --upgrade --user
pip install lz4 cloudpickle --upgrade --user
pip install uproot fnal-column-analysis-tools --upgrade --user

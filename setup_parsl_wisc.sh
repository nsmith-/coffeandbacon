!/usr/bin/env bash

echo 'source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc7-opt/setup.sh' >> ${HOME}/.bashrc
echo 'export PATH=${HOME}/.local/bin:$PATH' >> ${HOME}/.bashrc
echo 'export PYTHONPATH=${HOME}/.local/lib/python3.6/site-packages:$PYTHONPATH' >> ${HOME}/.bashrc

source ${HOME}/.bashrc

pip install --upgrade pip --user
pip install entrypoints --upgrade --user
pip install https://github.com/Parsl/parsl/zipball/master --user
pip install numpy scipy numba pycairo matplotlib pandas --upgrade --user
pip install lz4 cloudpickle --upgrade --user
pip install uproot fnal-column-analysis-tools --upgrade --user

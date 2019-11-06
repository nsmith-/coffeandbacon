# Coffea and bacon
A repository for using columnar analysis on bacon(bits) files.

## Setup
To use, first an environment must be setup.  If running on your laptop, running
`pip install coffea jupyter` should get most of the way there.
It might be preferable to install some of these things via `conda`.  Note, this is
a python 3 package; python 2 is end of life so get used to `print()`.
To install at the LPC, there is a script:
```bash
source setup_lcg.sh
```
On future use at LPC, run `source env_lcg.sh`.


### Setup using conda
#### Install conda if you don't have it
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Install the environment (might take a while, contains its own ROOT)
```bash
conda-env create -n coffea -f=coffea_env.yml
conda activate coffea
```

## Running the Hbb analysis on baconbits
The following recipe runs all the relevant code to produce templates similar to those of `sampleContainer`:
```bash
cd analysis
./make_pileup.py
./compile_corrections.py
./boostedHbbProcessor.py
./run_baconbits.py --executor futures --sample Hbb_2017
python baconbits-templates.py
python convert2d.py
ls hist_1DZbb*
```
#### Running Hcc at DESY
```bash
cd analysis
./make_pileup.py --samplejson samplesdesy.json --sample Hbb_create_2017
./compile_corrections.py
###
# If runnning futures or iterative
./boostedHbbProcessor.py
./run_baconbits.py --executor futures --samplejson samplesdesy.json --sample Hbb_create_2017 -j 10

# If runnning with parsl
./boostedHbbProcessor.py --externalSumW correction_files/sumw_mc.coffea
./run_baconbits.py --executor parsl --samplejson samplesdesy.json --sample Hbb_create_2017
###
python baconbits-templates.py --cc --split
python convert2d.py --cc --split
ls hist_1DZcc*
```

To just get your feet wet, look at `./run_baconbits.py --help`, then run
```bash
./download_testbits.sh
./run_baconbits.py --sample test_bits
```
which will not take much time.

## Plots with jupyter
To use jupyter on a laptop, `jupyter notebook` should work.
To use juptyer at LPC, start a notebook server with, e.g. `jupyter notebook --no-browser --port $PORT`,
substituing your favorite port: pick a random integer in (8000,65535).
Often, you'll need an ssh tunnel, which can be accomplished via, e.g. `ssh -L $PORT:localhost:$PORT server.address`

Check out some of the notebooks in the analysis directory.  Most will need the `hists.coffea` file as created by `run_baconbits.py`,
but the top background notebook can be run after downloading the test baconbits.

# Coffea and bacon
A repository for using columnar analysis on bacon(bits) files.

## Setup
To use, first an environment must be setup.  If running on your laptop, running
`pip install fnal-column-analysis-tools jupyter` should get most of the way there.
It might be preferable to install some of these things via `conda`.  Note, this is
a python 3 package, python 2 is end of life so get used to `print()`.
To install at the LPC, there is a script:
```bash
source setup_lcg.sh
```
On future use at LPC, run `source env_lcg.sh`.

## Running the Hbb analysis
The following recipe runs all the relevant code to produce templates similar to those of `sampleContainer`:
```python
cd analysis
./compile_corrections.py
./boostedHbbProcessor.py
./run_baconbits.py --executor futures
python baconbits-templates.py
python convert2d.py
ls hist_1DZbb*
```

## Plots with jupyter
To use jupyter on a laptop, `jupyter notebook` should work.
To use juptyer at LPC, start a notebook server with, e.g. `jupyter notebook --no-browser --port $PORT`,
substituing your favorite port: pick a random integer in (8000,65535).
Often, you'll need an ssh tunnel, which can be accomplished via, e.g. `ssh -L $PORT:localhost:$PORT server.address`

Check out some of the notebooks in the analysis directory.

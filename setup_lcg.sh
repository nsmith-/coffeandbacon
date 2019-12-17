#!/usr/bin/env bash

source env_lcg.sh

pip install --user coffea
pip install --user tqdm
pip install --user pycairo

# progressbar, sliders, etc.
jupyter nbextension enable --py widgetsnbextension


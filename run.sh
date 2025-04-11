#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate base

python3 ./Tmid_from_R.py
python3 ./SED.py
python3 ./flux_from_time.py

conda deactivate
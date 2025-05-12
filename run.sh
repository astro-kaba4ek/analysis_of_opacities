#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate base

echo "kappaNu_from_lambda"
python3 ./kappaNu_from_lambda.py

echo "kappaPR_from_T"
python3 ./kappaPR_from_T.py

echo "tauPR_from_R"
python3 ./tauPR_from_R.py

echo "tHC_from_R"
python3 ./tHC_from_R.py

echo "delta_tHC_from_R"
python3 ./delta_tHC_from_R.py

echo "Tmid_from_R"
python3 ./Tmid_from_R.py

echo "SED"
python3 ./SED.py

echo "flux_from_time"
python3 ./flux_from_time.py

conda deactivate
#! /bin/bash

./sfvp_t.py
./spectra_analysis_scalar_01Oct -field 4 -mt 49 sfvpData-0-49.nc
./spectra_analysis_scalar_01Oct -field 5 -mt 49 sfvpData-0-49.nc
./spectra-zono-temp-01Oct.py


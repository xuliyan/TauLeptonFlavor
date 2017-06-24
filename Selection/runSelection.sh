#!/bin/bash

# output ntuple directory
NTUPDIR=/tmp/arapyan/Flat

# integrated luminosity for data
LUMI=2215

root -l -q select3Mu.C+\(\"samples.conf\",\"${NTUPDIR}/3Mu\",0\)

rm *.so *.d

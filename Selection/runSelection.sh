#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/user/x/xuyan/3MuonProj/CMSSW_8_0_27/src/DataFlat

# integrated luminosity for data
LUMI=2215

root -l -q selectData.C+\(\"samples.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d

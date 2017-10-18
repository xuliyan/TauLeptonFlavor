#!/bin/bash

# output ntuple directory
NTUPDIR=/afs/cern.ch/user/x/xuyan/3MuonProj/CMSSW_8_0_27/src/GENFlat

# integrated luminosity for data
LUMI=2215

root -l -q select3Mu.C+\(\"samplesGEN.conf\",\"${NTUPDIR}\",0\)

rm *.so *.d

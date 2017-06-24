# TauLeptonFlavor

scram project CMSSW_8_0_27
cd CMSSW_8_0_27/src/
git clone git@github.com:arapyan/TauLeptonFlavor.git 
git clone git@github.com:MiT-HEP/BaconAna.git
scram b
cp /afs/cern.ch/work/a/arapyan/public/HLT_50nsGRun BaconAna/DataFormats/data/

To run the selection:
cd Selection
./runSelection.sh

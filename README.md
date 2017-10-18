# TauLeptonFlavor

scram project CMSSW_8_0_27
cd CMSSW_8_0_27/src/
git clone git@github.com:xuliyan/TauLeptonFlavor.git 
git clone git@github.com:MiT-HEP/BaconAna.git
scram b
cp /afs/cern.ch/work/a/arapyan/public/HLT_50nsGRun BaconAna/DataFormats/data/

To run the selection:
For MC, copy scripts from /Scripts/Tau3MuMC and rename them to select3Mu.C
cd Selection
./runSelectionGEN.sh
For Normalization, copy scripts from /Scripts/Ds2MuPiData and rename them to selectDsPhiPi.C
cd Selection
./runSelection.sh

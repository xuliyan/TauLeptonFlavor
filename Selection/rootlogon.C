{    
  if(gSystem->Getenv("CMSSW_VERSION")) {    
  
    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
    
    gROOT->Macro("$CMSSW_BASE/src/BaconAna/macros/setRootEnv.C+");
   }
  //gROOT->Macro("../Utils/CPlot.cc++");
  //gROOT->Macro("../Utils/MitStyleRemix.cc++");

  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}

//================================================================================================
//
// Select Z->mumu candidates
//
//  * outputs ROOT files of events passing selection
//
//________________________________________________________________________________________________
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "TLorentzVector.h"         // 4-vector class
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TVector3.h>
#include <TF1.h>
#include "TH1D.h"
#include "TRandom.h"
#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
// lumi section selection with JSON files
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"
#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

//=== MAIN MACRO ================================================================================================= 

void selectData(const TString conf="samples.conf", // input file
               const TString outputDir=".",   // output directory
	       const Bool_t  doScaleCorr=0,    // apply energy scale corrections
	       const Bool_t  doPU=0
	       ) {
  gBenchmark->Start("selectData");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  const Double_t MASS_LOW  = 40;
  const Double_t MASS_HIGH = 200;
  const Double_t PT_CUT    = 22;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;
  const Int_t LEPTON_ID = 13;
  gStyle->SetOptStat(111111);
  // load trigger menu                                                                                                  
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  //Initialize histograms
  TH1F* TriMuFra = new TH1F("TriMuFra","Fraction of trigger object",51,0,1.02);
  TH1F* NumGlobal = new TH1F("NUM","Number of global muon",10,0,10);
  TH1F* hist0 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 0","m_{#mu^{+}#mu^{-}} (Data)",100,0,5);
  TH1F* hist1 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 1","pT",100,0,20);
  TH1F* hist2 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 2","#mu^{-} pT",100,0,20);
  TH1F* hist3 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 3","pT",100,0,20);
  TH1F* hist4 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 4","#eta",100,-5,5);
  TH1F* hist5 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 5","#mu^{-} #eta",100,-5,5);
  TH1F* hist6 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 6","#pi^{-} #eta",100,-5,5);
  TH1F* hist7 = new TH1F("#tau^{#mp} -> #mu^{#mp} #mu^{#pm} #mu^{#mp} 7","m_{#mu^{#mp}#mu^{#pm}#mu^{#mp}} (Data)",200,0,10);
  Int_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0; 

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  //
  // parse .conf file
  //
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  
  //
  // Declare output ntuple variables
  //TLorentzVector *v1=0, *v2=0, *v3=0;
  Float_t sysinvmass;
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  //baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  //TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex"); 
  TFile *infile=0;
  TTree *eventTree=0;

    
  //
  // loop over samples
  //  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    //>>>>>>>>>>>>>>>>>>>>>
    
    CSample* samp = samplev[isam];
    
    //
    // Set up output ntuple
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_Tau3Mu_select.root");
    //if(isam!=0 && !doScaleCorr) outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("sysinvmass", &sysinvmass,"sysinvmass/F");
    //outTree->Branch("v1","TLorentzVector", &v1);
    //outTree->Branch("v2","TLorentzVector", &v2);
    //outTree->Branch("v3","TLorentzVector", &v3);
    //>>>>>>>>>>>>>>>>>>>>>
    
    //
    // loop through files
    //
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  
      
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();
      infile = TFile::Open(samp->fnamev[ifile]); 
      assert(infile);
      if (samp->fnamev[ifile] == "/dev/null") 
	{
	  cout <<"-> Ignoring null input "<<endl; 
	  continue;
	}


      Bool_t hasJSON = kFALSE;
      baconhep::RunLumiRangeMap rlrm;
      if(samp->jsonv[ifile].CompareTo("NONE")!=0) { 
        hasJSON = kTRUE;
	rlrm.addJSONFile(samp->jsonv[ifile].Data()); 
      }
  
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
      eventTree->SetBranchAddress("Info", &info);      TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr);   TBranch *muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("Vertex",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("Vertex");

      //GEN handle
      //>>>>>>>>>>>>>>>>>>>>>>>>
   
      //
      // loop over events
      //
      Double_t nsel=0, nselvar=0;
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<10000; ientry++) {
        infoBr->GetEntry(ientry);
	count1++;

	if(ientry%5000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

	//weight staff
	//>>>>>>>>>>>>>>>>>>>>>
     
        // check for certified lumi (if applicable)
        baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.hasRunLumi(rl)) continue;
	count2++;

        // trigger requirement               
        if(!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
	count3++;

        // good vertex requirement
        if(!(info->hasGoodPV)) continue;
	count4++;
	
	//Select muon system
	muonArr->Clear();
        muonBr->GetEntry(ientry);
	baconhep::TMuon* mu[6] = {NULL,NULL,NULL,NULL,NULL,NULL};
	baconhep::TMuon *mutemp1 = NULL;
	baconhep::TMuon *mutemp2 = NULL;
	baconhep::TMuon *mutemp3 = NULL;
	vector<baconhep::TMuon*> triplet;
	vector<vector<baconhep::TMuon*> > tripletcollect;
	TLorentzVector vmu[6];
	TLorentzVector vtemp1, vtemp2, vtemp3;
	Double_t massdiff = 0.06;
	for(int i=0; i<muonArr->GetEntriesFast(); i++){
	  mutemp1 = (baconhep::TMuon*)(*muonArr)[i];
	  if(!(mutemp1->typeBits & baconhep::EMuType::kTracker)) continue;
	  vtemp1.SetPtEtaPhiM(mutemp1->pt, mutemp1->eta, mutemp1->phi, MUON_MASS);
	  for(int j=i+1; j<muonArr->GetEntriesFast(); j++){
	    mutemp2 = (baconhep::TMuon*)(*muonArr)[j];
	    if(!(mutemp2->typeBits & baconhep::EMuType::kTracker)) continue;
	    vtemp2.SetPtEtaPhiM(mutemp2->pt, mutemp2->eta, mutemp2->phi, MUON_MASS);
	    for(int k=j+1; k<muonArr->GetEntriesFast(); k++){
	      mutemp3 = (baconhep::TMuon*)(*muonArr)[k];
	      if(!(mutemp3->typeBits & baconhep::EMuType::kTracker)) continue;
	      vtemp3.SetPtEtaPhiM(mutemp3->pt, mutemp3->eta, mutemp3->phi, MUON_MASS);

	      Int_t NumMu = 0;
	      if(mutemp1->q == mutemp2->q && mutemp1->q == mutemp3->q) continue;
	      if(mutemp1->typeBits & baconhep::EMuType::kTracker) NumMu++;
	      if(mutemp2->typeBits & baconhep::EMuType::kTracker) NumMu++;
	      if(mutemp3->typeBits & baconhep::EMuType::kTracker) NumMu++;
	      if(NumMu < 3) continue; //enforce the selected combination to have at least 3 tracker muons,
	      Double_t invmass = (vtemp1+vtemp2+vtemp3).M();
	      if(fabs(invmass - 1.77682) < massdiff){
		triplet.push_back(mutemp1);
		triplet.push_back(mutemp2);
		triplet.push_back(mutemp3);
		tripletcollect.push_back(triplet);
	      }
	    }
	  }
	}
	cout<<"number of triplet candidates: "<<tripletcollect.size()<<" evnets: "<<ientry<<endl;
      }//end of event loop
      delete infile;
      infile=0, eventTree=0;    
    }
    outFile->Write();
    outFile->Close(); 
  }
  delete info;
  delete muonArr;
  delete vertexArr;
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> mu mu" << endl;
  cout << "  Mass window: [" << MASS_LOW << ", " << MASS_HIGH << "]" << endl;
  cout << "  pT > " << PT_CUT << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;
  cout<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<endl;
}

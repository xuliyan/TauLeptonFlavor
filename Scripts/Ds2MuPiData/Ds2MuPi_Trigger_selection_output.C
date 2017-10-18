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

void selectDsPhiPi(const TString conf="samples.conf", // input file
               const TString outputDir=".",   // output directory
	       const Bool_t  doScaleCorr=0,    // apply energy scale corrections
	       const Bool_t  doPU=0
	       ) {
  gBenchmark->Start("selectDsPhiPi");

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
  TH1F* NumGlobal = new TH1F("NUM","Number of global muon",10,0,10);
  TH1F* hist0 = new TH1F("Ds^{-} -> #mu^{+} #mu^{-} #pi^{-} 0","m_{#mu^{+}#mu^{-}} (Data)",150,0,5);
  TH1F* hist1 = new TH1F("Ds^{-} -> #mu^{+} #mu^{-} #pi^{-} 1","#mu pT",100,0,20);
  TH1F* hist2 = new TH1F("Ds^{-} -> #mu^{+} #mu^{-} #pi^{-} 2","#mu^{-} pT",100,0,20);
  TH1F* hist3 = new TH1F("Ds^{-} -> #mu^{+} #mu^{-} #pi^{-} 3","#pi^{-} pT",75,0,15);
  TH1F* hist4 = new TH1F("Ds^{-} -> #mu^{+} #mu^{-} #pi^{-} 4","#mu #eta",100,-5,5);
  TH1F* hist5 = new TH1F("Ds^{-} -> #mu^{+} #mu^{-} #pi^{-} 5","#mu^{-} #eta",100,-5,5);
  TH1F* hist6 = new TH1F("Ds^{-} -> #mu^{+} #mu^{-} #pi^{-} 6","#pi^{-} #eta",100,-5,5);
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
  /*
  UInt_t runNum, evtNum, lumiSec, metFilterFailBits;
  UInt_t nPU, nPUm, nPUp;
  Float_t nPUmean, nPUmeanm, nPUmeanp;
  Float_t pvx, pvy, pvz, bsx, bsy, bsz;
  Float_t pfMET, pfMETphi, pfMETCov00, pfMETCov01, pfMETCov11, pfMETC, pfMETCphi, pfMETCCov00, pfMETCCov01, pfMETCCov11;
  Float_t mvaMET, mvaMETphi, mvaMETCov00, mvaMETCov01, mvaMETCov11, mvaMETU, mvaMETUphi, mvaMETUCov00, mvaMETUCov01, mvaMETUCov11, mvaMET0, mvaMET0phi, mvaMET0Cov00, mvaMET0Cov01, mvaMET0Cov11;
  Float_t puppET, puppETphi, puppETCov00, puppETCov01, puppETCov11;
  Float_t trkMET, trkMETphi;
  Float_t rhoIso, rhoJet;
  UInt_t triggerBits;
  //muon specific
  Float_t pt, eta, phi;
  Float_t ptErr, staPt, staEta, staPhi;
  Float_t pfPt, pfEta, pfPhi;
  Float_t trkIso, ecalIso, hcalIso, chHadIso, gammaIso, neuHadIso, puIso;
  Float_t d0, dz, sip3d;
  Float_t tkNchi2, muNchi2, trkKink, glbKink;
  */
    
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
    //
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_trigger_selection_57.root");
    //if(isam!=0 && !doScaleCorr) outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("Info","baconhep::TEventInfo",&info);      // event run number
    outTree->Branch("Muon","TClonesArray",&muonArr);        // event lumi section
    outTree->Branch("Vertex","TClonesArray",&vertexArr);  // event number
    
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
      //or(UInt_t ientry=0; ientry<100000; ientry++) {
        infoBr->GetEntry(ientry);
	count1++;

	if(ientry%20000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

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
	
	//Loop over muon collection
	muonArr->Clear();
        muonBr->GetEntry(ientry);
	vertexArr->Clear();
	vertexBr->GetEntry(ientry);
	
	outTree->Fill();
	  
      }//end of event loop
      outTree->Print();
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

  //Draw
  TCanvas *c0 = new TCanvas("1","1",1200,900);
  TAxis *xaxis = hist0->GetXaxis();
  TAxis *yaxis = hist0->GetYaxis();
  xaxis->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
  yaxis->SetTitle("Entries / 0.1 GeV");
  yaxis->SetTitleOffset(1.2);
  //c0->SetLogx();
  c0->SetLogy();
  c0->cd();
  hist0->Draw();
  auto legend = new TLegend(0.50,0.76,0.65,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist0,"#mu^{+}#mu^{-}","f");
  legend->Draw();
  c0->Print("invmass.png");

  TCanvas *c1 = new TCanvas("2","2",1200,900);
  xaxis = hist1->GetXaxis();
  yaxis = hist1->GetYaxis();
  xaxis->SetTitle("pT (GeV)");
  yaxis->SetTitle("Entries / 0.2 GeV");
  yaxis->SetTitleOffset(1.2);
  c1->cd();
  c1->SetLogy();
  hist1->SetLineColor(2);
  hist2->SetLineColor(4);
  hist1->SetFillStyle(0);
  hist2->SetFillStyle(0);
  hist1->Draw();
  hist2->Draw("same");
  legend = new TLegend(0.50,0.72,0.6,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist1,"#mu^{+}","f");
  legend->AddEntry(hist2,"#mu^{-}","f");
  legend->Draw();
  c1->Print("pt.png");

  TCanvas *c2 = new TCanvas("3","3",1200,900);
  xaxis = hist4->GetXaxis();
  yaxis = hist4->GetYaxis();
  xaxis->SetTitle("#eta");
  yaxis->SetTitle("Entries / 0.1");
  yaxis->SetTitleOffset(1.2);
  c2->cd();
  c2->SetLogy();
  hist4->SetLineColor(2);
  hist5->SetLineColor(4);
  hist4->SetFillStyle(0);
  hist5->SetFillStyle(0);
  hist4->Draw();
  hist5->Draw("same");
  legend = new TLegend(0.45,0.72,0.55,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist4,"#mu^{+}","f");
  legend->AddEntry(hist5,"#mu^{-}","f");
  legend->Draw();
  c2->Print("eta.png");
  gBenchmark->Show("selectDsPhiPi"); 
}

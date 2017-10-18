//================================================================================================
//
// Select tau->mumumu candidates
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

void select3Mu(const TString conf="samples.conf", // input file
	       const TString outputDir=".",  // output directory
	       const Bool_t  doScaleCorr=0   // apply energy scale corrections?
	       ) {
  gBenchmark->Start("select3Mu");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  const Bool_t isGENLevel = kFALSE;
  const Double_t PT_CUT_LEAD    = 2.5;
  const Double_t PT_CUT_SUB = 1;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;
  
  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 13;

  TH1F* NUM = new TH1F("NUM","Number of global muon",10,0,10);
  TH1F* MC = new TH1F("MC truth","MC truth",100,0,10);
  TH1F* R1 = new TH1F("DeltaR","DeltaR",100,0,0.03);
  TH1F* R2 = new TH1F("DeltaR1","DeltaR",100,0,0.03);
  TH1F* R3 = new TH1F("DeltaR2","DeltaR",100,0,0.03);

  TH1F* hist0 = new TH1F("tau -> 3 muon 0","RECO invariant mass",100,0,10);
  TH1F* hist1 = new TH1F("tau -> 3 muon 1","invariant mass",100,0,10);
  TH1F* hist2 = new TH1F("tau -> 3 muon 2","RECO pT",75,0,15);
  TH1F* hist3 = new TH1F("tau -> 3 muon 3","RECO Eta",50,-3,3);
  TH1F* hist4 = new TH1F("tau -> 3 muon 4","RECO Phi",25,-3.5,3.5);
  TH1F* hist5 = new TH1F("tau -> 3 muon 5","RECO pT",75,0,15);
  TH1F* hist6 = new TH1F("tau -> 3 muon 6","RECO Eta",50,-3,3);
  TH1F* hist7 = new TH1F("tau -> 3 muon 7","RECO Phi",25,-3.5,3.5);
  TH1F* hist8 = new TH1F("tau -> 3 muon 8","RECO pT",75,0,15);
  TH1F* hist9 = new TH1F("tau -> 3 muon 9","RECO Eta",50,-3,3);
  TH1F* hist10 = new TH1F("tau -> 3 muon 10","RECO Phi",25,-3.5,3.5);
  TH1F* hist11 = new TH1F("tau -> 3 muon 11","pT",25,0,10);
  TH1F* hist12 = new TH1F("tau -> 3 muon 12","Eta",50,-3,3);
  TH1F* hist13 = new TH1F("tau -> 3 muon 13","Phi",25,-3.5,3.5);
  TH1F* hist14 = new TH1F("tau -> 3 muon 14","pT",25,0,10);
  TH1F* hist15 = new TH1F("tau -> 3 muon 15","Eta",50,-3,3);
  TH1F* hist16 = new TH1F("tau -> 3 muon 16","Phi",25,-3.5,3.5);
  TH1F* hist17 = new TH1F("tau -> 3 muon 17","pT",25,0,10);
  TH1F* hist18 = new TH1F("tau -> 3 muon 18","Eta",50,-3,3);
  TH1F* hist19 = new TH1F("tau -> 3 muon 19","Phi",25,-3.5,3.5);
  TH1F* hist20 = new TH1F("tau -> 3 muon 20","GEN deltaR Sum",100,0,10);
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0, lessmuon=0, lowpt = 0, multidecay = 0;
  gStyle->SetOptStat(111111);

  
  // load trigger menu
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Utils/data/puWeights_76x.root", "read");


  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
  TH1D *h_rw_up = (TH1D*) f_rw->Get("puWeightsUp");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("puWeightsDown");



  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  // parse .conf file
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  
  // Declare output ntuple variables
  UInt_t  runNum, lumiSec, evtNum;
  UInt_t  npv, npu;
  Float_t scale1fb, scale1fbUp, scale1fbDown, puWeight,puWeightUp,puWeightDown;
  Float_t met, metPhi, sumEt;
  Float_t tkMet, tkMetPhi, tkSumEt;
  //Float_t mvaMet, mvaMetPhi, mvaSumEt, mvaMt, mvaU1, mvaU2;
  //Float_t puppiMet, puppiMetPhi, puppiSumEt, puppiMt, puppiU1, puppiU2;
  Int_t   q[3];
  TLorentzVector *lep1=0, *lep2=0, *lep3=0;
  Int_t lepID[3];
  ///// muon specific /////
  Float_t trkIso[3], emIso[3], hadIso[3];
  Float_t pfChIso[3], pfGamIso[3], pfNeuIso[3], pfCombIso[3];
  Float_t d0[3], dz[3];
  Float_t muNchi2[3];
  UInt_t nPixHits[3], nTkLayers[3], nValidHits[3], nMatch[3], typeBits[3];
  UInt_t TolMuNum;
  UInt_t TolPassPreSel=0;
  
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *gen  = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *muonArr    = new TClonesArray("baconhep::TMuon");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  
  TFile *infile=0;
  TTree *eventTree=0;

  // loop over samples
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    
    // Assume data sample is first sample in .conf file
    // If sample is empty (i.e. contains no ntuple files), skip to next sample
    Bool_t isData=kFALSE;
    if(isam==0 && !hasData) continue;
    else if (isam==0) isData=kTRUE;

    // Assume signal sample is given name "dstau" -- flag to store GEN Ds kinematics
    Bool_t isSignal = (snamev[isam].CompareTo("dstau",TString::kIgnoreCase)==0);
    CSample* samp = samplev[isam];
    
    cout<<"begin loop over files"<<endl;
    // loop through files
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  
      
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();
      infile = TFile::Open(samp->fnamev[ifile]); 
      assert(infile);

      // Access Event Tree
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);  
      eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PV",   &vertexArr); TBranch *vertexBr = eventTree->GetBranch("PV");
      Bool_t hasGen = eventTree->GetBranchStatus("GenEvtInfo");
      TBranch *genBr=0, *genPartBr=0;
      eventTree->SetBranchAddress("GenEvtInfo", &gen); genBr = eventTree->GetBranch("GenEvtInfo");
      eventTree->SetBranchAddress("GenParticle",&genPartArr); genPartBr = eventTree->GetBranch("GenParticle");

      //Loop over events
      int passpresel = 0;

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<5000; ientry++) {
        infoBr->GetEntry(ientry);
	
        if(ientry%20000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
   
        if(!isGENLevel){
	  // check for certified lumi (if applicable)
	  baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
	  //if(hasJSON && !rlrm.hasRunLumi(rl)) continue;  
	  
	  // trigger requirement               
	  //if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
	  
	  // good vertex requirement
	  //if(!(info->hasGoodPV)) continue;
	  passpresel++;
	}
	
	//////////////////////MATCHING/////////////////////////////////	
	genPartArr->Clear();
	genPartBr->GetEntry(ientry);
	muonArr->Clear();
	muonBr->GetEntry(ientry);

	//Store GEN level tau muon
	vector<baconhep::TGenParticle*> genmuonArr;
	for(int i=0; i<genPartArr->GetEntries(); i++){
	  baconhep::TGenParticle *genpar = (baconhep::TGenParticle*)((*genPartArr)[i]);
	  if(genpar->pdgId != 13 && genpar->pdgId != -13) continue;
	  if(genpar->status != 1) continue;
	  Int_t parentid1=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genpar->parent>-1 ? genpar->parent : 0))->pdgId;
	  if(parentid1 != 15 && parentid1 != -15) continue;
	  genmuonArr.push_back(genpar);
	}

	//Only study events contain single decay 229041
	if(genmuonArr.size() > 3) continue;
	count1++;

	/*
	if(genmuonArr[0]->eta > 2.4 ||
	   genmuonArr[1]->eta > 2.4 ||
	   genmuonArr[2]->eta > 2.4 ||
	   genmuonArr[0]->pt < 2.5 ||
	   genmuonArr[1]->pt < 2.5 ||
	   genmuonArr[2]->pt < 2.5) continue;
	count2++;
	*/
	
	//Find the number of global muons
	Int_t NUMGLOBAL = 0;
	for(int j=0; j<muonArr->GetEntries(); j++){
	  if(((baconhep::TMuon*)(*muonArr)[j])->typeBits & baconhep::EMuType::kGlobal) NUMGLOBAL++;
	}
	if(NUMGLOBAL < 2) {
	  lessmuon++;
	  continue;
	}
	count2++;

	//Find events contain low pt muons
	Bool_t isLowpt = kFALSE;
	for(int i=0; i<genmuonArr.size(); i++){
	  if(genmuonArr[i]->pt < 1) isLowpt = kTRUE;
	}
	if(isLowpt) lowpt++;
	
	//Initialize deltaR 2D array
	vector<vector<double> > deltaR(genmuonArr.size());
	for(int i=0; i<genmuonArr.size(); i++){
	  deltaR[i].resize(muonArr->GetEntries());
	  for(int j=0; j<muonArr->GetEntries(); j++){
	    deltaR[i][j] = 9999; //initialize
	  }
	}

	//Calculate deltaR
	for(int i=0; i<genmuonArr.size(); i++){
	  for(int j=0; j<muonArr->GetEntries(); j++){
	    //Only calculate global muon	    
	    if(!(((baconhep::TMuon*)(*muonArr)[j])->typeBits & baconhep::EMuType::kGlobal || ((baconhep::TMuon*)(*muonArr)[j])->typeBits & baconhep::EMuType::kTracker)) continue; //global muon
	    if(genmuonArr[i]->pdgId * ((baconhep::TMuon*)(*muonArr)[j])->q > 0) continue;//same sign
	    Double_t deltaR1 = toolbox::deltaR(genmuonArr[i]->eta, genmuonArr[i]->phi, ((baconhep::TMuon*)(*muonArr)[j])->eta, ((baconhep::TMuon*)(*muonArr)[j])->phi);
	    deltaR[i][j] = deltaR1;
	  }
	}

	//Find the smallest deltaR
	Double_t R[6] = {0,0,0,0,0,0};
	Double_t min = 8888, submin = 8888, subsubmin = 8888;
	Int_t rr[3][2];
	for(int k=0; k<3; k++){
	  for(int j=0; j<2; j++){
	    rr[k][j] = -99;
	  }
	}
	for(int k=0; k<genmuonArr.size(); k++){
	  for(int j=0; j<muonArr->GetEntries(); j++){
	    if(deltaR[k][j] < min){
	      min = deltaR[k][j];
	      rr[0][0] = k;
	      rr[0][1] = j;
	    }
	  }
	}
	for(int k=0; k<genmuonArr.size(); k++){
	  if(k==rr[0][0]) continue;
	  for(int j=0; j<muonArr->GetEntries(); j++){
	    if(j==rr[0][1]) continue;
	    if(deltaR[k][j] < submin){
	      submin = deltaR[k][j];
	      rr[1][0] = k;
	      rr[1][1] = j;
	    }
	  }
	}
	for(int k=0; k<genmuonArr.size(); k++){
	  if(k==rr[0][0] || k==rr[1][0]) continue;
	  for(int j=0; j<muonArr->GetEntries(); j++){
	    if(j==rr[0][1] || j==rr[1][1]) continue;
	    if(deltaR[k][j] < subsubmin){
	      subsubmin = deltaR[k][j];
	      rr[2][0] = k;
	      rr[2][1] = j;
	    }
	  }
	}

	//Initialize container
	const baconhep::TMuon *SelMuonArr[6];
	const baconhep::TGenParticle *SelMuonArrG[6];

	//Check if we formed 3 RECO-GEN muon pair candidates
	Bool_t is3pair = kTRUE;
	for(int k=0; k<3; k++){
	  for(int j=0; j<2; j++){
	    if(rr[k][j] == -99) is3pair = kFALSE;
	  }
	}

	//Exclude events not matched
	if(!is3pair) {
	  /*
	  if(NUMGLOBAL < 3){
	    for(int i=0; i<genmuonArr.size(); i++){
	    cout<<"ID: "<<genmuonArr[i]->pdgId<<" status: "<<genmuonArr[i]->status<<" parent ID: "<<dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genmuonArr[i]->parent>-1 ? genmuonArr[i]->parent : 0))->pdgId<<" pt: "<<genmuonArr[i]->pt<<endl;
	    }
	    for(int i=0; i<genPartArr->GetEntries(); i++){
	      baconhep::TGenParticle *genpar = (baconhep::TGenParticle*)((*genPartArr)[i]);
	      if(genpar->pdgId != 13 && genpar->pdgId != -13) continue;
	      cout<<"ID: "<<genpar->pdgId<<" status: "<<genpar->status<<" parent ID: "<<dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genpar->parent>-1 ? genpar->parent : 0))->pdgId<<" pt: "<<genpar->pt<<" eta: "<<genpar->eta<<endl;
	    }
	    for(int i=0; i<muonArr->GetEntries(); i++){
	    const baconhep::TMuon *muon = (baconhep::TMuon*)((*muonArr)[i]);
	    if(!(muon->typeBits & baconhep::EMuType::kGlobal)) continue; //global muon
	    cout<<"type: "<<muon->typeBits<<" q: "<<muon->q<<" pt: "<<muon->pt<<" eta: "<<muon->eta<<endl;
	    }
	    cout<<"end of event"<<endl;
	    cout<<endl;
	  }
	    */
	  continue;
	}
	count3++;

	//Store results
	SelMuonArrG[0] = genmuonArr[rr[0][0]];
	SelMuonArrG[1] = genmuonArr[rr[1][0]];
	SelMuonArrG[2] = genmuonArr[rr[2][0]];
	SelMuonArr[0] = ((baconhep::TMuon*)(*muonArr)[rr[0][1]]);
	SelMuonArr[1] = ((baconhep::TMuon*)(*muonArr)[rr[1][1]]);
	SelMuonArr[2] = ((baconhep::TMuon*)(*muonArr)[rr[2][1]]);
	R[0] = min;
	R[1] = submin;
	R[2] = subsubmin;

	
	//Exclude non 2 global 1 trk events
	Bool_t is2global1trk = kTRUE;
	Int_t glocount = 0;
	Int_t trkcount = 0;
	for(int i=0; i<3; i++){
	  if(SelMuonArr[i]->typeBits & baconhep::EMuType::kGlobal) glocount++;
	  else if
	    (SelMuonArr[i]->typeBits & baconhep::EMuType::kTracker) trkcount++;
	}
	if(glocount < 2) is2global1trk = kFALSE;
	if(!is2global1trk) continue;
	count4++;
	

	//Sort muon array
	Double_t maxpt=0, submaxpt=0;
	for(int l=0; l<3; l++){
	  if(SelMuonArr[l]->pt >= maxpt){
	    submaxpt = maxpt;
	    maxpt = SelMuonArr[l]->pt;	  
	    SelMuonArr[5]=SelMuonArr[4];
	    SelMuonArr[4]=SelMuonArr[3];
	    SelMuonArr[3]=SelMuonArr[l];
	    SelMuonArrG[5]=SelMuonArrG[4];
	    SelMuonArrG[4]=SelMuonArrG[3];
	    SelMuonArrG[3]=SelMuonArrG[l];
	    R[5] = R[4];
	    R[4] = R[3];
	    R[3] = R[l];
	  }
	  else if(SelMuonArr[l]->pt < maxpt && SelMuonArr[l]->pt >= submaxpt){
	    submaxpt = SelMuonArr[l]->pt;
	    SelMuonArr[5]=SelMuonArr[4];
	    SelMuonArr[4]=SelMuonArr[l];
	    SelMuonArrG[5]=SelMuonArrG[4];
	    SelMuonArrG[4]=SelMuonArrG[l];
	    R[5] = R[4];
	    R[4] = R[l];
	  }
	  else{
	    SelMuonArr[5]=SelMuonArr[l];
	    SelMuonArrG[5]=SelMuonArrG[l];
	    R[5] = R[l];
	  }
	}

	//

	//Match requirement
	Int_t passMatch = 0;
	Double_t deltaR1 = toolbox::deltaR(SelMuonArrG[3]->eta,SelMuonArrG[3]->phi,SelMuonArr[3]->eta,SelMuonArr[3]->phi);
	Double_t deltaR2 = toolbox::deltaR(SelMuonArrG[4]->eta,SelMuonArrG[4]->phi,SelMuonArr[4]->eta,SelMuonArr[4]->phi);
	Double_t deltaR3 = toolbox::deltaR(SelMuonArrG[5]->eta,SelMuonArrG[5]->phi,SelMuonArr[5]->eta,SelMuonArr[5]->phi);
	if(deltaR1 < 0.014) passMatch++;
	if(deltaR2 < 0.017) passMatch++;
	if(deltaR3 < 0.025) passMatch++;
	//if(1){
	if(passMatch == 3){
	  count5++;
	  R1->Fill(deltaR1);
	  R2->Fill(deltaR2);
	  R3->Fill(deltaR3);
	  TLorentzVector vtempG[3],vtemp[3];
	  for(int p=0; p<3; p++){
	    vtemp[p].SetPtEtaPhiM(SelMuonArr[p+3]->pt, SelMuonArr[p+3]->eta, SelMuonArr[p+3]->phi, MUON_MASS);
	    vtempG[p].SetPtEtaPhiM(SelMuonArrG[p+3]->pt, SelMuonArrG[p+3]->eta, SelMuonArrG[p+3]->phi, MUON_MASS);
	  }
	  Double_t invmass = (vtemp[0]+vtemp[1]+vtemp[2]).M();
	  Double_t invmassG = (vtempG[0]+vtempG[1]+vtempG[2]).M();
	  //if(invmassG > 1.77681 && invmassG < 1.77683){
	  if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
	  hist0->Fill(invmass);//store events come from tau
	  //if(deltaR1 < 0.03) hist2->Fill(SelMuonArrG[3]->pt); if(deltaR2 < 0.03) hist5->Fill(SelMuonArrG[4]->pt); if(deltaR3 < 0.03) hist8->Fill(SelMuonArrG[5]->pt);
	  hist2->Fill(SelMuonArr[3]->pt); hist5->Fill(SelMuonArr[4]->pt); hist8->Fill(SelMuonArr[5]->pt);
	  hist3->Fill(SelMuonArr[3]->eta); hist6->Fill(SelMuonArr[4]->eta); hist9->Fill(SelMuonArr[5]->eta);
	  hist4->Fill(SelMuonArr[3]->phi); hist7->Fill(SelMuonArr[4]->phi); hist10->Fill(SelMuonArr[5]->phi);	      
	  count6++;
	  /*
	    cout<<"GEN 1: pt: "<<SelMuonArrG[3]->pt<<" eta: "<<SelMuonArrG[3]->eta<<" phi: "<<SelMuonArrG[3]->phi<<" ID: "<<SelMuonArrG[3]->pdgId<<" parent: "<<SelMuonArrG[3]->parent<<" parent ID: "<<((baconhep::TGenParticle*)((*genPartArr)[(SelMuonArrG[3]->parent>-1 ? SelMuonArrG[3]->parent : 0)]))->pdgId<<endl;

	    cout<<"SIM 1: pt: "<<SelMuonArr[3]->pt<<" eta: "<<SelMuonArr[3]->eta<<" phi: "<<SelMuonArr[3]->phi<<" deltaR: "<<R[3]<<endl;

	    cout<<"GEN 2: pt: "<<SelMuonArrG[4]->pt<<" eta: "<<SelMuonArrG[4]->eta<<" phi: "<<SelMuonArrG[4]->phi<<" ID: "<<SelMuonArrG[4]->pdgId<<" parent: "<<SelMuonArrG[4]->parent<<" parent ID: "<<((baconhep::TGenParticle*)((*genPartArr)[(SelMuonArrG[4]->parent>-1 ? SelMuonArrG[4]->parent : 0)]))->pdgId<<endl;

	    cout<<"SIM 2: pt: "<<SelMuonArr[4]->pt<<" eta: "<<SelMuonArr[4]->eta<<" phi: "<<SelMuonArr[4]->phi<<" deltaR: "<<R[4]<<endl;

	    cout<<"GEN 3: pt: "<<SelMuonArrG[5]->pt<<" eta: "<<SelMuonArrG[5]->eta<<" phi: "<<SelMuonArrG[5]->phi<<" ID: "<<SelMuonArrG[5]->pdgId<<" parent: "<<SelMuonArrG[5]->parent<<" parent ID: "<<((baconhep::TGenParticle*)((*genPartArr)[(SelMuonArrG[5]->parent>-1 ? SelMuonArrG[5]->parent : 0)]))->pdgId<<endl;

	    cout<<"SIM 3: pt: "<<SelMuonArr[5]->pt<<" eta: "<<SelMuonArr[5]->eta<<" phi: "<<SelMuonArr[5]->phi<<" deltaR: "<<R[5]<<endl;

	    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	    cout<<endl;
	  */
	}
	
      }//end of event loop
      cout<<"There is "<<passpresel<<" events passed pre selection"<<endl;
      TolPassPreSel+=passpresel;
      delete infile;
      infile=0, eventTree=0;    
    }
  }
  delete h_rw;
  delete h_rw_up;
  delete h_rw_down;
  delete f_rw;
  delete info;
  delete gen;
  delete genPartArr;
  delete muonArr;
  delete vertexArr;
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " tau -> mu mu mu" << endl;
  cout << "  pT > " << PT_CUT_LEAD << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;

  //Draw Graph
  //THStack *all = new THStack("tau -> 3 muons","invariant mass");
  TCanvas *c0 = new TCanvas("c0","invariant mass",1200,900);
  TAxis *xaxis = hist0->GetXaxis();
  TAxis *yaxis = hist0->GetYaxis();

  xaxis->SetTitle("Invariant mass (GeV)");
  yaxis->SetTitle("Entries / 0.1 GeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.5,100000);
  c0->cd();
  

  hist0->SetFillColor(5);
  hist0->SetFillStyle(0);
  //hist1->SetFillColor(7);
  //hist3->SetFillColor(7);
  //hist4->SetFillColor(8);


  hist0->Draw();
  //hist1->Draw("SAME");
  
  // all->Add(hist1);
  //all->Add(hist2);
  //all->Add(hist3);
  //all->Add(hist4);
  //c0->cd();
  c0->SetLogy();
  /*
    all->Draw();
    TAxis *xaxis = all->GetXaxis();
    TAxis *yaxis = all->GetYaxis();
    xaxis->SetTitle("Invariant mass (GeV)");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.2);
    all->SetMinimum(8.);
    all->SetMaximum(120000.);
  */
  TF1 *f1 = new TF1("m1","gaus",1.5,2);
  TF1 *f2 = new TF1("m2","gaus",0,5);
  TF1 *total = new TF1("mstotal","gaus(0)+gaus(3)",0,5);
  Double_t par[6]={50,1.5,0.1,5,2.5,2};
  hist0->Fit(f1,"R0");
  //hist2->Fit(f2,"R0+");
  //f1->GetParameters(&par[0]);
  //f2->GetParameters(&par[3]);
  total->SetParameters(par);
  hist0->Fit(total,"R0+");
  total->SetLineColor(2);
  total->SetLineWidth(2);
  //total->Draw("SAME");
  //f1->SetLineColor(3);
  //f1->Draw("SAME");

  auto legend = new TLegend(0.5,0.7,0.7,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist0,"2 global + 1 tracker muons","f");
  legend->Draw();
  

  c0->Print("invariant mass.png");
  cout<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<count6<<" "<<lessmuon<<" "<<multidecay<<" "<<lowpt<<endl;

  TCanvas *c1 = new TCanvas("1","muon pT",1200,900);
  xaxis = hist2->GetXaxis();
  yaxis = hist2->GetYaxis();

  xaxis->SetTitle("pT (GeV)");
  yaxis->SetTitle("Entries / 0.2 GeV");
  yaxis->SetTitleOffset(1.5);
  yaxis->SetRangeUser(0.5,100000);
  c1->cd();
  c1->SetLogy();

  hist2->SetLineColor(2);
  hist5->SetLineColor(6);
  hist8->SetLineColor(4);
  hist2->SetFillStyle(0);
  hist5->SetFillStyle(0);
  hist8->SetFillStyle(0);
  //hist11->SetLineColor(2);
  //hist14->SetLineColor(6);
  //hist17->SetLineColor(4);
  //hist11->SetLineStyle(2);
  //hist14->SetLineStyle(2);
  //hist17->SetLineStyle(2);

  hist2->Draw();
  hist5->Draw("SAME");
  hist8->Draw("SAME");
  //hist11->Draw("SAME");
  //hist14->Draw("SAME");
  //hist17->Draw("SAME");

  legend = new TLegend(0.35,0.75,0.6,0.85);
  legend->AddEntry(hist2,"2 global + 1 tracker muon, lead","f");
  legend->AddEntry(hist5,"2 global + 1 tracker muon, sublead","f");
  legend->AddEntry(hist8,"2 global + 1 tracker muon, third","f");
  //legend->AddEntry(hist11,"2 muons 1 track with same signs, lead","f");
  //legend->AddEntry(hist14,"2 muons 1 track with same signs, sublead","f");
  //legend->AddEntry(hist17,"2 muons 1 track with same signs, track","f");
  legend->Draw();

  c1->Print("pt.png");

  TCanvas *c2 = new TCanvas("2","muon eta",1200,900);
  xaxis = hist3->GetXaxis();
  yaxis = hist3->GetYaxis();

  xaxis->SetTitle("eta");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0,100);
  c2->cd();

  hist3->SetLineColor(2);
  hist6->SetLineColor(6);
  hist9->SetLineColor(4);
  hist3->SetFillStyle(0);
  hist6->SetFillStyle(0);
  hist9->SetFillStyle(0);
  //hist12->SetLineColor(2);
  //hist15->SetLineColor(6);
  //hist18->SetLineColor(4);
  //hist12->SetLineStyle(2);
  //hist15->SetLineStyle(2);
  //hist18->SetLineStyle(2);

  hist3->Draw();
  hist6->Draw("SAME");
  hist9->Draw("SAME");
  //hist12->Draw("SAME");
  //hist15->Draw("SAME");
  //hist18->Draw("SAME");

  legend = new TLegend(0.15,0.75,0.4,0.85);
  legend->AddEntry(hist3,"2 global + 1 tracker muon, lead","f");
  legend->AddEntry(hist6,"2 global + 1 tracker muon, sublead","f");
  legend->AddEntry(hist9,"2 global + 1 tracker muon, third","f");
  //legend->AddEntry(hist12,"2 muons 1 track with same signs, lead","f");
  //legend->AddEntry(hist15,"2 muons 1 track with same signs, sublead","f");
  //legend->AddEntry(hist18,"2 muons 1 track with same signs, track","f");
  legend->Draw();

  c2->Print("eta.png");

  TCanvas *c3 = new TCanvas("3","muon phi",1200,900);
  xaxis = hist4->GetXaxis();
  yaxis = hist4->GetYaxis();

  xaxis->SetTitle("phi");
  yaxis->SetTitle("Entries / 0.28 rad");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0,250);
  c3->cd();

  hist4->SetLineColor(2);
  hist7->SetLineColor(6);
  hist10->SetLineColor(4);
  //hist13->SetLineColor(2);
  //hist16->SetLineColor(6);
  //hist19->SetLineColor(4);
  //hist13->SetLineStyle(2);
  //hist16->SetLineStyle(2);
  //hist19->SetLineStyle(2);

  hist4->Draw();
  hist7->Draw("SAME");
  hist10->Draw("SAME");
  //hist13->Draw("SAME");
  //hist16->Draw("SAME");
  //hist19->Draw("SAME");

  legend = new TLegend(0.15,0.75,0.4,0.85);
  legend->AddEntry(hist4,"2 global + 1 tracker muon, lead","f");
  legend->AddEntry(hist7,"2 global + 1 tracker muon, sublead","f");
  legend->AddEntry(hist10,"2 global + 1 tracker muon, third","f");
  //legend->AddEntry(hist13,"2 muons 1 track with same signs, lead","f");
  //legend->AddEntry(hist16,"2 muons 1 track with same signs, sublead","f");
  //legend->AddEntry(hist19,"2 muons 1 track with same signs, track","f");
  legend->Draw();

  c3->Print("phi.png");

  TCanvas *c4 = new TCanvas("4","deltaR",1200,900);
  xaxis = R1->GetXaxis();
  yaxis = R1->GetYaxis();

  xaxis->SetTitle("deltaR");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.2);
  //yaxis->SetRangeUser(0,1000);
  c4->cd();

  R1->SetLineColor(2);
  R2->SetLineColor(6);
  R3->SetLineColor(4);
  //hist13->SetLineColor(2);
  //hist16->SetLineColor(6);
  //hist19->SetLineColor(4);
  //hist13->SetLineStyle(2);
  //hist16->SetLineStyle(2);
  //hist19->SetLineStyle(2);

  R1->Draw();
  R2->Draw("SAME");
  R3->Draw("SAME");
  //hist13->Draw("SAME");
  //hist16->Draw("SAME");
  //hist19->Draw("SAME");

  legend = new TLegend(0.35,0.75,0.6,0.85);
  legend->AddEntry(R1,"2 global + 1 tracker muon, lead","f");
  legend->AddEntry(R2,"2 global + 1 tracker muon, sublead","f");
  legend->AddEntry(R3,"2 global + 1 tracker muon, third","f");
  //legend->AddEntry(hist13,"2 muons 1 track with same signs, lead","f");
  //legend->AddEntry(hist16,"2 muons 1 track with same signs, sublead","f");
  //legend->AddEntry(hist19,"2 muons 1 track with same signs, track","f");
  legend->Draw();
  c4->SetLogy();

  c4->Print("deltaR.png");

  TCanvas *c5 = new TCanvas("5","Numglobalmuon",1200,900);
  xaxis = NUM->GetXaxis();
  yaxis = NUM->GetYaxis();

  xaxis->SetTitle("Number of global muons");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.2);
  //yaxis->SetRangeUser(0,1000);
  c5->cd();

  NUM->SetLineColor(4);
  //hist13->SetLineColor(2);
  //hist16->SetLineColor(6);
  //hist19->SetLineColor(4);
  //hist13->SetLineStyle(2);
  //hist16->SetLineStyle(2);
  //hist19->SetLineStyle(2);

  NUM->Draw();
  //hist13->Draw("SAME");
  //hist16->Draw("SAME");
  //hist19->Draw("SAME");

  c5->Print("num.png");

  gBenchmark->Show("select3Mu");

}

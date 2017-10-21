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

bool Findtaudecay(TClonesArray *muonArr,vector<TLorentzVector> arr,Int_t* record,const int arr_len, TH1F* NUM)
{
  Int_t NUMGLOBAL = 0;
  for(int p=0; p<arr_len; p++){
    if(((baconhep::TMuon*)(*muonArr)[p])->typeBits & baconhep::EMuType::kGlobal) NUMGLOBAL++;
  }
  NUM->Fill(NUMGLOBAL);
  Int_t multicount = 0;
  Bool_t isTaudecay = kFALSE;
  Double_t diff = 65535;
  for(int i=0; i<arr_len; i++){
    for(int j=i+1; j<arr_len; j++){
      for(int k=j+1; k<arr_len; k++){
	//code here:
	Int_t Nglobalmuon = 0;
	if(((baconhep::TMuon*)(*muonArr)[i])->q == ((baconhep::TMuon*)(*muonArr)[j])->q && ((baconhep::TMuon*)(*muonArr)[i])->q ==((baconhep::TMuon*)(*muonArr)[k])->q) continue;
	//if(fabs(((baconhep::TMuon*)(*muonArr)[i])->eta) > 1.6) continue;
	//if(fabs(((baconhep::TMuon*)(*muonArr)[j])->eta) > 1.6) continue;
	//if(fabs(((baconhep::TMuon*)(*muonArr)[k])->eta) > 1.6) continue;
	//(((baconhep::TMuon*)(*muonArr)[i])->pt < 1) continue;
	//if(((baconhep::TMuon*)(*muonArr)[j])->pt < 1) continue;
	//if(((baconhep::TMuon*)(*muonArr)[k])->pt < 1) continue;
	if(((baconhep::TMuon*)(*muonArr)[i])->typeBits & baconhep::EMuType::kGlobal) Nglobalmuon++;
	if(((baconhep::TMuon*)(*muonArr)[j])->typeBits & baconhep::EMuType::kGlobal) Nglobalmuon++;
	if(((baconhep::TMuon*)(*muonArr)[k])->typeBits & baconhep::EMuType::kGlobal) Nglobalmuon++;
	if(Nglobalmuon < 3) continue; //enforce the selected combination to have at least 2 global muons, 
	TLorentzVector vcom,v1,v2,v3;
	v1=arr[i];
	v2=arr[j];
	v3=arr[k];
	vcom =v1+v2+v3;
	Double_t invmass = vcom.M();
	//if(invmass > 1.6768 && invmass < 1.8768){
	isTaudecay = kTRUE;
	if((abs(invmass-1.77682))<diff){
	  diff = abs(invmass-1.77682);
	  record[0]=i;
	  record[1]=j;
	  record[2]=k;
	}
	multicount++;
	//}
      }
    }
  }
  return isTaudecay;
}

bool Findtaudecay(vector<baconhep::TGenParticle *> muonarr,vector<TLorentzVector> arr,Int_t* record,const int arr_len, baconhep::TGenParticle *genPartArr)
{
  Int_t multicount = 0;
  Bool_t isTaudecay = kFALSE;
  for(int i=0; i<arr_len; i++){
    for(int j=i+1; j<arr_len; j++){
      for(int k=j+1; k<arr_len; k++){
	//code here:
	if(!(muonarr[i]->parent == muonarr[j]->parent && muonarr[i]->parent == muonarr[k]->parent)) continue;
	if(muonarr[0]->pdgId == muonarr[1]->pdgId && muonarr[0]->pdgId == muonarr[2]->pdgId) continue;
	TLorentzVector vcom,v1,v2,v3;
	isTaudecay = kTRUE;
	record[0]=i;
	record[1]=j;
	record[2]=k;
	multicount++;
      }	
    }
  }
  if(multicount > 1) cout<<"more than one system have come from same tau"<<endl;
  return isTaudecay;
}

void selectGEN(const TString conf="samples.conf", // input file
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
  TH1F* R1 = new TH1F("DeltaR","DeltaR",100,0,1);
  TH1F* R2 = new TH1F("DeltaR1","DeltaR",100,0,1);
  TH1F* R3 = new TH1F("DeltaR2","DeltaR",100,0,1);

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
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, lessmuon=0, lowpt = 0, multidecay = 0;
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

	//Sort muon array
	Double_t maxpt=0, submaxpt=0, subsubpt=0;
	for(int l=0; l<3; l++){
	  if(genmuonArr[l]->pt >= maxpt){
	    subsubpt = submaxpt;
	    submaxpt = maxpt;
	    maxpt = genmuonArr[l]->pt;	  
	  }
	  else if(genmuonArr[l]->pt < maxpt && genmuonArr[l]->pt >= submaxpt){
	    subsubpt = submaxpt;
	    submaxpt = genmuonArr[l]->pt;
	  }
	  else{
	    subsubpt = genmuonArr[l]->pt;
	  }
	}

	if(ientry <20) cout<<maxpt<<" "<<submaxpt<<" "<<subsubpt<<endl;

	Bool_t pass = kTRUE;
	if(genmuonArr[0]->eta > 2.4 ||
	   genmuonArr[1]->eta > 2.4 ||
	   genmuonArr[2]->eta > 2.4 ||
	   maxpt < 2.5 ||
	   submaxpt < 2.5 ||
	   subsubpt < 1) pass  =kFALSE;
	if(pass)
	  count3++;
	   
	
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
  yaxis->SetRangeUser(0.5,40000);
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
  legend->AddEntry(hist0,"3 global muons with opposite signs","f");
  legend->Draw();
  

  c0->Print("invariant mass.png");
  cout<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<lessmuon<<" "<<multidecay<<" "<<lowpt<<endl;

  TCanvas *c1 = new TCanvas("1","muon pT",1200,900);
  xaxis = hist2->GetXaxis();
  yaxis = hist2->GetYaxis();

  xaxis->SetTitle("pT (GeV)");
  yaxis->SetTitle("Entries / 0.2 GeV");
  yaxis->SetTitleOffset(1.5);
  yaxis->SetRangeUser(0.5,45000);
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
  legend->AddEntry(hist2,"3 global muons with opposite signs, lead","f");
  legend->AddEntry(hist5,"3 global muons with opposite signs, sublead","f");
  legend->AddEntry(hist8,"3 global muons with opposite signs, third","f");
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
  yaxis->SetRangeUser(0,250);
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
  legend->AddEntry(hist3,"3 global muons with opposite signs, lead","f");
  legend->AddEntry(hist6,"3 global muons with opposite signs, sublead","f");
  legend->AddEntry(hist9,"3 global muons with opposite signs, third","f");
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
  legend->AddEntry(hist4,"3 global muons with opposite signs, lead","f");
  legend->AddEntry(hist7,"3 global muons with opposite signs, sublead","f");
  legend->AddEntry(hist10,"3 global muons with opposite signs, third","f");
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

  legend = new TLegend(0.15,0.75,0.4,0.85);
  legend->AddEntry(R1,"3 global muons with opposite signs, lead","f");
  legend->AddEntry(R2,"3 global muons with opposite signs, sublead","f");
  legend->AddEntry(R3,"3 global muons with opposite signs, third","f");
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

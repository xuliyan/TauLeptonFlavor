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

bool Findtaudecay(TClonesArray *muonArr,vector<TLorentzVector> arr,Int_t* record,const int arr_len) //Find muon system has invariant mass closest to tau
{
  Int_t multicount = 0;
  Bool_t isTaudecay = kFALSE;
  Double_t diff = 9999;
  for(int i=0; i<arr_len; i++){
    for(int j=i+1; j<arr_len; j++){
      for(int k=j+1; k<arr_len; k++){
	//code here:
	Int_t Nglobalmuon = 0;
	if(((baconhep::TMuon*)(*muonArr)[i])->q == ((baconhep::TMuon*)(*muonArr)[j])->q && ((baconhep::TMuon*)(*muonArr)[i])->q ==((baconhep::TMuon*)(*muonArr)[k])->q) continue;
	//if(fabs(((baconhep::TMuon*)(*muonArr)[i])->eta) > 1.6) continue;
	//if(fabs(((baconhep::TMuon*)(*muonArr)[j])->eta) > 1.6) continue;
	//if(fabs(((baconhep::TMuon*)(*muonArr)[k])->eta) > 1.6) continue;
	//if(((baconhep::TMuon*)(*muonArr)[i])->pt < 1) continue;
	//if(((baconhep::TMuon*)(*muonArr)[j])->pt < 1) continue;
	//if(((baconhep::TMuon*)(*muonArr)[k])->pt < 1) continue;
	if(((baconhep::TMuon*)(*muonArr)[i])->typeBits & baconhep::EMuType::kGlobal) Nglobalmuon++;
	if(((baconhep::TMuon*)(*muonArr)[j])->typeBits & baconhep::EMuType::kGlobal) Nglobalmuon++;
	if(((baconhep::TMuon*)(*muonArr)[k])->typeBits & baconhep::EMuType::kGlobal) Nglobalmuon++;
	if(Nglobalmuon < 2) continue; //enforce the selected combination to have at least 2 global muons, 
	TLorentzVector vcom,v1,v2,v3;
	v1=arr[i];
	v2=arr[j];
	v3=arr[k];
	vcom =v1+v2+v3;
	Double_t invmass = vcom.M();
	if(invmass > 1.67 && invmass < 1.87){
	  isTaudecay = kTRUE;
	  if((abs(invmass-1.77682))<diff){
	    diff = abs(invmass-1.77682);
	    record[0]=i;
	    record[1]=j;
	    record[2]=k;
	  }
	}	
      }
    }
  }
  return isTaudecay;
}

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
  
  const Int_t LEPTON_ID = 13;

  TH1F* parent_h = new TH1F("muon parent","parent id",505,-5,500);
  TH1F* hista = new TH1F("tau -> 3 muon 3B","invariant mass",25,1.5,2);
  TH1F* histb = new TH1F("tau -> 3 muon 2B1M","invariant mass",25,1.5,2);
  TH1F* histc = new TH1F("tau -> 3 muon 2B1E","invariant mass",25,1.5,2);
  TH1F* histd = new TH1F("tau -> 3 muon 3M","invariant mass",25,1.5,2);
  TH1F* histe = new TH1F("tau -> 3 muon 2M1B","invariant mass",25,1.5,2);
  TH1F* histf = new TH1F("tau -> 3 muon 2M1E","invariant mass",25,1.5,2);
  TH1F* histg = new TH1F("tau -> 3 muon 3E","invariant mass",25,1.5,2);
  TH1F* histh = new TH1F("tau -> 3 muon 2E1B","invariant mass",25,1.5,2);
  TH1F* histi = new TH1F("tau -> 3 muon 2E1M","invariant mass",25,1.5,2);
  TH1F* histj = new TH1F("tau -> 3 muon 1B1M1E","invariant mass",25,1.5,2);
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0, count7=0, count8=0, count9=0, count10=0, count11=0, count12=0;

  gStyle->SetOptStat(111111);

  
  // load trigger menu
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  // parse .conf file
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);
  
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
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<10000; ientry++) {
        infoBr->GetEntry(ientry);
	
        if(ientry%50000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
   
        if(!isGENLevel){
	// check for certified lumi (if applicable)
        //baconhep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        //if(hasJSON && !rlrm.hasRunLumi(rl)) continue;  

        // trigger requirement               
        //if (!isMuonTrigger(triggerMenu, info->triggerBits)) continue;
      
        // good vertex requirement
        //if(!(info->hasGoodPV)) continue;
	}
	
	// SELECTION PROCEDURE:
	Int_t globalmuon=0;
	Int_t trk=0;
	Int_t Numtrk[200];
	Int_t runcount=0;
	Int_t Numglomuon[50];
	Bool_t passSel=kFALSE;	
	//=================================BEGIN OF SIM LEVEL=====================================
	if(!isGENLevel){
	  muonArr->Clear();
	  muonBr->GetEntry(ientry);
	  // Looking for two gloabl muons
	  for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
	    const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[i]);
	    if(mu->typeBits & baconhep::EMuType::kGlobal){
	    globalmuon++;
	    Numglomuon[runcount] = i;
	    runcount++;
	    }
	  }//End of muon array loop
	}
	if (globalmuon > 1){
	  count1++;
	  passSel = kTRUE;
	}	
	if(1){//set kTRUE or pssSel here will not affect the result since Findtaudecay enforced selecred combination to have more than 2 global muons.
	  if (!isGENLevel){
	    const baconhep::TMuon *SelMuonArr[6];
	    Int_t record[3]={0,0,0};
	    vector<TLorentzVector> vLep;
	    TLorentzVector vLepRe[3], vCom;

	    //store global muon into array
	    for(int p=0; p<muonArr->GetEntriesFast(); p++){
	      TLorentzVector vtemp;
	      vtemp.SetPtEtaPhiM(((baconhep::TMuon*)(*muonArr)[p])->pt, ((baconhep::TMuon*)(*muonArr)[p])->eta, ((baconhep::TMuon*)(*muonArr)[p])->phi, MUON_MASS);
	      vLep.push_back(vtemp);//store corresponding 4-vector
	    }

	    if(Findtaudecay(muonArr,vLep,record,muonArr->GetEntriesFast())){    
	      for(int p=0; p<3; p++){
	      SelMuonArr[p] = (baconhep::TMuon*)((*muonArr)[record[p]]);
	      }	    
	      Double_t maxpt=0, submaxpt=0;
	      for(int l=0; l<3; l++){
		if(SelMuonArr[l]->pt >= maxpt){
		  submaxpt = maxpt;
		  maxpt = SelMuonArr[l]->pt;	  
		  SelMuonArr[5]=SelMuonArr[4];
		  SelMuonArr[4]=SelMuonArr[3];
		  SelMuonArr[3]=SelMuonArr[l];
		}
		else if(SelMuonArr[l]->pt < maxpt && SelMuonArr[l]->pt >= submaxpt){
		  submaxpt = SelMuonArr[l]->pt;
		  SelMuonArr[5]=SelMuonArr[4];
		  SelMuonArr[4]=SelMuonArr[l];
		}
		else{
		  SelMuonArr[5]=SelMuonArr[l];
		}
	      }
	      for(int p=3; p<6; p++){
		vLepRe[p-3].SetPtEtaPhiM(SelMuonArr[p]->pt, SelMuonArr[p]->eta, SelMuonArr[p]->phi, MUON_MASS);
	      }
	      vCom = vLepRe[0]+vLepRe[1]+vLepRe[2];
	      Double_t invmass = vCom.M();
	      count2++;
	      
	      Int_t B = 0;
	      Int_t M = 0;
	      Int_t E = 0;
	      if (fabs(SelMuonArr[3]->eta)<0.8) B++;
	      if (fabs(SelMuonArr[4]->eta)<0.8) B++;
	      if (fabs(SelMuonArr[5]->eta)<0.8) B++;
	      if (fabs(SelMuonArr[3]->eta)>=0.8 && fabs(SelMuonArr[3]->eta)<1.6) M++;
	      if (fabs(SelMuonArr[4]->eta)>=0.8 && fabs(SelMuonArr[4]->eta)<1.6) M++;
	      if (fabs(SelMuonArr[5]->eta)>=0.8 && fabs(SelMuonArr[5]->eta)<1.6) M++;
	      if (fabs(SelMuonArr[3]->eta)>=1.6) E++;
	      if (fabs(SelMuonArr[4]->eta)>=1.6) E++;
	      if (fabs(SelMuonArr[5]->eta)>=1.6) E++;
	      if (B==3){
		count3++;
		hista->Fill(invmass);}//store events come from tau
	      else if(B==2&&M==1){
		count3++;
		hista->Fill(invmass);}
	      else if(B==2&&E==1){
		count3++;
		hista->Fill(invmass);}
	      else if(M==3){
		count6++;
		histd->Fill(invmass);}
	      else if(M==2&&B==1){
		count6++;
		histd->Fill(invmass);}
	      else if(M==2&&E==1){
		count6++;
		histd->Fill(invmass);}
	      else if(E==3){
		count9++;
		histg->Fill(invmass);}
	      else if(E==2&&B==1){
		count9++;
		histg->Fill(invmass);}
	      else if(E==2&&M==1){
		count9++;
		histg->Fill(invmass);}
	      else if(B==1&&M==1&&E==1){
		count12++;
		histj->Fill(invmass);}	      
	      else{
		cout<<B<<M<<E<<endl;
	      }
	    }
	  }	  
	}
      }//end of event loop
      delete infile;
      infile=0, eventTree=0;    
    }
  }
  
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
  cout << " Lead pT > " << PT_CUT_LEAD << endl;
  cout << " Third pT > " << PT_CUT_SUB << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  //Draw Graph
  gStyle->SetOptStat(0);
  TCanvas *c0 = new TCanvas("c0","invariant mass",1200,900);
  TAxis *xaxis = histg->GetXaxis();
  TAxis *yaxis = histg->GetYaxis();

  xaxis->SetTitle("Invariant mass (GeV)");
  yaxis->SetTitle("Entries / 20 Mev");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.5,60000);
  c0->cd();
  c0->SetLogy();

  histg->SetLineColor(46);
  histj->SetLineColor(28);
  histi->SetLineColor(9);
  histh->SetLineColor(8);
  histf->SetLineColor(7);
  histe->SetLineColor(6);
  histd->SetLineColor(4);
  histc->SetLineColor(3);
  histb->SetLineColor(2);
  hista->SetLineColor(1);
  histg->SetFillStyle(0);
  hista->SetFillStyle(0);
  histb->SetFillStyle(0);
  histc->SetFillStyle(0);
  histd->SetFillStyle(0);
  histe->SetFillStyle(0);
  histf->SetFillStyle(0);
  histh->SetFillStyle(0);
  histi->SetFillStyle(0);
  histj->SetFillStyle(0);

  histg->Draw();
  hista->Draw("same");
  histb->Draw("same");
  histc->Draw("same");
  histd->Draw("same");
  histe->Draw("same");
  histf->Draw("same");
  histh->Draw("same");
  histi->Draw("same");
  histj->Draw("same");


  //histd->Draw();
  //gPad->Update();
  /* collect stat of the first histogram (h1) */
  /*
  TPaveStats *tps1 = (TPaveStats*) histd->FindObject("stats");
  tps1->SetName("Hist1 Stats");
  tps1->SetTextColor(4);
  tps1->SetLineColor(4);
  double X1 = tps1->GetX1NDC();
  double X2 = tps1->GetX2NDC();

  double Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8;
  Y1=0.735;
  Y2=0.935;
  Y3=0.535;
  Y4=Y1;
  Y5=0.335;
  Y6=Y3;
  Y7=0.135;
  Y8=Y5;
  

  histc->Draw();
  gPad->Update();
  TPaveStats *tps2 = (TPaveStats*) histc->FindObject("stats");
  tps2->SetTextColor(6);
  tps2->SetLineColor(6);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y3);
  tps2->SetY2NDC(Y4);
  
  histb->Draw();
  gPad->Update();
  TPaveStats *tps3 = (TPaveStats*) histb->FindObject("stats");
  tps3->SetTextColor(8);
  tps3->SetLineColor(8);
  tps3->SetX1NDC(X1);
  tps3->SetX2NDC(X2);
  tps3->SetY1NDC(Y5);
  tps3->SetY2NDC(Y6);
  
  hista->Draw();
  gPad->Update();
  TPaveStats *tps4 = (TPaveStats*) hista->FindObject("stats");
  tps4->SetTextColor(2);
  tps4->SetLineColor(2);
  tps4->SetX1NDC(X1);
  tps4->SetX2NDC(X2);
  tps4->SetY1NDC(Y7);
  tps4->SetY2NDC(Y8);
  

  histd->Draw();
  histc->Draw("same");
  histb->Draw("same");
  hista->Draw("same");
  tps1->Draw("same");
  tps2->Draw("same");
  tps3->Draw("same");
  tps4->Draw("same");
  */
  
  auto legend = new TLegend(0.15,0.78,0.35,0.88);
  legend->AddEntry(hista,"3 global muons 3B","f");
  //legend->AddEntry(histb,"3 global muons 2B1M","f");
  //legend->AddEntry(histc,"3 global muons 2B1E","f");
  legend->AddEntry(histd,"3 global muons 3M","f");
  //legend->AddEntry(histe,"3 global muons 2M1B","f");
  //legend->AddEntry(histf,"3 global muons 2M1E","f");
  legend->AddEntry(histg,"3 global muons 3E","f");
  //legend->AddEntry(histh,"3 global muons 2E1B","f");
  //legend->AddEntry(histi,"3 global muons 2E1M","f");
  legend->AddEntry(histj,"3 global muons 1B1M1E","f");
  legend->SetTextSize(0.02);
  legend->Draw();
  

  c0->Print("invariant mass.png");
  cout<<count1<<" passed all cuts" <<count2<<" BBB "<<count3<<" BBM "<<count4<<" BBE "<<count5<<" MMM "<<count6<<" MMB "<<count7<<" MME "<<count8<<" EEE "<<count9<<" EEB "<<count10<<" EEM "<<count11<<" BME "<<count12<<endl;

  gBenchmark->Show("select3Mu");
}

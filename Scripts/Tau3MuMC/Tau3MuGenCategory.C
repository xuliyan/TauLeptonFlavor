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
#include <TPaveStats.h>
#include <TPad.h>
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
  //if(multicount > 1) cout<<"more than one system have come from same tau"<<endl;
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
  const Bool_t isGENLevel = kTRUE;
  const Double_t PT_CUT_LEAD    = 2.5;
  const Double_t PT_CUT_SUB = 1;
  const Double_t ETA_CUT   = 2.4;
  const Double_t MUON_MASS = 0.105658369;
  
  const Int_t BOSON_ID  = 24;
  const Int_t LEPTON_ID = 13;

  TH1F* hista = new TH1F("tau -> 3 muon BBB","invariant mass",200,0,20);
  TH1F* histb = new TH1F("tau -> 3 muon BBE","invariant mass",200,0,20);
  TH1F* histc = new TH1F("tau -> 3 muon BEE","invariant mass",200,0,20);
  TH1F* histd = new TH1F("tau -> 3 muon EEE","invariant mass",200,0,20);
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0, count7=0;
  gStyle->SetOptStat(111111);

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
        infoBr->GetEntry(ientry);
	
        if(ientry%50000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;
	Int_t taumuon = 0;
	Int_t Nummuon[50];
	Int_t runcount = 0;
	Bool_t passSel=kFALSE;

	//===========================BEGIN OF GEN LEVEL======================
	if(isGENLevel){
	  genPartArr->Clear();
	  genPartBr->GetEntry(ientry);
	  //Loop over genParArr to find muon
	  for(Int_t i=0; i<genPartArr->GetEntriesFast(); i++) {
	    const baconhep::TGenParticle *genpar = (baconhep::TGenParticle*)((*genPartArr)[i]);
	    
	    //Select muon
	    if(genpar->pdgId != 13 && genpar->pdgId != -13) continue;
	    if(genpar->status != 1) continue;
            Int_t parentid1=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genpar->parent>-1 ? genpar->parent : 0))->pdgId;
	    if(parentid1 != 15 && parentid1 != -15) continue;

	    taumuon++;
	    Nummuon[runcount] = i;
	    runcount++;
	    
	  }//End of muon array loop 1
	}//============================END OF GEN LEVEL================================
	
	if (taumuon > 2){
	  count1++;
	  passSel=kTRUE;//only events with 3 muons can pass the selection
	}
	
	
	if(passSel){	  
	  if (isGENLevel){
	    vector<baconhep::TGenParticle *> TauMuonArr;
	    //const baconhep::TGenParticle* GenMuonArr[muon];
	    const baconhep::TGenParticle *SelMuonArr[6];
	    Int_t record[3]={0,0,0};
	    vector<TLorentzVector> vLep;
	    TLorentzVector vLepRe[3], vCom;

	    //store global muon into array
	    for(int p=0; p<taumuon; p++){
	      TLorentzVector vtemp;
	      TauMuonArr.push_back((baconhep::TGenParticle*)((*genPartArr)[Nummuon[p]])); //store muon
	      vtemp.SetPtEtaPhiM(TauMuonArr[p]->pt, TauMuonArr[p]->eta, TauMuonArr[p]->phi, MUON_MASS);
	      vLep.push_back(vtemp);//store corresponding 4-vector
	    }
	    
	    baconhep::TGenParticle *genpartarr = (baconhep::TGenParticle*)(genPartArr);
	    if(Findtaudecay(TauMuonArr,vLep,record,taumuon,genpartarr)){     
	      for(int p=0; p<3; p++){
	      SelMuonArr[p] = TauMuonArr[record[p]];
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
	      if(ientry <250)
		cout<<SelMuonArr[3]->pt<<" "<<SelMuonArr[4]->pt<<" "<<SelMuonArr[5]->pt<<endl;
	      for(int p=3; p<6; p++){
		vLepRe[p-3].SetPtEtaPhiM(SelMuonArr[p]->pt, SelMuonArr[p]->eta, SelMuonArr[p]->phi, MUON_MASS);
	      }
	      vCom = vLepRe[0]+vLepRe[1]+vLepRe[2];
	      Double_t invmass = vCom.M();
	      count2++;

	      bool passcut = kTRUE;
	      /*
	      if(fabs(SelMuonArr[0]->eta) > ETA_CUT) passcut=kFALSE;  // lepton |eta| cut
	      if(fabs(SelMuonArr[1]->eta) > ETA_CUT) passcut=kFALSE;  // lepton |eta| cut
	      if(fabs(SelMuonArr[2]->eta) > ETA_CUT) passcut=kFALSE;  // lepton |eta| cut
	      if(passcut) count2++;
	      if(GenMuonArr[0]->pt<PT_CUT_LEAD || GenMuonArr[1]->pt<PT_CUT_LEAD || GenMuonArr[2]->pt<PT_CUT_SUB) passcut=kFALSE;  // lepton pT cut
	      */
	     
	      if(passcut){
		Int_t Bcount = 0;
		if (fabs(SelMuonArr[3]->eta)<0.8) Bcount++;
		if (fabs(SelMuonArr[4]->eta)<0.8) Bcount++;
		if (fabs(SelMuonArr[5]->eta)<0.8) Bcount++;
		if (Bcount == 3){
		  count3++;
		  hista->Fill(invmass);}//store events come from tau
		else if(Bcount == 2){
		  count4++;
		  histb->Fill(invmass);}
		else if(Bcount == 1){
		  count5++;
		  histc->Fill(invmass);}
		else if(Bcount == 0){
		  count6++;
		  histd->Fill(invmass);}
	      }
	    }
	  }
        }
      }//end of event loop
       cout<<"There is "<<passpresel<<" events passed pre selection"<<endl;
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
  cout << "  pT > " << PT_CUT_LEAD << endl;
  cout << "  |eta| < " << ETA_CUT << endl;
  cout << endl;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;

   //Draw Graph
  TCanvas *c0 = new TCanvas("c0","invariant mass",1200,900);
  TAxis *xaxis = histd->GetXaxis();
  TAxis *yaxis = histd->GetYaxis();

  xaxis->SetTitle("Invariant mass (GeV)");
  yaxis->SetTitle("Entries / 0.1 GeV");
  yaxis->SetTitleOffset(1.2);
  //yaxis->SetRangeUser(0.5,150);
  c0->cd();
  c0->SetLogy();
  

  histd->SetLineColor(4);
  histc->SetLineColor(6);
  histb->SetLineColor(8);
  hista->SetLineColor(2);


  histd->Draw();
  gPad->Update();
  /* collect stat of the first histogram (h1) */
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
  
  auto legend = new TLegend(0.5,0.7,0.7,0.8);
  legend->AddEntry(hista,"3 global muons BBB","f");
  legend->AddEntry(histb,"3 global muons BBE","f");
  legend->AddEntry(histc,"3 global muons BEE","f");
  legend->AddEntry(histd,"3 global muons EEE","f");
  legend->Draw();
  

  c0->Print("invariant mass.png");
  cout<<count1<<" passed all cuts" <<count2<<" BBB "<<count3<<" BBE "<<count4<<" BEE"<<count5<<" EEE "<<count6<<endl;

  gBenchmark->Show("select3Mu");

}

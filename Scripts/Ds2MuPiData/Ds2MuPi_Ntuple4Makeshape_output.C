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
  TH1F* TriMuFra = new TH1F("TriMuFra","Fraction of trigger object",51,0,1.02);
  TH1F* NumGlobal = new TH1F("NUM","Number of global muon",10,0,10);
  TH1F* hist0 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #pi^{#pm} 0","m_{#mu^{+}#mu^{-}} (Data)",100,0,5);
  TH1F* hist1 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #pi^{#pm} 1","pT",100,0,20);
  TH1F* hist2 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #pi^{#pm} 2","#mu^{-} pT",100,0,20);
  TH1F* hist3 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #pi^{#pm} 3","pT",100,0,20);
  TH1F* hist4 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #pi^{#pm} 4","#eta",100,-5,5);
  TH1F* hist5 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #pi^{#pm} 5","#mu^{-} #eta",100,-5,5);
  TH1F* hist6 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #pi^{#pm} 6","#pi^{-} #eta",100,-5,5);
  TH1F* hist7 = new TH1F("Ds^{#pm} -> #mu^{+} #mu^{-} #Chi^{#pm} 7","m_{#mu^{+}#mu^{-}#Chi^{#pm}} (Data)",200,0,10);
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
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_Ds2MuPi_select.root");
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
      //or(UInt_t ientry=0; ientry<100000; ientry++) {
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
	
	//Loop over muon collection
	muonArr->Clear();
        muonBr->GetEntry(ientry);
	vector<baconhep::TMuon*> MuArr;
	vector<baconhep::TMuon*> AntiMuArr;
	vector<Int_t> SeqnumMu; //num in muonArr
	vector<Int_t> SeqnumAnti;
	Int_t TriMuNum = 0;
	Int_t TrkMuNum = 0;
        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          baconhep::TMuon *muon = (baconhep::TMuon*)((*muonArr)[i]);
	  if(!(muon->typeBits & baconhep::EMuType::kTracker)) continue; //tracker muon
	  TrkMuNum++;
	  if(!isMuonTriggerObj(triggerMenu, muon->hltMatchBits, kFALSE)) continue; //trigger obj
	  TriMuNum++;
	  if(muon->q == -1) {
	    MuArr.push_back(muon);
	    SeqnumMu.push_back(i);
	  }
	  if(muon->q == 1) {
	    AntiMuArr.push_back(muon);
	    SeqnumAnti.push_back(i);
	  }
	}//end of muon loop;

	if(MuArr.size() < 1 || AntiMuArr.size() < 1) continue;
	count5++;

        TriMuFra->Fill((Float_t)TriMuNum/(Float_t)TrkMuNum); //tracker muon trigger object fraction

	//Select muon pair
	baconhep::TMuon *mu1 = NULL;
	baconhep::TMuon *mu2 = NULL;
	TLorentzVector vtemp1, vtemp2, vmu1, vmu2;
	Double_t massdiffm = 9999;
	Int_t seqnum[2] = {0,0};
	for(int i=0; i<MuArr.size(); i++){
	  vtemp1.SetPtEtaPhiM(MuArr[i]->pt, MuArr[i]->eta, MuArr[i]->phi, MUON_MASS);
	  for(int j=0; j<AntiMuArr.size(); j++){
	    vtemp2.SetPtEtaPhiM(AntiMuArr[j]->pt, AntiMuArr[j]->eta, AntiMuArr[j]->phi, MUON_MASS);
	    Double_t invmass = (vtemp1+vtemp2).M();
	    if(fabs(invmass - 1.019445) < massdiffm){
	      massdiffm = fabs(invmass - 1.019445);
	      mu1 = MuArr[i];
	      mu2 = AntiMuArr[j];
	      vmu1 = vtemp1;
	      vmu2 = vtemp2;
	      seqnum[0] = SeqnumMu[i];
	      seqnum[1] = SeqnumAnti[j];
	    }
	  }
	}
	if(massdiffm > 0.06) continue; //exclude muon paires that have invmass being outside of [0.959,1.079]

	//Select track
	baconhep::TMuon *trk1 = NULL;
	TLorentzVector vtemp3, vtrk1;
	Double_t massdiffm3 = 9999;	
	for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
	  if(i == seqnum[0] || i == seqnum[1]) continue;
          baconhep::TMuon *trk = (baconhep::TMuon*)((*muonArr)[i]);
	  if(trk->typeBits != 0) continue; //Exclude muon
	  vtemp3.SetPtEtaPhiM(trk->pt, trk->eta, trk->phi, 0.13957);
	  Double_t invmass3 = (vtemp1+vtemp2+vtemp3).M();
	  if(fabs(invmass3 - 1.96847) < massdiffm3){
	    massdiffm3 = fabs(invmass3 - 1.96847);
	    trk1 = trk;
	    vtrk1 = vtemp3;
	  }
	}
	hist0->Fill((vmu1+vmu2).M());
	hist7->Fill((vmu1+vmu2+vtrk1).M());
	hist1->Fill(mu1->pt);
	hist2->Fill(mu2->pt);
	hist3->Fill(trk1->pt);
	hist4->Fill(mu1->eta);
	hist5->Fill(mu2->eta);
	hist6->Fill(trk1->eta);

	//Fill tree
	sysinvmass = (vmu1+vmu2+vtrk1).M();
	//*v1 = vmu1;
	//*v2 = vmu2;
	//*v3 = vtrk1;
	if(sysinvmass > 1.67682 || sysinvmass < 1.87682) continue; //Exclude signal region 5 sigma -- 100MeV around tau mass
	outTree->Fill();
	
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

  //Draw
  //gStyle->SetOptStat(0);
  TCanvas *c0 = new TCanvas("1","1",1200,900);
  TAxis *xaxis = hist0->GetXaxis();
  TAxis *yaxis = hist0->GetYaxis();
  xaxis->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
  yaxis->SetTitle("Entries / 50 MeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.5,50000);
  //c0->SetLogx();
  c0->SetLogy();
  c0->cd();
  hist0->SetFillColor(38);
  hist0->Draw();
  auto legend = new TLegend(0.50,0.76,0.65,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist0,"#mu^{+}#mu^{-}","f");
  legend->Draw();
  c0->Print("invmassdimu.png");

  TCanvas *c1 = new TCanvas("2","2",1200,900);
  xaxis = hist3->GetXaxis();
  yaxis = hist3->GetYaxis();
  xaxis->SetTitle("pT (GeV)");
  yaxis->SetTitle("Entries / 0.2 GeV");
  yaxis->SetTitleOffset(1.2);
  c1->cd();
  c1->SetLogy();
  hist1->SetLineColor(2);
  hist2->SetLineColor(4);
  hist3->SetLineColor(6);
  hist1->SetFillStyle(0);
  hist2->SetFillStyle(0);
  hist3->SetFillStyle(0);
  hist3->Draw();
  hist1->Draw("same");
  hist2->Draw("same");
  legend = new TLegend(0.50,0.72,0.6,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist1,"#mu^{+}","f");
  legend->AddEntry(hist2,"#mu^{-}","f");
  legend->AddEntry(hist3,"#pi^{#pm}","f");
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
  hist6->SetLineColor(6);
  hist4->SetFillStyle(0);
  hist5->SetFillStyle(0);
  hist6->SetFillStyle(0);
  hist4->Draw();
  hist5->Draw("same");
  hist6->Draw("same");
  legend = new TLegend(0.45,0.50,0.55,0.62);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist4,"#mu^{+}","f");
  legend->AddEntry(hist5,"#mu^{-}","f");
  legend->AddEntry(hist6,"#pi^{#pm}","f");
  legend->Draw();
  c2->Print("eta.png");

  TCanvas *c3 = new TCanvas("4","4",1200,900);
  xaxis = hist7->GetXaxis();
  yaxis = hist7->GetYaxis();
  xaxis->SetTitle("m_{#mu^{+}#mu^{-}#Chi^{#pm}} (GeV)");
  yaxis->SetTitle("Entries / 50 MeV");
  yaxis->SetTitleOffset(1.2);
  yaxis->SetRangeUser(0.5,50000);
  //c0->SetLogx();
  c3->SetLogy();
  c3->cd();
  hist7->SetFillColor(38);
  hist7->Draw();
  legend = new TLegend(0.50,0.76,0.65,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist7,"#mu^{+}#mu^{-} #Chi^{#pm}","f");
  legend->Draw();
  c3->Print("invmassmutrk.png");

  TCanvas *c4 = new TCanvas("5","5",1200,900);
  xaxis = TriMuFra->GetXaxis();
  yaxis = TriMuFra->GetYaxis();
  xaxis->SetTitle("Fraction of trigger object");
  yaxis->SetTitle("Entries");
  yaxis->SetTitleOffset(1.4);
  //yaxis->SetRangeUser(0.5,50000);
  //c0->SetLogx();
  //c4->SetLogy();
  c4->cd();
  TriMuFra->Draw();
  legend = new TLegend(0.11,0.81,0.47,0.87);
  legend->AddEntry(TriMuFra,"hltL2fL1sL1DoubleMuorTripleMuL1f0L2PreFiltered0","f");
  legend->SetTextSize(0.02);
  legend->Draw();
  c4->Print("triobjfra.png");
  
  gBenchmark->Show("selectDsPhiPi"); 
}

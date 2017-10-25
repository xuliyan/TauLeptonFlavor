#include <TMath.h>
void make_Tau3Mu_shapes(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  // --- Create workspace --- 
  RooWorkspace w("w","w");
  RooRealVar sysinvmass("sysinvmass","invmass",1.57682,1.97682,""); //the name "sysinvmass" will be used by RooDataSet to import data
  sysinvmass.setBins(40);
  //signal PDF
  RooRealVar* mu1 = new RooRealVar("DG_mu1","#mu1",1.96847,1.8,2.2,"");
  RooRealVar* mu2 = new RooRealVar("DG_mu2","#mu2",1.96847,1,3,"");
  RooRealVar* sigma1 = new RooRealVar("DG_sigma1","#sigma_{1}",0,0.05,"");
  RooRealVar* sigma2 = new RooRealVar("DG_sigma2","#sigma_{2}",0,0.1,"");
  mu1->setConstant(kFALSE);
  mu2->setConstant(kFALSE);
  sigma1->setConstant(kFALSE);
  sigma2->setConstant(kFALSE);
  RooRealVar* frac = new RooRealVar("DG_frac","frac",0.5,0,1); 
  RooRealVar* Ns = new RooRealVar("DG_Ns","N_{s}",1000,0,10000,"events");
  Ns->setConstant(kFALSE);
  RooGaussian* gauss1 = new RooGaussian("G1","",sysinvmass,*mu1,*sigma1);
  RooGaussian* gauss2 = new RooGaussian("G2","",sysinvmass,*mu2,*sigma2);
  RooAddPdf* doublegauss = new RooAddPdf("SummedG1G2","",RooArgList(*gauss1,*gauss2),*frac);
  RooAddPdf* ex_doublegauss = new RooAddPdf("DG_ex","extDgauss",RooArgList(*doublegauss),RooArgList(*Ns));
  w.import(*ex_doublegauss);
  //background
  RooRealVar *pC = new RooRealVar("pol2_pC","C",0.5,0,1,"");
  pC->setConstant(kFALSE);
  RooRealVar *p0 = new RooRealVar("pol2_p0","p_0",0.3,0,1,"");
  p0->setConstant(kFALSE);
  RooRealVar *p1 = new RooRealVar("pol2_p1","p_1",0.27,0,1,"");
  p1->setConstant(kFALSE);
  RooRealVar *p2 = new RooRealVar("pol2_p2","p_2",0.5,0,1,"");
  p1->setConstant(kFALSE);
  RooRealVar *p3 = new RooRealVar("pol2_p3","p_3",0.5,0,1,"");
  p1->setConstant(kFALSE);
  RooRealVar *Nbkg   = new RooRealVar("pol2_Nbkg","N_{bkg}",1000,0,10000,"events");
  Nbkg->setConstant(kFALSE);
  RooBernstein* bern = new RooBernstein("Bern","",sysinvmass,RooArgList(*pC,*p0,*p1));
  RooAddPdf* ex_bern = new RooAddPdf("Bern_ex","",RooArgList(*bern),RooArgList(*Nbkg));
  w.import(*ex_bern);
  w.factory("SUM::model_s(nB[1000,0,10000]*Bern, nS[1000,0,10000]*SummedG1G2)");

  // --- Import unbinned dataset ---
  TFile file("/afs/cern.ch/user/x/xuyan/3MuonProj/CMSSW_8_0_27/src/DataFlat/ntuples/data_Tau3Mu_select.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(sysinvmass);
  RooArgSet imargset(imarglist);
  RooDataSet data("data","data",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform extended ML fit of composite PDF to toy data ---
  //w.pdf("model_s")->fitTo(data);

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = sysinvmass.frame();
  data.plotOn(frame);
  frame->SetTitle("m_{#mu^{+}#mu^{-}#Chi^{#pm}}");
  //w.pdf("model_s")->plotOn(frame,LineColor(kRed));
  //w.pdf("model_s")->plotOn(frame,Components("Bern"),LineStyle(kDashed));
  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{#mu^{+}#mu^{-}#Chi^{#pm}} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 10 MeV");
  yaxis->SetTitleOffset(1.2);
  frame->SetMaximum(3000);
  frame->SetMinimum(5);
  c1->SetLogy();
  frame->Draw();
  c1->Print("data.png");
  
  // --- Output root file ---
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sysinvmass[1.57682,1.97682]");
  wUP->import(data,Rename("data_obs"));
  wUP->import(*doublegauss);
  wUP->import(*bern);
  wUP->writeToFile("Tau3Mu-shapes-UnbinnedParam.root");
}

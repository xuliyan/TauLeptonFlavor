#include <TMath.h>
#include <TLegend.h>
void make_tau3mu_shapes(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  // --- Create workspace --- 
  RooWorkspace w("w","w");
  RooRealVar sysinvmass("sysinvmass","invmass",1.07682,2.47682,""); //the name "sysinvmass" will be used by RooDataSet to import data
  sysinvmass.setBins(70);

  //--- signal PDF ---
  RooRealVar* mu1 = new RooRealVar("DG_mu1","#mu1",1.77682,1.75,1.8,"");
  RooRealVar* mu2 = new RooRealVar("DG_mu2","#mu2",1.77682,1.5,2,"");
  RooRealVar* sigma1 = new RooRealVar("DG_sigma1","#sigma_{1}",0,0.05,"");
  RooRealVar* sigma2 = new RooRealVar("DG_sigma2","#sigma_{2}",0,0.1,"");
  mu1->setConstant(kFALSE);
  mu2->setConstant(kFALSE);
  sigma1->setConstant(kFALSE);
  sigma2->setConstant(kFALSE);
  RooRealVar* frac = new RooRealVar("DG_frac","frac",0.5,0,1); 
  RooRealVar* Ns = new RooRealVar("DG_Ns","N_{s}",3000,0,10000,"events");
  Ns->setConstant(kFALSE);
  RooGaussian* gauss1 = new RooGaussian("G1","",sysinvmass,*mu1,*sigma1);
  RooGaussian* gauss2 = new RooGaussian("G2","",sysinvmass,*mu2,*sigma2);
  RooAddPdf* doublegauss = new RooAddPdf("SummedG1G2","",RooArgList(*gauss1,*gauss2),*frac);
  RooAddPdf* ex_doublegauss = new RooAddPdf("DG_ex","extDgauss",RooArgList(*doublegauss),RooArgList(*Ns));
  w.import(*ex_doublegauss);

  //--- background PDF ---
  RooRealVar *pC = new RooRealVar("pol2_pC","C",0.5,0,1,"");
  pC->setConstant(kFALSE);
  RooRealVar *p0 = new RooRealVar("pol2_p0","p_0",0.3,0,1,"");
  p0->setConstant(kFALSE);
  RooRealVar *p1 = new RooRealVar("pol2_p1","p_1",0.27,0,1,"");
  p1->setConstant(kFALSE);
  RooRealVar *p2 = new RooRealVar("pol2_p2","p_2",0.5,0,1,"");
  p2->setConstant(kFALSE);
  RooRealVar *p3 = new RooRealVar("pol2_p3","p_3",0.5,0,1,"");
  p3->setConstant(kFALSE);
  RooRealVar *p4 = new RooRealVar("pol2_p4","p_4",0.5,0,1,"");
  p4->setConstant(kFALSE);
  RooRealVar *Nbkg   = new RooRealVar("pol2_Nbkg","N_{bkg}",2000,0,10000,"events");
  Nbkg->setConstant(kFALSE);
  RooBernstein* bern = new RooBernstein("Bern","",sysinvmass,RooArgList(*pC,*p0,*p1));
  RooAddPdf* ex_bern = new RooAddPdf("Bern_ex","",RooArgList(*bern),RooArgList(*Nbkg));
  w.import(*ex_bern);

  //--- combine PDF ---
  RooAddPdf* ex_sum = new RooAddPdf("SUM_ex","",RooArgList(*bern,*doublegauss),RooArgList(*Nbkg,*Ns));
  w.import(*ex_sum);
  
  // --- Import unbinned dataset ---
  TFile file("/afs/cern.ch/user/x/xuyan/3MuonProj/CMSSW_8_0_27/src/Ntuples/data_Tau3Mu_select_highpT_triobjcut.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(sysinvmass);
  RooArgSet imargset(imarglist);
  RooDataSet data("data","data",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform extended ML fit of composite PDF to toy data ---
  //w.pdf("model_s")->fitTo(data);
  //sysinvmass.setRange("SB1",1.47682,1.67682);
  //sysinvmass.setRange("SB2",1.87682,2.07682);
  //RooFitResult *r = ex_bern->fitTo(data,Range("SB1,SB2"),Save());
  //RooFitResult *r = ex_doublegauss->fitTo(data,Save());
  //RooFitResult *r = gauss1->fitTo(data,Save());

  // --- Plot ---
  gStyle->SetOptStat(111111);
  RooPlot *frame = sysinvmass.frame();
  data.plotOn(frame);
  frame->SetTitle("m_{#mu^{#mp}#mu^{#pm}#mu^{#mp}} Data");
  //w.pdf("model_s")->plotOn(frame,LineColor(kRed));
  //ex_doublegauss->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //gauss1->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //ex_bern->plotOn(frame,Components("Bern"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  //ex_doublegauss->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  //gauss1->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  //ex_doublegauss->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //gauss1->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //ex_doublegauss->plotOn(frame,Components("Bern"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  //data.plotOn(frame);
  //w.pdf("model_s")->plotOn(frame,Components("Bern"),LineStyle(kDashed));

  //axis,log scale and range setting functions must be called after all plotOn functions being called
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("m_{#mu^{#mp}#mu^{#pm}#mu^{#mp}} (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 20 MeV");
  yaxis->SetTitleOffset(1.2);
  frame->SetMaximum(800);
  frame->SetMinimum(0.5);
  //c1->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
  l->AddEntry(frame->findObject("sigfun"),"Signal Fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err"),"Fit Error 1 #sigma","f");
  //l->Draw("same");
  c1->Print("data.png");
  
  // --- Output root file ---
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sysinvmass[1.07682,2.47682]");
  wUP->import(data,Rename("data_obs"));
  wUP->import(*doublegauss);
  wUP->import(*bern);
  wUP->writeToFile("Tau3Mu-shapes-UnbinnedParam.root");
}

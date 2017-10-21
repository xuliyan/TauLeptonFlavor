#include <TMath.h>
void make_Ds2mupi_shapes_1(int seed=37)
{
  using namespace RooFit;

  // ---Define PDFs ---
  RooRealVar sysinvmass("sysinvmass","invmass",1.5,2.5,""); 
  //signal PDF
  RooRealVar* mu1 = new RooRealVar("DG_mu1","#mu1",1.96847,1.5,2.5,"");
  RooRealVar* mu2 = new RooRealVar("DG_mu2","#mu2",1.96847,1,3,"");
  RooRealVar* sigma1 = new RooRealVar("DG_sigma1","#sigma_{1}",0.005,0.05,"");
  RooRealVar* sigma2 = new RooRealVar("DG_sigma2","#sigma_{2}",0.1,0.5,"");
  mu1->setConstant(kFALSE);
  mu2->setConstant(kFALSE);
  sigma1->setConstant(kFALSE);
  sigma2->setConstant(kFALSE);
  RooRealVar* frac = new RooRealVar("DG_frac","frac",0.5,0,1); 
  RooGaussian* gauss1 = new RooGaussian("G1","",sysinvmass,*mu1,*sigma1);
  RooGaussian* gauss2 = new RooGaussian("G2","",sysinvmass,*mu2,*sigma2);
  RooAddPdf* doublegauss = new RooAddPdf("SummedG1G2","",RooArgList(*gauss1,*gauss2),*frac);
  //background
  RooRealVar *pC = new RooRealVar("pol2_pC","C",0.5,0,1,"");
  pC->setConstant(kFALSE);
  RooRealVar *p0 = new RooRealVar("pol2_p0","p_0",0.3,0,1,"");
  p0->setConstant(kFALSE);
  RooRealVar *p1 = new RooRealVar("pol2_p1","p_1",0.27,0,1,"");
  p1->setConstant(kFALSE);
  RooBernstein* bern = new RooBernstein("Bern","",sysinvmass,RooArgList(*pC,*p0,*p1));
  
  // --- Import unbinned dataset ---
  TFile file("/afs/cern.ch/user/x/xuyan/3MuonProj/CMSSW_8_0_27/src/DataFlat/ntuples/data_Ds2MuPi_select.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(sysinvmass);
  RooArgSet imargset(imarglist);
  RooDataSet data("data","data",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset
  
  // --- Output root file ---
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("sysinvmass[1.5,2.5]");
  wUP->import(data,Rename("data_obs"));
  wUP->import(*doublegauss);
  wUP->import(*bern);
  wUP->writeToFile("Ds2mupi-shapes-UnbinnedParam.root");
}

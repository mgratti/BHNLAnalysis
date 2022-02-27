#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
using namespace RooFit;

// Create a workspace with toy data, which is readible by the flashgg final fit package
// Had to use C++ because, model.generate() does not seem to work in ROOT 6.12 in python
// Other functions were working fine

void makeExampleWS(){

  // Set the seed
  RooRandom::randomGenerator()->SetSeed(1234);

  //
  int n = 2000;
  TString foutName = "ws/allData.root";
  double mass_low = 2.75;  // 10 sigma for mass=3 GeV
  double mass_high = 3.25; //

  // output vars 
  RooRealVar CMS_hgg_mass("CMS_hgg_mass","CMS_hgg_mass",mass_low,mass_high);
  //RooRealVar intLumi("intLumi","intLumi",0.,1e+09);
  //RooRealVar weight("weight","weight",-1e+09,1e+09); 

  std::vector<TString> cats = {"Lxy0to1_OS", "Lxy1to5_OS", "LxyGt5_OS"};
  
  // Create output
  TFile *fout = new TFile(foutName, "RECREATE"); 
  fout->mkdir("tagsDumper");
  fout->cd("tagsDumper");
  RooWorkspace *w = new RooWorkspace("cms_hgg_13TeV", "cms_hgg_13TeV");

  for (int i=0; i<cats.size(); i++){
    // Toy data generation
    RooRealVar pow1("pow1","exponent of power law",-3,-6,-0.0001);
    RooGenericPdf bkg("powerlaw","TMath::Power(@0,@1)",RooArgList(CMS_hgg_mass,pow1));
    RooDataSet *data = bkg.generate(CMS_hgg_mass, n);

    //RooRealVar a0("a0", "a0", -0.5-0.1*i, 0.0, 1.0);
    //RooRealVar a1("a1", "a1",  0.2,       0.0, 1.0);
    //RooChebychev bkg("bkg", "Background", CMS_hgg_mass, RooArgList(a0, a1));
    //model = ROOT.RooAddPdf("bkg","bkg",ROOT.RooArgList(bkg));
    //model = ROOT.RooAddPdf("model", "g1+g2+a", [bkg, sig], [bkgfrac])
 
    TString dname = "Data_13TeV_" + cats.at(i);
    data->SetName(dname);
    data->SetNameTitle(dname,dname);
    data->Print("v");

    w->import(*data);

    // also create binned datasets
    CMS_hgg_mass.setBins(1000);
    RooDataHist hdata_1000("hdata_1000","binned version of data, 1000 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_1000);
    CMS_hgg_mass.setBins(500);
    RooDataHist hdata_500("hdata_500","binned version of data, 500 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_500);
    CMS_hgg_mass.setBins(320);
    RooDataHist hdata_320("hdata_320","binned version of data, 320 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_320);
    CMS_hgg_mass.setBins(300);
    RooDataHist hdata_300("hdata_300","binned version of data, 300 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_300);
    CMS_hgg_mass.setBins(100);
    RooDataHist hdata_100("hdata_100","binned version of data, 100 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_100);
    CMS_hgg_mass.setBins(10);
    RooDataHist hdata_10("hdata_10","binned version of data, 10 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_10);
    CMS_hgg_mass.setBins(2);
    RooDataHist hdata_2("hdata_2","binned version of data, 2 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_2);
    CMS_hgg_mass.setBins(2);
    RooDataHist hdata_1("hdata_1","binned version of data, 1 bins",RooArgSet(CMS_hgg_mass),*data) ;
    w->import(hdata_1);

  }
  
  w->Print();
  w->Write();
  fout->Close();

  w->Delete();
  fout->Delete();

}

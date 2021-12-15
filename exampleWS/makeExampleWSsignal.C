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

// Had to use C++ because, model.generate() does not seem to work in ROOT 6.12 in python
// Other functions were working fine


// TODO: add dZ variable...

void makeExampleWSsignal(){

  //
  int n = 500;
  TString foutName = "ws/output_M3_BHNL.root";
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
    
    // Fake data generation
    RooRealVar mean("mean", "mean", 3., 0.0, 6.0);
    RooRealVar sigma("sigma", "sigma",  0.025,0.0, 0.5);
    RooGaussian signal("signal", "signal", CMS_hgg_mass, mean, sigma);
  
    RooDataSet *data = signal.generate(CMS_hgg_mass, n);
    TString dname = "process_signal_" + cats.at(i);
    data->SetName(dname);
    data->SetNameTitle(dname,dname);
    data->Print("v");

    w->import(*data);

  }
  
  w->Print();
  w->Write();
  fout->Close();

  w->Delete();
  fout->Delete();

}

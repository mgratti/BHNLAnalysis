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

void makeExampleWSsignal(){

  //
  double mass_low = 2.625;  // 5 sigma for mass=2.75 GeV
  double mass_high = 2.875;  //
  double m = 2.75;
  int nb = 100;
  TString foutName = "ws/output_M2p75_BHNL.root";

  // output vars 
  RooRealVar hnl_mass("hnl_mass","hnl_mass",mass_low,mass_high);
  //RooRealVar dZ("dZ","dZ",-20,20);
  //RooRealVar intLumi("intLumi","intLumi",0.,1e+09);
  //RooRealVar weight("weight","weight",-1e+09,1e+09); 

  //std::vector<TString> cats = {"Lxy0to1_OS", "Lxy1to5_OS", "LxyGt5_OS"};
  
  // Create signal model 
  RooRealVar mean("mean", "mean", m, 0.0, 6.0);
  mean.setConstant(kTRUE);
  RooRealVar sigma("sigma", "sigma",  0.025,0.0, 0.5);
  sigma.setConstant(kTRUE);
  RooGaussian gauss("gauss", "gauss", hnl_mass, mean, sigma);
  
  // Create output
  TFile *fout = new TFile(foutName, "RECREATE"); 
  RooWorkspace *w = new RooWorkspace("ws", "ws");

  /*
  for (int i=0; i<cats.size(); i++){
    
    // Fake data generation
    RooRealVar mean("mean", "mean", 3., 0.0, 6.0);
    RooRealVar sigma("sigma", "sigma",  0.025,0.0, 0.5);
    RooGaussian signal("signal", "signal", hnl_mass, mean, sigma);
  
    RooDataSet *data = signal.generate(hnl_mass, n);
    //RooDataSet *data = signal.generate(RooArgSet(hnl_mass,dZ),n);  

    // Fit data with Gaussian, and save the PDF

    TString dname = "process_signal_" + cats.at(i);
    data->SetName(dname);
    data->SetNameTitle(dname,dname);
    data->Print("v");

    w->import(*data);

  }
  */

  w->import(gauss); 
  w->import(mean);
  w->import(sigma);
  
  w->Print();
  w->Write();
  fout->Close();

  w->Delete();
  fout->Delete();

}

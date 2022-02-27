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

// Final script usef to create workspace compatible with flashggfinalfit from data trees

void makeAnalysisWS(float mHNL=2.75,float sigma=0.025,int nsigma=5){

  //double mass_low = 2.625;  // 5 sigma for mass=2.75 GeV
  //double mass_high = 2.875;  //
  //double mHNL = 2.75;
  float mass_low = mHNL - nsigma*sigma;
  float mass_high = mHNL + nsigma*sigma;
  int nbins = 100;

  // input
  TString fname = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/flat_bparknano_30Dec21.root";
  //TString fname = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/Chunk68_n500/flat/flat_bparknano_30Dec21.root";
  //TString fname = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/Chunk68_n500/flat/flat_bparknano_30Dec21_nj496.root";
  TString tname = "signal_tree";
  TFile f = TFile(fname, "READ");
  TTree *tree = (TTree*)f.Get(tname);

  TString preselection = "trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge))";

  std::map<TString,TString> cats_def;
  std::map<TString,TString> cats_sel;

  std::vector<TString> cats = {"lxy0to1_OS", "lxy1to5_OS", "lxygt5_OS", "lxy0to1_SS", "lxy1to5_SS", "lxygt5_SS" };

  cats_def["lxy0to1_OS"] = "sv_lxy<=1 && trgmu_charge!=mu_charge";
  cats_sel["lxy0to1_OS"] = "hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10";
  
  cats_def["lxy1to5_OS"] = "sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge";
  cats_sel["lxy1to5_OS"] = "hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25";
  
  cats_def["lxygt5_OS"] = "sv_lxy>5 && trgmu_charge!=mu_charge";
  cats_sel["lxygt5_OS"] = "hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20";
  
  cats_def["lxy0to1_SS"] = "sv_lxy<=1 && trgmu_charge==mu_charge";
  cats_sel["lxy0to1_SS"] = "hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10";
  
  cats_def["lxy1to5_SS"] = "sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge";
  cats_sel["lxy1to5_SS"] = "hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25";
  
  cats_def["lxygt5_SS"] = "sv_lxy>5 && trgmu_charge==mu_charge";
  cats_sel["lxygt5_SS"] = "hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20";

  // input vars 
  RooRealVar hnl_mass = RooRealVar("hnl_mass","hnl_mass",mass_low,mass_high);
  RooRealVar trgmu_softid = RooRealVar("trgmu_softid","trgmu_softid",0,1);
  RooRealVar mu_looseid = RooRealVar("mu_looseid","mu_looseid",0,1);
  RooRealVar pi_packedcandhashighpurity = RooRealVar("pi_packedcandhashighpurity","pi_packedcandhashighpurity",0,1);
  RooRealVar trgmu_mu_mass = RooRealVar("trgmu_mu_mass","trgmu_mu_mass",0,13000);
  RooRealVar sv_lxy = RooRealVar("sv_lxy","sv_lxy",0.,1500.);
  RooRealVar trgmu_charge = RooRealVar("trgmu_charge","trgmu_charge",-1,1);
  RooRealVar mu_charge = RooRealVar("mu_charge","mu_charge",-1,1);
  RooRealVar hnl_charge = RooRealVar("hnl_charge","hnl_charge",-1,1);
  RooRealVar pi_pt = RooRealVar("pi_pt","pi_pt",0.,13000.);
  RooRealVar sv_lxysig = RooRealVar("sv_lxysig","sv_lxysig",0.,150000.);
  RooRealVar mu_dxysig = RooRealVar("mu_dxysig","mu_dxysig",0.,150000.);
  RooRealVar pi_dxysig = RooRealVar("pi_dxysig","pi_dxysig",0.,150000.);

  RooArgSet mvars;
  mvars.add(hnl_mass);
  mvars.add(trgmu_softid);
  mvars.add(mu_looseid);
  mvars.add(pi_packedcandhashighpurity);
  mvars.add(trgmu_mu_mass);
  mvars.add(sv_lxy);
  mvars.add(trgmu_charge);
  mvars.add(mu_charge);
  mvars.add(hnl_charge);
  mvars.add(pi_pt);
  mvars.add(sv_lxysig);
  mvars.add(mu_dxysig);
  mvars.add(pi_dxysig);

  // create output
  TString foutName = TString::Format("wsAnalysis/allData_m%.2f_s%.3f_ns%d.root",mHNL,sigma,nsigma);
  
  TFile *fout = new TFile(foutName, "RECREATE"); 
  RooWorkspace *w = new RooWorkspace("cms_hgg_13TeV", "cms_hgg_13TeV");
  fout->mkdir("tagsDumper");
  fout->cd("tagsDumper");

  // Now do the hard part, fill the dataset from the tree
  for( auto cat: cats){

    TString selection = preselection + " && " + cats_def[cat] + " && " + cats_sel[cat];

    // Dataset approach: it loads all the dataset in memory and takes something like one hour (maybe do it in the evening)
    /*
    TString dname = "Data_13TeV_" + cat;
    RooDataSet data(dname,dname, tree, mvars, selection);
    data.Print("v");
    // TODO: add RooDataHist
  
    data.get()->find("hnl_mass")->SetName("CMS_hgg_mass");
    data.get()->find("CMS_hgg_mass")->SetNameTitle("CMS_hgg_mass","CMS_hgg_mass");
    data.Print("v");
    w->import(data);
    */

    // Histogram approach
    TH1D *h = new TH1D("h", "h", nbins, mass_low, mass_high);
    tree->Draw("hnl_mass>>h", selection, "goff");
    TString dname = "Data_13TeV_" + cat;
    RooDataHist rdh = RooDataHist(dname, dname, RooArgList(hnl_mass), h);
    //rdh.get()->find("hnl_mass")->SetName("CMS_hgg_mass");
    //rdh.get()->find("CMS_hgg_mass")->SetNameTitle("CMS_hgg_mass","CMS_hgg_mass");
    rdh.Print("v");
    w->import(rdh);
    


  }

  
  w->Print();
  w->Write();
  fout->Close();


}

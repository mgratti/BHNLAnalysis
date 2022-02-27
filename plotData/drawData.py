'''
Script to plot the data in appropriate windows as shown in 
https://indico.cern.ch/event/1120717/contributions/4705643/attachments/2383967/4075277/2022_02_03_Bpark_v2.pdf
'''

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import RooFit as RF
from glob import glob
from collections import OrderedDict
import sys
sys.path.append('/work/mratti/plotting/myplotting')
from spares import defaultLabels

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.StatOverflows(ROOT.kTRUE) # consider overflows for mean and rms calculation
ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
ROOT.gROOT.ProcessLine('setTDRStyle()')
ROOT.gStyle.SetTitleXOffset(1.1);
ROOT.gStyle.SetTitleYOffset(1.45);


fname = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/merged/flat_bparknano_30Dec21.root"
#fname = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/Chunk68_n500/flat/flat_bparknano_30Dec21.root"
#fname = "/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/V10_30Dec21/ParkingBPH1_Run2018D/Chunk68_n500/flat/flat_bparknano_30Dec21_nj496.root"
tname = "signal_tree"
f=ROOT.TFile(fname, 'READ')
t=f.Get(tname)

preselection = "trgmu_softid == 1 && mu_looseid == 1 && pi_packedcandhashighpurity == 1 && ((trgmu_charge!=mu_charge && (trgmu_mu_mass < 2.9 || trgmu_mu_mass > 3.3)) || (trgmu_charge==mu_charge))"

cats_def =  OrderedDict() 
cats_sel =  OrderedDict() 

cats_def["lxy0to1_OS"] = "sv_lxy<=1 && trgmu_charge!=mu_charge"
cats_sel["lxy0to1_OS"] = "hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10"

cats_def["lxy1to5_OS"] = "sv_lxy>1 && sv_lxy<=5 && trgmu_charge!=mu_charge"
cats_sel["lxy1to5_OS"] = "hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25"

cats_def["lxygt5_OS"] = "sv_lxy>5 && trgmu_charge!=mu_charge"
cats_sel["lxygt5_OS"] = "hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20"

cats_def["lxy0to1_SS"] = "sv_lxy<=1 && trgmu_charge==mu_charge"
cats_sel["lxy0to1_SS"] = "hnl_charge==0 && pi_pt>1.1 && sv_lxysig>30 && abs(mu_dxysig)>5 && abs(pi_dxysig)>10"

cats_def["lxy1to5_SS"] = "sv_lxy>1 && sv_lxy<=5 && trgmu_charge==mu_charge"
cats_sel["lxy1to5_SS"] = "hnl_charge==0 && pi_pt>1.2 && sv_lxysig>100 && abs(mu_dxysig)>12 && abs(pi_dxysig)>25"

cats_def["lxygt5_SS"] = "sv_lxy>5 && trgmu_charge==mu_charge"
cats_sel["lxygt5_SS"] = "hnl_charge==0 && pi_pt>1.3 && sv_lxysig>100 && abs(mu_dxysig)>15 && abs(pi_dxysig)>20"



mass_windows = [
#nbins, xlow, xup
(20,0.8,1.2),
(40,1.2,2.0),
(40,2.0,2.8),
(15,2.8,3.4),
(55,3.4,5.6),
]


for icat in cats_def.keys():

  selection = preselection + ' && ' + cats_def[icat] + ' && ' + cats_sel[icat] 

  # full range plot
  h = ROOT.TH1F("h", "h", 200, 0.8, 5.6)
  t.Draw('hnl_mass>>h', selection, 'goff')

  c2 = ROOT.TCanvas("","", 400,400)
  c2.cd()
  h.GetXaxis().SetTitle('m_{#mu#pi} (GeV)')
  h.GetYaxis().SetTitle('Entries / ({} GeV)'.format((5.6-0.8)/200.))
  h.Draw("histE")
  h.SetMarkerStyle(1)
  label='L = 5.3 fb^{-1}'
  defaultLabels([label],0.80,0.94)
  c2.SaveAs("./plots/data_fullRange_{}.pdf".format(icat))
  c2.SaveAs("./plots/data_fullRange_{}.png".format(icat))
  c2.SaveAs("./plots/data_fullRange_{}.root".format(icat))

  # mass window plots
  for (nbins,mlow,mhigh) in mass_windows:

    k = ROOT.TH1F("k", "k", nbins, mlow, mhigh)
    t.Draw('hnl_mass>>k', selection, 'goff')
  
    c = ROOT.TCanvas("c", "c", 400, 400)
    k.Draw("histE")
    k.SetMarkerStyle(1)
    defaultLabels([label],0.80,0.94)
    k.GetXaxis().SetTitle('m_{#mu#pi} (GeV)')
    k.GetYaxis().SetTitle('Entries / ({} GeV)'.format((mhigh-mlow)/float(nbins)))
    c.SaveAs('./plots/data_m{}To{}_{}.pdf'.format(mlow,mhigh,icat))
    c.SaveAs('./plots/data_m{}To{}_{}.png'.format(mlow,mhigh,icat))


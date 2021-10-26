import ROOT
import os
import sys
import ratioplot as RP
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *
# Enable multi-threading
#ROOT.ROOT.EnableImplicitMT()
from quantity import Quantity,quantities_to_plot_all,quantities_to_plot_gen



if __name__ == "__main__":

  #### ROOT Options
  gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()
  ROOT.TH1.StatOverflows(ROOT.kTRUE) # consider overflows for mean and rms calculation
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gStyle.SetTitleXOffset(1.1);
  gStyle.SetTitleYOffset(1.45);


  path_low = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/mc_central/test/BToNMuX_NToEMuPi_SoftQCD_b_mN1p0_ctau10p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_n1_f1.root'
  path_high = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/mc_central/test/BToNMuX_NToEMuPi_SoftQCD_b_mN1p5_ctau100p0mm_TuneCP5_13TeV-pythia8-evtgen/merged/flat_bparknano_n1_f1.root'
  f_low = ROOT.TFile.Open(path_low)
  f_high = ROOT.TFile.Open(path_high)
  treename = 'signal_tree'
  runtname = 'run_tree'
  outDir = 'plots/1708_central_mainGrid'
  try: os.mkdir(outDir)
  except: pass
  
  ctau_low = 1.0 # cm
  mass_low = 1.0 # GeV
  ctau_high = 10.0 # cm
  mass_high = 1.5 # GeV

  t_low = f_low.Get(treename) 
  t_high = f_high.Get(treename)
  rt_low = f_low.Get(runtname)
  rt_high = f_high.Get(runtname)


  qs = quantities_to_plot_all + quantities_to_plot_gen

#  qs = [
#   { 'what':'sv_lxy', 'name':'Lxy', 'title':'L_{xy} (cm)', 'binning':'(50,0.,100.)', 'log':True },
#   { 'what':'sv_lxy', 'name':'Lxy', 'title':'L_{xy} (cm)', 'binning':'(50,0.,100.)', 'log':False },
#   { 'what':'hnl_mass', 'name':'hnlmass', 'title':'m_{#mu #pi} (GeV)', 'binning':'(100,2.,4.)', 'log': False},
#   #{ 'what':'hnl_fitted_mass', 'name':'hnlmass', 'title':'fitted m_{#mu #pi} (GeV)', 'binning':'(50,0.,5.)', 'log': False},
#  ]

  histo_saver = []

  # selections ? no
  # weights ? no

  for q in qs:
    
    t_high.Draw('%s>>%s%s' % (q.name_flat, q.label+'_high',      '({},{},{})'.format(q.nbins,q.bin_min,q.bin_max)),  '(1)' ,  'goff')
    h_high = t_high.GetHistogram()
    t_low.Draw( '%s>>%s%s' % (q.name_flat, q.label+'_low',       '({},{},{})'.format(q.nbins,q.bin_min,q.bin_max)),  '(1)' ,  'goff')
    h_low = t_low.GetHistogram()

    if h_low: print 'Entries, h_low={}'.format(h_low.GetEntries())
    if h_high: print 'Entries, h_high={}'.format(h_high.GetEntries())

    histo_saver.append(h_low) 
    histo_saver.append(h_high) 

    RP.makeRatioPlot( hNum=h_low, 
                      hDen=h_high, 
                      hDen2="", 
                      nameNum="m={}GeV,c#tau={}cm".format(mass_low, ctau_low), 
                      nameDen="m={}GeV,c#tau={}cm".format(mass_high,ctau_high), 
                      nameDen2="", 
                      xtitle=q.title,
                      ytitle="a.u.", 
                      ratiotitle="Ratio", 
                      norm=True, 
                      log=False, 
                      plotName="ratio_{}_woLog_wNorm".format(q.label), 
                      outDir=outDir,ratioyrange=(0.,4.) )
#    RP.makeRatioPlot( hNum=h_low, 
#                      hDen=h_high, 
#                      hDen2="", 
#                      nameNum="m={}GeV,c#tau={}cm".format(mass_low, ctau_low), 
#                      nameDen="m={}GeV,c#tau={}cm".format(mass_high,ctau_high), 
#                      nameDen2="", 
#                      xtitle=q.title,
#                      ytitle="a.u.", 
#                      ratiotitle="Ratio", 
#                      norm=True, 
#                      log=True, 
#                      plotName="ratio_{}_wLog_wNorm".format(q.label), 
#                      outDir=outDir,ratioyrange=(0.,4.) )


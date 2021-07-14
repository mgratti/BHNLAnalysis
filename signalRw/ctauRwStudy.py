import ROOT
import os
import sys
import ratioplot as RP
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *
# Enable multi-threading
#ROOT.ROOT.EnableImplicitMT()
from quantity import Quantity,quantities_to_plot_all



if __name__ == "__main__":

  #### ROOT Options
  gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()
  ROOT.TH1.StatOverflows(ROOT.kTRUE) # consider overflows for mean and rms calculation
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gStyle.SetTitleXOffset(1.1);
  gStyle.SetTitleYOffset(1.45);


  path_low  = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V20_emu/mass3.0_ctau184.0/nanoFiles/merged/flat_bparknano_n1_f1.root'
  path_high = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V25/mass3.0_ctau2000.0/nanoFiles/merged/flat_bparknano_n1_f1.root'
  f_low = ROOT.TFile.Open(path_low)
  f_high = ROOT.TFile.Open(path_high)
  treename = 'signal_tree'
  runtname = 'run_tree'
  outDir = 'plots/V20_emu_vs_V25'
  try: os.mkdir(outDir)
  except: pass
  
  ctau_low = 18.4 # cm
  ctau_high = 200. # cm

  t_low = f_low.Get(treename) 
  t_high = f_high.Get(treename)
  rt_low = f_low.Get(runtname)
  rt_high = f_high.Get(runtname)


  qs = quantities_to_plot_all

#  qs = [
#   { 'what':'sv_lxy', 'name':'Lxy', 'title':'L_{xy} (cm)', 'binning':'(50,0.,100.)', 'log':True },
#   { 'what':'sv_lxy', 'name':'Lxy', 'title':'L_{xy} (cm)', 'binning':'(50,0.,100.)', 'log':False },
#   { 'what':'hnl_mass', 'name':'hnlmass', 'title':'m_{#mu #pi} (GeV)', 'binning':'(100,2.,4.)', 'log': False},
#   #{ 'what':'hnl_fitted_mass', 'name':'hnlmass', 'title':'fitted m_{#mu #pi} (GeV)', 'binning':'(50,0.,5.)', 'log': False},
#  ]

  histo_saver = []

  # selections ? no
  # weights ? yes, either no weight or cta weight 

  # determine the n gen
  feff_low = 3.93e-03
  feff_high = 1.02e-03

  h_ngen_low = ROOT.TH1D('h_ngen_low', 'h_ngen_low', 1,0.,1.E12)
  rt_low.Draw('geneventsumw>>h_ngen_low', '', 'goff')
  ngen_low = h_ngen_low.GetEntries() * h_ngen_low.GetMean()
  ngen_low /= feff_low
  ngen_low /= 2. # this is due to the fact that low was generated 50% - 50% e, mu, while high was generated 100% muons
  h_ngen_high = ROOT.TH1D('h_ngen_high', 'h_ngen_high', 1,0.,1.E12)
  rt_high.Draw('geneventsumw>>h_ngen_high', '', 'goff')
  ngen_high = h_ngen_high.GetEntries() * h_ngen_high.GetMean()
  ngen_high /= feff_high
  print 'ngen_low={}, ngen_high={}'.format(ngen_low,ngen_high)

  ctau_weight = '({old_ctau}/{new_ctau} * exp( (1./{old_ctau} - 1./{new_ctau}) * gen_hnl_ct))'

  for q in qs:

    weight_to_low  = ctau_weight.format(old_ctau=ctau_high, new_ctau=ctau_low)
    weight_to_high = ctau_weight.format(old_ctau=ctau_low, new_ctau=ctau_high)
    
    t_high.Draw('%s>>%s%s' % (q.name_flat, q.label+'_highToLow', '({},{},{})'.format(q.nbins/2,q.bin_min,q.bin_max)),  weight_to_low + '/' + str(ngen_high),  'goff')
    h_highToLow = t_high.GetHistogram()
    t_high.Draw('%s>>%s%s' % (q.name_flat, q.label+'_high',      '({},{},{})'.format(q.nbins/2,q.bin_min,q.bin_max)),  '(1)'         + '/' + str(ngen_high),  'goff')
    h_high = t_high.GetHistogram()
    t_low.Draw( '%s>>%s%s' % (q.name_flat, q.label+'_lowToHigh', '({},{},{})'.format(q.nbins/2,q.bin_min,q.bin_max)),  weight_to_high+ '/' + str(ngen_low),   'goff')
    h_lowToHigh = t_low.GetHistogram()
    t_low.Draw( '%s>>%s%s' % (q.name_flat, q.label+'_low',       '({},{},{})'.format(q.nbins/2,q.bin_min,q.bin_max)),  '(1)'         + '/' + str(ngen_low),   'goff')
    h_low = t_low.GetHistogram()

    if h_low: print 'Entries, h_low={}'.format(h_low.GetEntries())
    if h_high: print 'Entries, h_high={}'.format(h_high.GetEntries())

    # print the yields
    print '\nHigh to Low ctau'
    print 'rw ={:.2e}'.format(h_highToLow.Integral())
    print 'gen={:.2e}'.format(h_low.Integral())
    print 'df ={:.1f}%'.format(abs(h_highToLow.Integral()/h_low.Integral()-1)*100)
    
    print '\nLow to High ctau'
    print 'rw ={:.2e}'.format(h_lowToHigh.Integral())
    print 'gen={:.2e}'.format(h_high.Integral())
    print 'df ={:.1f}%'.format(abs(h_lowToHigh.Integral()/h_high.Integral()-1)*100)


    histo_saver.append(h_low) 
    histo_saver.append(h_lowToHigh) 
    histo_saver.append(h_high) 
    histo_saver.append(h_highToLow) 

    ## ratio plot for reweighting from High To Low 
#    RP.makeRatioPlot( hNum=h_low, 
#                      hDen=h_highToLow, 
#                      hDen2="", 
#                      nameNum="gen c#tau={}cm".format(ctau_low), 
#                      nameDen="rw  c#tau={}cm, orig. c#tau={}cm".format(ctau_low,ctau_high), 
#                      nameDen2="", 
#                      xtitle=q.title,
#                      ytitle="a.u.", 
#                      ratiotitle="gen/rw", 
#                      norm=True, 
#                      log=q['log'], 
#                      plotName="ratio_{}_HighToLow_{}Log_wNorm".format(q.label, 'w' if q['log'] else 'wo'), 
#                      outDir=outDir,ratioyrange=(0.4,1.6) )

    RP.makeRatioPlot( hNum=h_low, 
                      hDen=h_highToLow, 
                      hDen2="", 
                      nameNum="gen c#tau={}cm".format(ctau_low), 
                      nameDen="rw  c#tau={}cm, orig. c#tau={}cm".format(ctau_low,ctau_high), 
                      nameDen2="", 
                      xtitle=q.title,
                      ytitle="a.u.", 
                      ratiotitle="gen/rw", 
                      norm=False, 
                      log=False, 
                      plotName="ratio_{}_HighToLow_woLog_woNorm".format(q.label), 
                      outDir=outDir,ratioyrange=(0.4,1.6) )

    ## ratio plot for reweighting from Low To High
#    RP.makeRatioPlot( hNum=h_high, 
#                      hDen=h_lowToHigh, 
#                      hDen2="", 
#                      nameNum="gen c#tau={}cm".format(ctau_high), 
#                      nameDen="rw  c#tau={}cm, orig. c#tau={}cm".format(ctau_high,ctau_low), 
#                      nameDen2="", 
#                      xtitle=q.title,
#                      ytitle="a.u.", 
#                      ratiotitle="gen/rw", 
#                      norm=True, 
#                      log=q['log'], 
#                      plotName="ratio_{}_LowToHigh_{}Log_wNorm".format(q.label, 'w' if q['log'] else 'wo'), 
#                      outDir=outDir,ratioyrange=(0.4,1.6) )
    RP.makeRatioPlot( hNum=h_high, 
                      hDen=h_lowToHigh, 
                      hDen2="", 
                      nameNum="gen c#tau={}cm".format(ctau_high), 
                      nameDen="rw  c#tau={}cm, orig. c#tau={}cm".format(ctau_high,ctau_low), 
                      nameDen2="", 
                      xtitle=q.title,
                      ytitle="a.u.", 
                      ratiotitle="gen/rw", 
                      norm=False, 
                      log=False, 
                      plotName="ratio_{}_LowToHigh_woLog_woNorm".format(q.label), 
                      outDir=outDir,ratioyrange=(0.4,1.6) )


    # compare low and high directly: how different are they in the first place?
    RP.makeRatioPlot( hNum=h_low,
                      hDen=h_high,
                      hDen2="",
                      nameNum="gen c#tau={}cm".format(ctau_low), 
                      nameDen="gen c#tau={}cm".format(ctau_high), 
                      nameDen2="", 
                      xtitle=q.title,
                      ytitle="a.u.", 
                      ratiotitle="gen_low/gen_high", 
                      norm=True, 
                      log=False, 
                      plotName="ratio_{}_LowFracHigh_woLog_wNorm".format(q.label), 
                      outDir=outDir,ratioyrange=(0.4,1.6) )


    break

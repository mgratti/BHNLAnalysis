'''
Script to calculate the inclusive B->Jpsi(mumu)K channel cross-section
Please see https://indico.cern.ch/event/1099857/contributions/4627416/attachments/2353189/4020897/2021_11_30_Bpark_SignalNorm.pdf
'''


import ROOT
from ROOT import RooFit as RF
import sys
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *
from glob import glob
#import sys
#sys.path.append('/work/mratti/plotting/myplotting')
#from spares import *


from argparse import ArgumentParser
parser = ArgumentParser(description='Production helper for B-initiated HNL signals', add_help=True)
parser.add_argument('--whichSignalModel', type=str, dest='whichSignalModel', help='', default='voigt')
parser.add_argument('--whichPRecBkgModel', type=str, dest='whichPRecBkgModel', help='', default='cb')
parser.add_argument('--doNewSelection', dest='doNewSelection', help='', action='store_true', default=False)
parser.add_argument('--doFiducial', dest='doFiducial', help='', action='store_true', default=False)
parser.add_argument('--doFidOneBin', dest='doFidOneBin', help='', action='store_true', default=False)
parser.add_argument('--whichCutRemoved', type=int, dest='whichCutRemoved', help='', default=None)
opt = parser.parse_args()



mass_xlow   = 5.1 # for signal fit only
mass_xhigh  = 5.6 # for signal fit only
doSignalFit = False
doFullFit   = True
whichSignalModel  = opt.whichSignalModel #'voigt' # gauss, voigt, dscb, gauss2Sum
whichCombBkgModel = 'exp' # poly1, exp
whichPRecBkgModel = opt.whichPRecBkgModel # argus argusPlusGauss argusConvGauss exp cb
doAddJpsiPi = False
doNewFiles = True
doNewSelection = opt.doNewSelection
doHLTWeights = True
doFiducial = opt.doFiducial
doFidOneBin = opt.doFidOneBin
deltaDiMu = 0.05
doDrawParameters = False
doDrawLegend = True


labelmodel={}
labelmodel['gauss2Sum'] = 'Sum of two Gaussians'
labelmodel['gauss'] = 'Gaussian'
labelmodel['voigt'] = 'Voigtian'
labelmodel['dscb'] = 'Double-sided Crystal Ball'

def getTotalEffFiducial(nSelMC):

  b_fiducial_cuts = 'b_pt > 10 && b_pt < 17 && abs(b_y) < 1.45'
  if doFidOneBin:
    b_fiducial_cuts = 'b_pt > 10 && b_pt < 13 && abs(b_y) < 1.45'
  gen_filter_cuts = 'l1_pt > 6.8 && abs(l1_eta) < 1.55'
  
  ### Epsilon 2 
  #### from 'unfiltered' sample (i.e. no muon filter, but B filter with fiducial)
  path = '/work/mratti/GEN_HNL/CMSSW_10_2_15/src/HNLsGen/genLevelAnalysis/outputfiles/'
  pl = 'V15_control_BfilterNoMufilter/'
  fname = 'mass999_ctau999_miniGenTree.root'
  f = ROOT.TFile.Open(path+pl+fname,'read')
  t = f.Get('tree')
  
  num_cuts = b_fiducial_cuts + ' && ' + gen_filter_cuts
  den_cuts = b_fiducial_cuts
  
  hnum = ROOT.TH1F('hnum', 'hnum', 1, 0., 20.)
  hden = ROOT.TH1F('hden', 'hden', 1, 0., 20.)
  
  t.Draw('b_pt>>hnum', num_cuts, 'goff')
  t.Draw('b_pt>>hden', den_cuts, 'goff')
  
  
  peff = ROOT.TEfficiency(hnum,hden)
  eff2,eff2_errup,eff2_errdn = peff.GetEfficiency(1), peff.GetEfficiencyErrorUp(1), peff.GetEfficiencyErrorLow(1)
  
  print('===> Resulting efficiency Epsilon2')
  print('eff2 = {:.4f} + {:.4f} - {:.4f}'.format(eff2,eff2_errup,eff2_errdn))
  
  
  ### Epsilon 1, 
  #### denominator from filtered sample (i.e. with muon filter)
  #### numerator from filtered sample, but at RECO level
  pl = 'V15_control/'
  fname = 'mass999_ctau999_miniGenTree.root'
  f = ROOT.TFile.Open(path+pl+fname,'read')
  t = f.Get('tree')
   
  den_cuts = b_fiducial_cuts + ' && ' + gen_filter_cuts
  
  hnum = ROOT.TH1F('hnum', 'hnum', 1, 0., 20.)
  hden = ROOT.TH1F('hden', 'hden', 1, 0., 20.)
  
  ##treco.Draw('b_pt>>hnum', num_cuts, 'godff')
  t.Draw('b_pt>>hden', den_cuts, 'goff')
  hnum.SetBinContent(1,nSelMC)
  hnum.SetBinError(1,ROOT.TMath.Sqrt(nSelMC))  

  peff = ROOT.TEfficiency(hnum,hden)
  eff1,eff1_errup,eff1_errdn = peff.GetEfficiency(1), peff.GetEfficiencyErrorUp(1), peff.GetEfficiencyErrorLow(1) 
 
  print('===> Resulting efficiency Epsilon1')
  print('eff1 = {:.4f} + {:.4f} - {:.4f}'.format(eff1,eff1_errup,eff1_errdn))

  toteff = eff1*eff2
  toteff_err = ROOT.TMath.Sqrt((eff1_errup/eff1)**2 + (eff2_errup/eff2)**2) * toteff

  return toteff,toteff_err



def getChiSquare(fitmodel,RDSet):
 
  flparams = fitmodel.getParameters(RDSet).selectByAttrib('Constant',ROOT.kFALSE) 
  nflparams = flparams.getSize()
  my_frame = mass.frame()
  RDSet.plotOn(my_frame, RF.DataError(ROOT.RooAbsData.Poisson))
  nmeas = my_frame.GetNbinsX()
  ndof = - nflparams + nmeas  

  fitmodel.plotOn(my_frame) 
  my_chi2 = my_frame.chiSquare(nflparams)

  from scipy.stats import chi2
  prob = 1-chi2.cdf(my_chi2*ndof, ndof)

  #chi2_var = ROOT.RooChi2Var('chi2_var', 'chi2_var', fitmodel, RDSet, RF.DataError(ROOT.RooAbsData.Poisson))
  #print 'Chi2 from chi var', chi2_var.printValue()
  

  print ('\n****************************\n')
  print ('Fit quality information')
  print ('nflparams = {}'.format(nflparams))
  print ('nmeas     = {}'.format(nmeas))
  print ('ndof      = {}'.format(ndof))
  print ('chi2      = {:.1f}'.format(my_chi2))
  print ('prob      = {:.4f}'.format(prob))
  print ('\n****************************\n')


  return my_chi2,prob

def drawPlot(frame,frame2,chisq,prob,sigmaBpm,sigmaBpmErr,lumi,label='',leg=None):

  hpull = frame.pullHist()
  hpull.SetMarkerSize(0.7)
  frame2.addPlotable(hpull,'P');

  c = ROOT.TCanvas('', '', 400, 400)
  c.Draw()
  c.Divide(1,2)
  #ROOT.gPad.SetLeftMargin(0.15)
  #frame.GetYaxis().SetTitleOffset(1.6)
  c.cd(1)
  ROOT.gPad.SetLeftMargin(0.15)
  ROOT.gPad.SetBottomMargin(0.10)
  ROOT.gPad.SetPad(0.01,0.2,0.99,0.99)
  frame.GetXaxis().SetTitleSize(0.04);
  frame.GetXaxis().SetTitle('m_{K J/#psi} (GeV)')
  frame.GetYaxis().SetTitleSize(0.04);
  frame.GetYaxis().SetTitle('Entries')
  frame.GetYaxis().SetTitleOffset(1.35)
  frame.Draw()

  if doDrawLegend and leg and not doDrawParameters:
    leg.Draw("same")     

  labelchi2 = '#chi^{{2}}/n_{{dof}} = {:.1f}'.format(chisq)
  labelprob = '  p-value = {:.3f}'.format(prob)
  if sigmaBpm is not None and lumi is not None:
    labelsigma = '#sigma(B^{{\pm}}) = ({:.1f} #pm {:.1f}) x 10^{{9}} fb ({})'.format(sigmaBpm/1E9,sigmaBpmErr/1E9, 'inclusive' if not doFiducial else 'fiducial')
    labellumi = 'L = {:.3f} fb^{{-1}}'.format(lumi)
    #defaultLabels([labelchi2+labelprob,labelsigma,labellumi], 0.55, 0.23) #spacing = 0.04, size = 0.027, dx = 0.12)
    defaultLabels([labelchi2+labelprob,labelsigma], 0.55, 0.23) #spacing = 0.04, size = 0.027, dx = 0.12)
    defaultLabels([labellumi], 0.82,0.93)
  else: 
    defaultLabels([labelmodel[whichSignalModel],labelchi2+labelprob], 0.60, 0.25)

  c.cd(2)
  frame2.Draw()
  ROOT.gPad.SetLeftMargin(0.15)
  ROOT.gPad.SetPad(0.01,0.01,0.99,0.2)
  ROOT.gPad.SetGridy()
  frame2.GetYaxis().SetNdivisions(504)
  frame2.GetYaxis().SetLabelSize(0.17)
  frame2.GetYaxis().SetTitleSize(0.17)
  frame2.GetYaxis().SetTitleOffset(0.24)
  frame2.GetYaxis().SetRangeUser(-3,3)        
  frame2.GetYaxis().SetTitle('Pulls') 
  frame2.GetXaxis().SetTitle('')      
  frame2.GetXaxis().SetLabelOffset(5)          
           
  c.SaveAs('fit{}.pdf'.format(label))




if __name__ == '__main__':

  #### ROOT Options
  ROOT.gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()
  ROOT.TH1.StatOverflows(ROOT.kTRUE) # consider overflows for mean and rms calculation
  #ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  ROOT.gROOT.ProcessLine('.L ./tdrstyle.C')
  ROOT.gROOT.ProcessLine('setTDRStyle()')
  ROOT.gStyle.SetTitleXOffset(1.1);
  ROOT.gStyle.SetTitleYOffset(1.45);


  #############
  # Data import
  #############

  if doFiducial:
    mult = 1./25. # scale the starting points for nsig and nbkg, when applying the fiducial cuts
  else: 
    mult = 1.


  if doNewFiles:
    # version of October 
    file_mc = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_Oct20_Oct20.root'
    files_data_periodA = ['/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/data/Control_Oct20/ParkingBPH1_Run2018A/merged/flat_bparknano_Oct20_TEST.root']
    # for data 1.61% of the jobs failed
 
  else:
    # these are the files that were shared with Ludovico
    file_mc =  '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_forRenormalisationStudy.root'
    files_data_periodA = glob('/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/data/F1/ParkingBPH1_Run2018A/merged/flat_bparknano_forRenormalisationStudy.root')

  
  # tree
  treename = 'control_tree'
  
  # define the chains from where to import the data
  tree_data = ROOT.TChain(treename)
  for fname in files_data_periodA:
    tree_data.AddFile(fname)
  tree_mc = ROOT.TChain(treename)
  tree_mc.AddFile(file_mc)
  
  # vars
  mass = ROOT.RooRealVar('b_mass','m_{#mu#mu K}',5.0,6.0,'GeV')
  b_cos2d = ROOT.RooRealVar('b_cos2d', 'b_cos2d', 0,1)
  b_pt = ROOT.RooRealVar('b_pt', 'b_pt', 0.,13000.)
  k_pt = ROOT.RooRealVar('k_pt', 'k_pt', 0.,13000.)
  l1_pt = ROOT.RooRealVar('l1_pt', 'l1_pt', 0.,13000.)
  l2_pt = ROOT.RooRealVar('l2_pt', 'l2_pt', 0.,13000.)
  b_eta = ROOT.RooRealVar('b_eta', 'b_eta', -10., 10.)
  k_eta = ROOT.RooRealVar('k_eta', 'k_eta', -10., 10.)
  l1_eta = ROOT.RooRealVar('l1_eta', 'l1_eta', -10., 10.)
  l2_eta = ROOT.RooRealVar('l2_eta', 'l2_eta', -10., 10.)
  dimu_mass = ROOT.RooRealVar('dimu_mass', 'dimu_mass',0,100,'GeV')
  sv_lxy = ROOT.RooRealVar('sv_lxy', 'sv_lxy', 0, 200, 'cm')
  sv_prob = ROOT.RooRealVar('sv_prob', 'sv_prob', 0., 1.)
  hlt_mu9_ip6 = ROOT.RooRealVar('hlt_mu9_ip6', 'hlt_mu9_ip6', 0, 1)
  l1_softid = ROOT.RooRealVar('l1_softid', 'l1_softid', 0,1)
  l2_softid = ROOT.RooRealVar('l2_softid', 'l2_softid', 0,1)
  sv_lxysig = ROOT.RooRealVar('sv_lxysig', 'sv_lxysig', 0,1E09)
  #l1_softid = ROOT.ROOT.RooRealVar('', '', 0,1)
  # nice to haves: impact parameters (?), significance of the displacement
  if doNewFiles:
    b_y = ROOT.RooRealVar('b_y', 'b_y', -10., 10.)
    weight_hlt = ROOT.RooRealVar('weight_hlt', 'weight_hlt', 0., 1.)
    #matched_b_y = ROOT.RooRealVar('matched_b_y', 'matched_b_y', -10., 10.) # do not uncomment, still need to understand why, FIXME
    #matched_b_pt = ROOT.RooRealVar('matched_b_pt', 'matched_b_pt', 0.,13000.)
    #sv_lxysig = ROOT.RooRealVar('sv_lxysig', 'sv_lxysig', 0., 100.)
    #dimu_sv_prob = ROOT.RooRealVar('dimu_sv_prob', 'dimu_sv_prob', 0., 1.)
    
  
  # mc
  ismatched = ROOT.RooRealVar('ismatched', 'ismatched', 0,1)
  
  # set of vars
  mvars = ROOT.RooArgSet()
  mvars.add(mass)
  mvars.add(b_cos2d)
  mvars.add(b_pt)
  mvars.add(k_pt)
  mvars.add(l1_pt)
  mvars.add(l2_pt)
  mvars.add(b_eta)
  mvars.add(k_eta)
  mvars.add(l1_eta)
  mvars.add(l2_eta)
  mvars.add(dimu_mass)
  mvars.add(sv_lxy)
  mvars.add(sv_prob)
  mvars.add(hlt_mu9_ip6)
  if doNewFiles:
    mvars.add(b_y)
    mvars.add(l1_softid)
    mvars.add(l2_softid)
    mvars.add(sv_lxysig)
    #mvars.add(matched_b_y) # do not uncomment!
    #mvars.add(matched_b_y)
  
  mvars_mc = mvars
  if doNewFiles and doHLTWeights:
    mvars_mc.add(weight_hlt)
  #mvars_mc.add(ismatched)
  
  # nice to haves: impact parameters (?), significance of the displacement
  
  # mc
  ismatched = ROOT.RooRealVar('ismatched', 'ismatched', 0,1)
  
  # selection
  sel = '(1)'
  #selBase = 'b_pt>30 && abs(b_eta)<1.45 && b_cos2d>0.995 && dimu_mass > (3.097-0.20) && dimu_mass < (3.097+0.20) && b_mass>5 && sv_lxy>0.1 && k_pt>1.5 && abs(k_eta)<2.0 && abs(l2_eta)<2.0'
  #selBase = 'b_cos2d>0.995 && dimu_mass > (3.097-0.20) && dimu_mass < (3.097+0.20) && b_mass>5 && sv_lxy>0.1 && k_pt>1.5 && abs(k_eta)<2.0 && abs(l2_eta)<2.0'

  ## baseline (dimu window 0.20 GeV: can go tighter?, like 0.05?)
  #selBase = 'b_cos2d>0.995 && dimu_mass > (3.097-0.20) && dimu_mass < (3.097+0.20) && b_mass>5 && sv_lxy>0.1 && k_pt>1.5 && abs(k_eta)<2.0 && abs(l2_eta)<2.0 && l1_pt > 4. && l2_pt > 4.'
  
  ## baseline + trigger cut (only for cross-checking purposes)
  my_old_selections = [
    'b_mass > 5 ',
    'b_cos2d > 0.995', 
    'abs(dimu_mass-3.097) < {} '.format(deltaDiMu),
    'sv_lxy > 0.1 ', 
    'sv_prob > 0.', 
    'k_pt > 1.5 ',
    'abs(k_eta) < 2. ',
    #'l1_pt > 9. ',
    'l2_pt > 2. ', 
    'abs(l2_eta) < 2 ',
    'l1_softid==1',
    #'l2_softid==1',
    #'sv_lxysig > 5',
    #'hlt_mu9_ip6==1'
  ]
  ## start from Ludovico's selections and add soft id
  my_new_selections = [
    'b_mass > 5 ',
    'b_cos2d > 0.9995 ',
    'abs(dimu_mass-3.097) < {} '.format(deltaDiMu),
    'sv_lxy > 0.035 ',
    'sv_prob > 0.08 ',
    'k_pt > 1.0 ',
    'abs(k_eta) < 1.6 ',
    ##'l1_pt > 9. ',
    'l2_pt > 2. ',
    'abs(l2_eta) < 1.8 ',
    'l1_softid==1',
    ##'hlt_mu9_ip6==1',
  ]

  label_new_selections = [
    'bMass',
    'cos2D',
    'dimuMass',
    'Lxy',
    'svProb',
    'kPt',
    'kEta',
    ## 'l1Pt',
    'l2Pt',
    'l2Eta',
    'l1Soft',
    ## 'hltMu9',
  ]

  if opt.doNewSelection and opt.whichCutRemoved is not None:
    my_new_selections = my_new_selections[0:opt.whichCutRemoved] + my_new_selections[opt.whichCutRemoved+1:]
    removed_selection_label = label_new_selections[opt.whichCutRemoved]

  ## ludovico's
  ludo_selections = [
    'b_cos2d > 0.9995 ',
    'dimu_mass < (3.097+0.15) ',
    'dimu_mass > (3.097-0.15) ',
    'b_mass > 5 ',
    'sv_lxy > 0.035 ',
    'sv_prob > 0.08 ',
    'k_pt > 1.0 ',
    'abs(k_eta) < 1.6 ',
    #'l1_pt > 9. ',
    'l2_pt > 2. ',
    'abs(l2_eta) < 1.8 ',
    #'hlt_mu9_ip6==1',
  ]

  if doNewFiles and doNewSelection:
    selBase = ' && '.join(my_new_selections)
  else:
    selBase = ' && '.join(my_old_selections)

  ## fiducial cuts
  selFiducial =   'b_pt > 10 && b_pt < 17 && abs(b_y) < 1.45'  
  if doFidOneBin:
    selFiducial = 'b_pt > 10 && b_pt < 13 && abs(b_y) < 1.45' 

  if doFiducial:
    sel = selBase +  ' && ' + selFiducial 
  else:
    sel = selBase

  print sel
 
  # define the datasets
  Rdata = ROOT.RooDataSet('Rdata', 'data', tree_data, mvars, sel) # no weight
  Rmc = ROOT.RooDataSet('Rmc', 'MC', tree_mc, mvars_mc, sel) ## 

  print('before weighting')
  Rmc.Print()
  if doNewFiles and doHLTWeights:
    Rmc = ROOT.RooDataSet('Rmc', 'MC', tree_mc, mvars_mc, sel, weight_hlt.GetName())

  print('after weighting')
  Rmc.Print()

  ## test
  frametest = mass.frame(RF.Title(''))
  Rdata.plotOn(frametest, RF.MarkerSize(0.5), RF.XErrorSize(0), RF.Name('data'))
  ctest = TCanvas()
  frametest.Draw()
  ctest.SaveAs('test.pdf')

  
  #############
  # Fit models
  #############
  
  ## bkg
  nbkg_comb = ROOT.RooRealVar('nbkg_comb','number of bkg events (combinatorial)',            1000*mult,0,1.0E07) 
  nbkg_prec = ROOT.RooRealVar('nbkg_prec','number of bkg events (partially reconstructed)', 40000*mult,0,1.0E07) 
  
  ### 1) bkg from ~5 GeV (JPsi + hadrons), misreconstruction 
  argus_shape = ROOT.RooRealVar('argus_shape','argus shape parameter',-3,-10.,-1.) # -5 ok
  argus_fallm = ROOT.RooRealVar('argus_fallm','argus falling mass',5.279, 5.16,5.32,)
  #argus_fallm = ROOT.RooRealVar('argus_fallm','argus falling mass',5.16, 5.16,5.16,)
  #argus_fallm = ROOT.RooRealVar('argus_fallm','argus falling mass',5.279, 5.279,5.279)
  argus = ROOT.RooArgusBG('argus','argus',mass,argus_fallm,argus_shape)

  # argus with gaussian resolution
  gaussReso_mean =  ROOT.RooRealVar('gaussReso_mean', 'gaussReso_mean', 5, 4., 5.3)
  gaussReso_sigma = ROOT.RooRealVar('gaussReso_sigma', 'gaussReso_sigma', 1., 0., 1.)
  gaussReso = ROOT.RooGaussian('gaussReso', 'gaussReso', mass,gaussReso_mean,gaussReso_sigma) 
  if whichPRecBkgModel == 'argusConvGauss':
    mass.setBins(10000, 'cache')
  argusConvGauss = ROOT.RooFFTConvPdf('argusConvGauss', 'argusConvGauss', mass, argus, gaussReso)
  # argus + gaussian
  argusPlusGauss_frac = ROOT.RooRealVar('argusPlusGauss_frac','argusPlusGauss_frac', 0.5, 0.0 ,1.)
  argusPlusGauss = ROOT.RooAddPdf('argusPlusGauss', 'argusPlusGauss', ROOT.RooArgList(argus,gaussReso), ROOT.RooArgList(argusPlusGauss_frac))
  #argusConvGauss = ROOT.RooNumConvPdf('argusConvGauss', 'argusConvGauss', mass, argus, gaussReso)
  #RooFFTConvPdf lxg("lxg","landau (X) gauss",t,landau,gauss) ;
  
  # crystal ball
  cb_mean   = ROOT.RooRealVar('cb_mean','cb_mean',  5.0, 4.9, 5.1)
  cb_sigma  = ROOT.RooRealVar('cb_sigma','cb_sigma', 0.1, 0., 0.3) 
  cb_alpha  = ROOT.RooRealVar('cb_alpha', 'cb_alpha', 0.5, 0., 2.) # alpha: where the Gaussian and exponential connect, at alpha * sigma of the Gaussian (pos)
  cb_n      = ROOT.RooRealVar('cb_n', 'cb_n', 1., 0., 10.)         # the slope of the exponential 
  cb        = ROOT.RooCBShape('cb', 'cb', mass, cb_mean, cb_sigma, cb_alpha, cb_n)

  ##erf_xshift = ROOT.RooRealVar('erf_xshift', 'xshift', 6., -1000., 1000.) 
  ##erf_yshift = ROOT.RooRealVar('erf_yshift', 'yshift', 1., -1000., 1000.) 
  ##erf = ROOT.RooGenericPdf('erf', 'erf', '-1.*TMath::Erf(mass-erf_xshift)+erf_yshift', ROOT.RooArgList(mass,erf_xshift,erf_yshift))
  # careful, this is not positive definite...
  
  exp1_a = ROOT.RooRealVar('exp1_a', 'exp1_a', -1., -10.,0.)
  exp1 = ROOT.RooExponential('exp1', 'Exponential', mass, exp1_a)
  
  ### 2) combinatorial
  #### exponential
  exp_a = ROOT.RooRealVar('exp_a', 'exp_a', -0.0001, -5.,0.)
  exp = ROOT.RooExponential('exp', 'Exponential', mass, exp_a)
  
  #### polynomials
  poly_a0 = ROOT.RooRealVar('poly_a0', 'poly_a0',   2., -10000., 10000.)
  poly_a1 = ROOT.RooRealVar('poly_a1', 'poly_a1',   1., -10000., 10000.)
  poly_a2 = ROOT.RooRealVar('poly_a2', 'poly_a2', 0.01, -10000., 10000.)
  poly1 = ROOT.RooPolynomial('poly1', 'poly1', mass, ROOT.RooArgList(poly_a0, poly_a1))
  poly2 = ROOT.RooPolynomial('poly2', 'poly2', mass, ROOT.RooArgList(poly_a0, poly_a1, poly_a2))

  ### 3) JPsi pi
  #### gauss
  gaussJpsiPi_mean  = ROOT.RooRealVar('gaussJpsiPi_mean', 'gaussJpsiPi_mean', 5.35, 5.3, 5.45  )
  gaussJpsiPi_sigma  = ROOT.RooRealVar('gaussJpsiPi_sigma', 'gaussJpsiPi_sigma', 0.05, 0.02, 0.1  )
  gaussJpsiPi = ROOT.RooGaussian('gaussJpsiPi', 'gaussJpsiPi', mass, gaussJpsiPi_mean, gaussJpsiPi_sigma)
  nbkg_peak = ROOT.RooRealVar('nbkg_peak','number of bkg events (J/#psi #pi)', 200,0,1.0E07) 
  
  ## signal
  nsig  = ROOT.RooRealVar('nsig','number of signal events', 20000*mult,0,1E07)
  
  ### gauss
  gauss_mean  = ROOT.RooRealVar('gauss_mean', 'gauss_mean', 5.28, 5., 6.)
  gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'gauss_sigma', 0.06, 0., 0.08)
  gauss = ROOT.RooGaussian('gauss', 'gauss', mass,gauss_mean,gauss_sigma)
 
  ### sum of two gaussians
  gauss2_meanDisp = ROOT.RooRealVar('gauss2_meanDisp', 'gauss2_meanDisp', 0., 0., 0.2)
  gauss2Sum_frac01 = ROOT.RooRealVar('gauss2Sum_frac01','gauss2Sum_frac01', 0.8, 0.0 ,0.95)
  gauss2_sigmaMult = ROOT.RooRealVar('gauss2_sigmaMult', 'gauss2_sigmaMult', 0.5, 0.01, 1.)

  #if not doFullFit:
    ##gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'gauss_sigma', 0.06, 0., 0.08)
    #gauss2_sigmaMult = ROOT.RooRealVar('gauss2_sigmaMult', 'gauss2_sigmaMult', 0.5, 0.01, 1.)
    #gauss2_meanDisp = ROOT.RooRealVar('gauss2_meanDisp', 'gauss2_meanDisp', 0., 0., 0.2)
  #  gauss2Sum_frac01 = ROOT.RooRealVar('gauss2Sum_frac01','gauss2Sum_frac01', 0.8, 0.0 ,0.95)
  #else:
    ##gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'gauss_sigma', 0.0375 , 0.0375 , 0.0375 )
    #gauss2_sigmaMult = ROOT.RooRealVar('gauss2_sigmaMult', 'gauss2_sigmaMult', 0.686, 0.686, 0.686)
    #gauss2_meanDisp = ROOT.RooRealVar('gauss2_meanDisp', 'gauss2_meanDisp', 0.0044, 0.0044, 0.0044)
  #  gauss2Sum_frac01 = ROOT.RooRealVar('gauss2Sum_frac01','gauss2Sum_frac01', 0.51, 0.51 ,0.51)
  gauss2_mean =  ROOT.RooFormulaVar('gauss2_mean', 'gauss_mean+gauss2_meanDisp', ROOT.RooArgList(gauss_mean,gauss2_meanDisp) )
  gauss2_sigma = ROOT.RooFormulaVar('gauss2_sigma','gauss_sigma*gauss2_sigmaMult', ROOT.RooArgList(gauss_sigma,gauss2_sigmaMult))
  gauss2 = ROOT.RooGaussian('gauss2', 'gauss2', mass, gauss2_mean, gauss2_sigma)
  gauss2Sum = ROOT.RooAddPdf('gauss2Sum', 'gauss2Sum', ROOT.RooArgList(gauss, gauss2), ROOT.RooArgList(gauss2Sum_frac01))

  ### convolution of three gaussians 
  # TODO
 
  ### RooVoigtian is an efficient implementation of the convolution of a Breit-Wigner with a Gaussian
  voigt_mean  = ROOT.RooRealVar('voigt_mean', 'voigt_mean', 5.28, 5., 6.)
  voigt_sigma = ROOT.RooRealVar('voigt_sigma', 'voigt_sigma', 0.02, 0., 0.06)
  voigt_width = ROOT.RooRealVar('voigt_width', 'voigt_width', 0.002, 0., 0.06)
  voigt       = ROOT.RooVoigtian('voigt', 'voigt', mass, voigt_mean, voigt_sigma, voigt_width)
  
  ### double-sided crystal ball (gaussian + asymmetric exponential tails) 
  dscb1_mean   = ROOT.RooRealVar('dscb1_mean','dscb1_mean',  5.28, 5., 6.)
  dscb1_sigma  = ROOT.RooRealVar('dscb1_sigma','dscb1_sigma', 0.02, 0., 0.05) 
  dscb1_alpha  = ROOT.RooRealVar('dscb1_alpha', 'dscb1_alpha', 0.5, 0., 5.) # alpha: where the Gaussian and exponential connect, at alpha * sigma of the Gaussian (pos)
  dscb1_n      = ROOT.RooRealVar('dscb1_n', 'dscb1_n', 1., 0., 10.)         # the slope of the exponential 
  dscb2_alpha  = ROOT.RooRealVar('dscb2_alpha', 'dscb2_alpha', 0.5, 0., 5.)
  dscb2_n      = ROOT.RooRealVar('dscb2_n', 'dscb2_n', 1., 0., 10.)
  dbscb_frac12 = ROOT.RooRealVar('dbscb_frac12','dbscb_frac12', 0.0 ,1.0)
  
  dscb1        = ROOT.RooCBShape('dscb1', 'dscb1', mass, dscb1_mean, dscb1_sigma, dscb1_alpha, dscb1_n)
  dscb2        = ROOT.RooCBShape('dscb2', 'dscb2', mass, dscb1_mean, dscb1_sigma, dscb2_alpha, dscb2_n)
  dscb         = ROOT.RooAddPdf('dscb', 'dscb', ROOT.RooArgList(dscb1, dscb2), ROOT.RooArgList(dbscb_frac12)) 
  
  # choose signal model
  if whichSignalModel == 'gauss':
    fitmodel_signal = ROOT.RooAddPdf('fitmodel_signal','J/#psi K signal',ROOT.RooArgList(gauss),ROOT.RooArgList(nsig)) 
  elif whichSignalModel == 'voigt':
    fitmodel_signal = ROOT.RooAddPdf('fitmodel_signal','J/#psi K signal',ROOT.RooArgList(voigt),ROOT.RooArgList(nsig))
  elif whichSignalModel == 'dscb':
    fitmodel_signal = ROOT.RooAddPdf('fitmodel_signal','J/#psi K signal',ROOT.RooArgList(dscb),ROOT.RooArgList(nsig))
  elif whichSignalModel == 'gauss2Sum':
    fitmodel_signal = ROOT.RooAddPdf('fitmodel_signal','J/#psi K signal',ROOT.RooArgList(gauss2Sum),ROOT.RooArgList(nsig))
  else:
    print('Warning: fit model for signal was not defined')
    
  
  # choose bkg model for the full fit
  if doFullFit:
    if whichCombBkgModel == 'exp':
      fitmodel_bkg_comb = ROOT.RooAddPdf('fitmodel_bkg_comb', 'Combinatorial bkg', ROOT.RooArgList(exp), ROOT.RooArgList(nbkg_comb))
    elif whichCombBkgModel == 'poly1':
      fitmodel_bkg_comb = ROOT.RooAddPdf('fitmodel_bkg_comb', 'Combinatorial bkg', ROOT.RooArgList(poly1), ROOT.RooArgList(nbkg_comb))
    else: 
      print('Warning: fit model for combinatorial background was not defined')
    if whichPRecBkgModel == 'argus':
      fitmodel_bkg_prec = ROOT.RooAddPdf('fitmodel_bkg_prec', 'Partially reco bkg', ROOT.RooArgList(argus), ROOT.RooArgList(nbkg_prec))
    elif whichPRecBkgModel == 'argusConvGauss':
      fitmodel_bkg_prec = ROOT.RooAddPdf('fitmodel_bkg_prec', 'Partially reco bkg', ROOT.RooArgList(argusConvGauss), ROOT.RooArgList(nbkg_prec))
    elif whichPRecBkgModel == 'argusPlusGauss':
      fitmodel_bkg_prec = ROOT.RooAddPdf('fitmodel_bkg_prec', 'Partially reco bkg', ROOT.RooArgList(argusPlusGauss), ROOT.RooArgList(nbkg_prec))
    elif whichPRecBkgModel == 'cb':
      fitmodel_bkg_prec = ROOT.RooAddPdf('fitmodel_bkg_prec', 'Partially reco bkg', ROOT.RooArgList(cb), ROOT.RooArgList(nbkg_prec))
    #elif whichPRecBkgModel == 'erf':
    #  fitmodel_bkg_prec = ROOT.RooAddPdf('fitmodel_bkg_prec', 'Partially reco bkg', ROOT.RooArgList(erf), ROOT.RooArgList(nbkg_prec))
    elif whichPRecBkgModel == 'exp':
      fitmodel_bkg_prec = ROOT.RooAddPdf('fitmodel_bkg_prec', 'Partially reco bkg', ROOT.RooArgList(exp1), ROOT.RooArgList(nbkg_prec))
    else:
      print('Warning: fit model for partially reconstructed background was not defined')
      
    fitmodel_bkg_peak = ROOT.RooAddPdf('fitmodel_bkg_peak', 'J/#psi #pi bkg', ROOT.RooArgList(gaussJpsiPi), ROOT.RooArgList(nbkg_peak))

  
    if not doAddJpsiPi:
      fitmodel = ROOT.RooAddPdf('fitmodel',
                                'signal + bkg(comb) + bkg(prec)',
                                ROOT.RooArgList(fitmodel_signal,fitmodel_bkg_comb,fitmodel_bkg_prec),
                                ROOT.RooArgList(nsig,nbkg_comb,nbkg_prec)) 
      
    else:
      fitmodel = ROOT.RooAddPdf('fitmodel',
                                'signal + bkg(comb) + bkg(prec) + bkg(peak)',
                                ROOT.RooArgList(fitmodel_signal,fitmodel_bkg_comb,fitmodel_bkg_prec,fitmodel_bkg_peak),
                                ROOT.RooArgList(nsig,nbkg_comb,nbkg_prec,nbkg_peak)) 
      
  # choose bkg model for the MC-only fit
  if doSignalFit:
    fitmodel_bkg_comb = ROOT.RooAddPdf('fitmodel_bkg_comb', 'Combinatorial bkg', ROOT.RooArgList(exp), ROOT.RooArgList(nbkg_comb))
    fitmodel =  ROOT.RooAddPdf('fitmodel', 
                               'signal + bkg(comb)',  
                               ROOT.RooArgList(fitmodel_signal,fitmodel_bkg_comb), 
                               ROOT.RooArgList(nsig,nbkg_comb))


  #############
  # Signal fit 
  #############
  if doSignalFit: 
    mass.setRange(mass_xlow,mass_xhigh)
    results = fitmodel.fitTo(Rmc, RF.Extended(True), RF.Save(), RF.SumW2Error(ROOT.kTRUE)) 
    frame = mass.frame(RF.Title(''), RF.XErrorSize(0), RF.Name('MC'))
    Rmc.plotOn(frame, RF.MarkerSize(0.5)) # weirdly, it's important to plot Rmc first...
    fitmodel.plotOn(frame, RF.Components('fitmodel_signal'),RF.LineColor(ROOT.kRed), RF.LineStyle(ROOT.kDashed), RF.Name('signal'))
    fitmodel.plotOn(frame, RF.Components('fitmodel_bkg_comb'),RF.LineColor(ROOT.kOrange), RF.LineStyle(ROOT.kDashed), RF.Name('bkg_comb'))
    fitmodel.plotOn(frame, RF.LineColor(ROOT.kBlue)) # full model at the end !
    if doDrawParameters: 
      fitmodel.paramOn(frame, RF.ShowConstants(ROOT.kTRUE),RF.Format('NEU',RF.AutoPrecision()),RF.Layout(0.65,0.93,0.92))
      frame.getAttText().SetTextSize(0.02) #
    frame2 = mass.frame(RF.Title(' '))
    chisq,prob=getChiSquare(fitmodel,Rmc)
    frame2 = mass.frame(RF.Title(' '))
    drawPlot(frame,frame2,chisq,prob,None,None,label='MC_{}'.format(whichSignalModel))

  
  #############
  # Full fit 
  #############
  if doFullFit:
    results = fitmodel.fitTo(Rdata, RF.Extended(True), RF.Save()) 
    frame = mass.frame(RF.Title(''))
    Rdata.plotOn(frame, RF.MarkerSize(0.5), RF.XErrorSize(0), RF.Name('data')) 
    fitmodel.plotOn(frame, RF.Components('fitmodel_signal'),RF.LineColor(ROOT.kRed), RF.LineStyle(ROOT.kDashed), RF.Name('signal'))
    fitmodel.plotOn(frame, RF.Components('fitmodel_bkg_comb'),RF.LineColor(ROOT.kOrange), RF.LineStyle(ROOT.kDashed), RF.Name('bkg_comb'))
    fitmodel.plotOn(frame, RF.Components('fitmodel_bkg_prec'),RF.LineColor(ROOT.kGreen), RF.LineStyle(ROOT.kDashed), RF.Name('bkg_prec'))
    if doAddJpsiPi: 
      fitmodel.plotOn(frame, RF.Components('fitmodel_bkg_peak'), RF.LineColor(ROOT.kMagenta), RF.LineStyle(ROOT.kDashed), RF.Name('bkg_peak'))
    fitmodel.plotOn(frame, RF.LineColor(ROOT.kBlue)) # full model at the end !
    if doDrawParameters:
      fitmodel.paramOn(frame, RF.ShowConstants(ROOT.kTRUE),RF.Format('NEU',RF.AutoPrecision()),RF.Layout(0.65,0.93,0.92))
      frame.getAttText().SetTextSize(0.02) #
    leg = defaultLegend(0.50,0.60,0.86,0.86) # ROOT.TLegend(0.65,0.73,0.86,0.87)
    leg.AddEntry(frame.findObject("Rdata"),"Data", "EP")
    leg.AddEntry(frame.findObject("fitmodel_signal"),"Signal","L")
    leg.AddEntry(frame.findObject("fitmodel_bkg_comb"),"Combinatorial bkg", "L")
    leg.AddEntry(frame.findObject("fitmodel_bkg_prec"),"Partially reconstructed bkg", "L")
    
    frame2 = mass.frame(RF.Title(' '))
  
    chisq,prob=getChiSquare(fitmodel,Rdata)
  
    ###########
    # Calculate normalisation factor
    ##########
    # get the parameter and the parameters error 
    nSig = nsig.getVal()
    if doFiducial:
      nSelMC = Rmc.sumEntries()
      totalEffMC,totalEffMC_err = getTotalEffFiducial(nSelMC)
  
    else:
      nGenTot = 13787971.0 # before filter, number of Bbar events
      nSelMC = Rmc.sumEntries()  # get the entries in the tree 
      totalEffMC = nSelMC / nGenTot
      totalEffMC_err = ROOT.TMath.Sqrt( 1/nSelMC + 1/nGenTot ) * totalEffMC
      
    BR_JpsiK =    10.20E-04 #pm 0.19  our fit  #Gamma274/Gamma, https://pdglive.lbl.gov/BranchingRatio.action?desig=3&parCode=S041&home=MXXX045, 2021 
    BR_JpsiMuMu = 5.961E-02 #pm 0.033          #Gamma7/Gamma, https://pdglive.lbl.gov/BranchingRatio.action?desig=2&parCode=M070&home=MXXX025, 2021
    ###print ('nSel={}'.format(nSel))
    lumi = 0.774 # /fb 
  
    sigmaBpm = nSig / (BR_JpsiK*BR_JpsiMuMu) / totalEffMC / lumi
    statErr_data = 1./ROOT.TMath.Sqrt(nSig)  * sigmaBpm
    statErr_MC   = totalEffMC_err/totalEffMC * sigmaBpm
    sigmaBpmErr = ROOT.TMath.Sqrt(statErr_data**2+statErr_MC**2)

    print ('Extracted value of sigma (fb) = {:.4e} +/- {:.4e} (stat data) +/- {:.4e} (stat MC)'.format(sigmaBpm,statErr_data,statErr_MC))
  
  
    #######
    # Draw and save result
    #######
    fiduciallabel = '_inclusive' if not doFiducial else ('_fiducial' if not doFidOneBin else '_fiducialOneBin')
    selectionlabel = 'newSel' if doNewSelection else 'oldSel'
    if doNewSelection and opt.whichCutRemoved is not None:
      selectionlabel += '_no{}'.format(removed_selection_label) 
    label='data_sig{}_bkg{}_{}{}'.format(whichSignalModel,whichPRecBkgModel,selectionlabel,fiduciallabel) 
    drawPlot(frame,frame2,chisq,prob,sigmaBpm,sigmaBpmErr,lumi,label,leg)

    with open('fit{}.txt'.format(label), 'w') as fout:
      fout.write('Extracted value of sigma (fb) = {:.2e}\n'.format(sigmaBpm))
      fout.write('Chi2/ndof={:.3f}\n'.format(chisq))
      fout.write('Prob={:.3f}\n'.format(prob))
      fout.write('Selection:\n')
      for cut in sel.split('&&'):
        fout.write('  '+cut+'\n')


  #############
  # Sidebands fit
  ############
#  mass.setRange('low' , 5., 5.16)
#  mass.setRange('high', 5.5,6.)
#  fitmodel_bkg.fitTo(Rdata, RF.Extended(True), RF.Range('low,high'))
#  #fitmodel_bkg_comb.fitTo(Rdata, RF.Extended(True), RF.Range('high'))
#  
#  #plot
#  frame = mass.frame(RF.Title(''))
#  Rdata.plotOn(frame, RF.MarkerSize(0.5)) 
#  fitmodel_bkg.plotOn(frame, RF.Components('fitmodel_bkg_comb'),RF.LineColor(ROOT.kOrange))
#  fitmodel_bkg.plotOn(frame, RF.Components('fitmodel_bkg_prec'),RF.LineColor(ROOT.kGreen))
#  fitmodel_bkg.plotOn(frame, RF.LineColor(ROOT.kBlue)) # total model as the last one, so that pulls are correct...
#  #fitmodel_bkg_comb.plotOn(frame, RF.LineColor(ROOT.kBlue)) # total model as the last one, so that pulls are correct...
#  frame2 = mass.frame(RF.Title(' '))
#  drawPlot(frame,frame2,'bkg') 



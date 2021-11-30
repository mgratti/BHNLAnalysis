import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.StatOverflows(ROOT.kTRUE) # consider overflows for mean and rms calculation
ROOT.gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
ROOT.gROOT.ProcessLine('setTDRStyle()')
ROOT.gStyle.SetTitleXOffset(1.1);
ROOT.gStyle.SetTitleYOffset(1.45);
import sys
sys.path.append('/work/mratti/plotting/myplotting')
sys.path.append('../signalRw/')
import ratioplot as RP
from spares import getOverflowedHisto
import array

signalfile = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_Oct20_Oct20.root'
datafile   = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/data/Control_Oct20/ParkingBPH1_Run2018A/merged/flat_bparknano_Oct20_TEST.root'

treeData = ROOT.TChain('control_tree')
treeData.AddFile(datafile)
treeMC = ROOT.TChain('control_tree')
treeMC.AddFile(signalfile)

sel = 'b_cos2d > 0.9995  && dimu_mass < (3.097+0.05)  && dimu_mass > (3.097-0.05)  && b_mass > 5 && sv_prob > 0.08  && k_pt > 1.0  && abs(k_eta) < 1.6  && l2_pt > 2.  && abs(l2_eta) < 1.8  && l1_softid==1'
  
sel_sigWindow = ' && abs(b_mass - 5.279) < 0.07'

sel_sig = sel + sel_sigWindow 
sel_sig_w = '(' + sel_sig + ')*(weight_hlt)'  


bins = [0., 0.1, 0.2, 0.5, 1., 2., 3.]

hdata = ROOT.TH1F('hdata', 'hdata', 6, array.array('d', bins))
hmc_w = ROOT.TH1F('hmc_w', 'hmc_w', 6, array.array('d', bins))

treeData.Draw("sv_lxy>>hdata", sel_sig)
treeMC.Draw("sv_lxy>>hmc_w", sel_sig_w)

histo_saver = []
histo_saver.append(hdata)
histo_saver.append(hmc_w)


RP.makeRatioPlot( hNum=hdata, 
                  hDen=hmc_w, 
                  hDen2="", 
                  nameNum="data         (|m_{J/#psi K} - 5.279 | < 0.07)", 
                  nameDen="signal MC (|m_{J/#psi K} - 5.279 | < 0.07)", 
                  nameDen2="", 
                  xtitle='L_{xy} (cm)',
                  ytitle="a.u.", 
                  ratiotitle="data/MC", 
                  norm=True, 
                  log=False, 
                  plotName="ratio_Lxy", 
                  outDir='./',
                  ratioyrange=(0.,2.0) )


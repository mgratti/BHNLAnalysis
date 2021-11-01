import ROOT
ROOT.gROOT.SetBatch(True)
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

signalfile = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V15_control/mass999_ctau999/nanoFiles/merged/flat_bparknano_Oct20_Oct20.root'
datafile   = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/data/Control_Oct20/ParkingBPH1_Run2018A/merged/flat_bparknano_Oct20_TEST.root'



treeData = ROOT.TChain('control_tree')
treeData.AddFile(datafile)
treeMC = ROOT.TChain('control_tree')
treeMC.AddFile(signalfile)

sig_cuts = [0,1,5,10,15,20,25,30,35,40,50,60,70,80,90,100]
#sig_cuts = [0,1,5,10]


effs_sig = []
effs_bkg = []
significances = []

for lxysig in sig_cuts:

  #sel = 'b_cos2d > 0.995 && dimu_mass > (3.097-0.05) && dimu_mass < (3.097+0.05) && b_mass > 5  && sv_lxy > 0.1  && sv_prob > 0. && k_pt > 1.5  && abs(k_eta) < 2.  && l2_pt > 4.  && abs(l2_eta) < 2  && l1_softid==1'
  sel = 'b_cos2d > 0.9995  && dimu_mass < (3.097+0.05)  && dimu_mass > (3.097-0.05)  && b_mass > 5  && sv_lxy > 0.035  && sv_prob > 0.08  && k_pt > 1.0  && abs(k_eta) < 1.6  && l2_pt > 2.  && abs(l2_eta) < 1.8  && l1_softid==1'
  
  sel_testedCut = sel + ' && sv_lxysig > {}'.format(lxysig)

  sel_bkgWindow = ' && abs(b_mass - 5.279) > 0.15'

  sel_testedCut_sig = sel_testedCut
  sel_testedCut_bkg = sel_testedCut + sel_bkgWindow

  sel_sig = sel
  sel_bkg = sel + sel_bkgWindow


  hpass = ROOT.TH1F('hpass', 'hpass', 1, 5.,6.)
  htot = ROOT.TH1F('htot', 'htot', 1, 5.,6.)

  treeData.Draw("b_mass>>hpass", sel_testedCut_bkg)
  treeData.Draw("b_mass>>htot" , sel_bkg)
  peff = ROOT.TEfficiency(hpass,htot)
  eff_bkg,eff_bkg_errup,eff_bkg_errdn = peff.GetEfficiency(1), peff.GetEfficiencyErrorUp(1), peff.GetEfficiencyErrorLow(1)
  
  treeMC.Draw("b_mass>>hpass", sel_testedCut_sig)
  treeMC.Draw("b_mass>>htot", sel_sig)

  peff = ROOT.TEfficiency(hpass,htot)
  eff_sig,eff_sig_errup,eff_sig_errdn = peff.GetEfficiency(1), peff.GetEfficiencyErrorUp(1), peff.GetEfficiencyErrorLow(1)
  pass_sig = hpass.GetBinContent(1)

  effs_sig.append(eff_sig)
  effs_bkg.append(eff_bkg)
  #effs_sig.append((eff_sig,eff_sig_errup,eff_sig_errdn))
  #effs_bkg.append((eff_bkg,eff_bkg_errup,eff_bkg_errdn))
  pass_bkg = hpass.GetBinContent(1)

  significance = pass_sig / ROOT.TMath.Sqrt( pass_bkg ) if pass_bkg!=0 else -1
  significances.append(significance)

### print

header = '{:10} {:10} {:10}'.format('cut value', 'sig eff', 'bkg eff')
print header

for i,lxysig in enumerate(sig_cuts):
  line = '{:10} {:10.3f} {:10.3f}'.format( lxysig, effs_sig[i], effs_bkg[i])
  print line



fig,ax = plt.subplots()


ax.plot(np.array(effs_sig), 1-np.array(effs_bkg), '*')
ax.set_xlabel('Signal Efficiency')
ax.set_xlim(0,1)
ax.set_ylabel('Background Rejection')
ax.set_ylim(0,1)
#ax.set_xticks(minor=True)
ax.set_xticks(np.arange(0., 1.1, 0.1 ))
ax.set_yticks(np.arange(0., 1.1, 0.1 ))
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
#ax.set_yticks(minor=True)
xs = np.linspace(0,1,200)
ys = -1 * xs + 1
ax.plot(xs,ys)
fig.savefig('roc_lxysig.pdf')
fig.savefig('roc_lxysig.png')








'''
Script to produce figure 9 of the AN
'''

from ROOT import *
import numpy as np
import array
import sys

sys.path.append('/work/mratti/plotting/myplotting')
from spares import *

gStyle.SetOptStat(000000)
gStyle.SetOptTitle(0)
TH1.SetDefaultSumw2()
gROOT.SetBatch(True)

gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
gROOT.ProcessLine('setTDRStyle()')


start = 'std2018-like' #  'std2018-like' 'data-like' 'poisson'
# start distribution must be something that I can sample from, so gaussian or poisson, cannot take a histogram easily
target = 'data' # 'data' 'data-like'
# there is no such a limitation for target distribution
nExp=1

# RDstar profile (obtained from RDstar McM request)
rdst_pu = [1.286e-05,4.360e-05,1.258e-04,2.721e-04,4.548e-04,7.077e-04,1.074e-03,1.582e-03,2.286e-03,3.264e-03,4.607e-03,6.389e-03,8.650e-03,1.139e-02,1.456e-02,1.809e-02,2.190e-02,2.589e-02,2.987e-02,3.362e-02,3.686e-02,3.938e-02,4.100e-02,4.173e-02,4.178e-02,4.183e-02,4.189e-02,4.194e-02,4.199e-02,4.205e-02,4.210e-02,4.178e-02,4.098e-02,3.960e-02,3.761e-02,3.504e-02,3.193e-02,2.840e-02,2.458e-02,2.066e-02,1.680e-02,1.320e-02,9.997e-03,7.299e-03,5.139e-03,3.496e-03,2.305e-03,1.479e-03,9.280e-04,5.729e-04,3.498e-04,2.120e-04,1.280e-04,7.702e-05,4.618e-05,2.758e-05,1.641e-05,9.741e-06,5.783e-06,3.446e-06,2.066e-06,1.248e-06,7.594e-07,4.643e-07,2.842e-07,1.734e-07,1.051e-07,6.304e-08,3.733e-08,2.179e-08,1.251e-08,7.064e-09,3.920e-09,2.137e-09,1.144e-09,6.020e-10,3.111e-10,1.579e-10,7.880e-11,3.866e-11,1.866e-11,8.864e-12,4.148e-12,1.914e-12,8.721e-13,3.928e-13,1.753e-13,7.757e-14,3.413e-14,1.496e-14,6.545e-15,2.862e-15,1.253e-15,5.493e-16,2.412e-16,1.060e-16,4.658e-17,2.045e-17,8.949e-18,3.899e-18]
bins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99]
pu_rdst_mc = TH1F('pu_rdst_mc', 'MC, like R(D*)', 100, 0, 100)
for i,pu in enumerate(rdst_pu):
	pu_rdst_mc.SetBinContent(i+1,pu)

# Profiles from files
_file0 = TFile.Open("MC_PileUp_2018_Autumn18.root")
_file1 = TFile.Open("pileup_2018.root")

pu_std_mc =   _file0.Get("pileup")           # from Riccardo HNL analysis # 100 bins 1-100  # no norm
pu_bph_data = _file1.Get("pileup_2018total") # from sara's calculation, inherited through Anne-Mazarine  #200 bins 1-200 # norm
pu_std_mc.SetTitle('MC, std 2018')
pu_bph_data.SetTitle('Data, b-park')
print 'nbins pu_std_mc=', pu_std_mc.GetNbinsX()
print 'nbins pu_bph_data=', pu_bph_data.GetNbinsX()
#########pu_bph_data.Scale(100000000) # scale like the other fake ones

# Profile from Poisson mean 20 , from sampling
pu_bph_mc = TH1F('pu_bph_mc', 'like MC, pre-legacy b-park, Poisson #mu=20', 200, 0, 200)
mypois = TF1("mypois","TMath::Poisson(x,20)",0,200)
pu_bph_mc.FillRandom("mypois", 100000000) # 100M   (with 500M or more, the distirbutions gets capped... boh)

# Profile close  to pu_std_mc, Gaussian
pu_std_mc_gaus = TH1F('pu_std_mc_gaus', 'like MC, std 2018, Gaussian #mu=35,#sigma=11', 200, 0, 200)
#LHErootfile = 'root://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/BHNL_Bc_LHEGEN_v0_mgSecondRun/BHNL_Bc_LHEtoRoot_step0_njAll.root'
mygaus = TF1("mygaus","TMath::Gaus(x,35,11)",0,200)
pu_std_mc_gaus.FillRandom("mygaus", 100000000) # 100M   (with 500M or more, the distirbutions gets capped... boh)

# Profile close to pu_bph_data => pu_bph_data_gaus, Gaussian
pu_bph_data_gaus = TH1F('pu_bph_data_gaus', 'like Data b-park, Gaussian #mu=29,#sigma=8', 200, 0, 200)
mygaus2 = TF1("mygaus2","TMath::Gaus(x,29,8)",0,200)
pu_bph_data_gaus.FillRandom("mygaus2", 100000000) # 100M   (with 500M or more, the distirbutions gets capped... boh)


# Plot them all normalised together
c = TCanvas()

pu_bph_data.GetXaxis().SetTitle('Pile-up')
pu_bph_data.GetYaxis().SetTitle('a.u.')
pu_bph_data.GetXaxis().SetRangeUser(0,100)
pu_bph_mc.SetLineColor(kRed)
pu_rdst_mc.SetLineColor(kRed)
pu_std_mc.SetLineColor(kBlue)
pu_std_mc_gaus.SetLineColor(kCyan)
pu_bph_data.SetLineColor(kGreen+3)
pu_bph_data.SetLineWidth(1)
pu_bph_data_gaus.SetLineColor(kGreen)

for i in range(1,pu_bph_data.GetNbinsX()+1):
	print(i,pu_bph_data.GetBinContent(i),pu_bph_data_gaus.GetBinContent(i))

#pu_bph_mc.DrawNormalized('hist')
pu_bph_data.DrawNormalized('hist')
pu_rdst_mc.DrawNormalized('samehist')
pu_std_mc.DrawNormalized('samehist')
#pu_std_mc_gaus.DrawNormalized('samehist')
#pu_bph_data_gaus.DrawNormalized('samehist')

leg = TLegend()
leg = defaultLegend(x1=0.6,y1=0.7,x2=0.95,y2=0.90)
leg.AddEntry(pu_bph_data, pu_bph_data.GetTitle(), 'L')
leg.AddEntry(pu_rdst_mc, pu_rdst_mc.GetTitle(), 'L')
leg.AddEntry(pu_std_mc, pu_std_mc.GetTitle(), 'L')
leg.Draw('same')


#c.BuildLegend(0.6,0.6,0.9,0.9, "", "L")
c.SaveAs('pu_comparison.pdf')
c.SaveAs('pu_comparison.C')
c.SaveAs('pu_comparison.root')



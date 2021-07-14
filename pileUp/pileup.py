'''
Study different PU options and relative consequences in terms of PU-reweighting
'''

from ROOT import *
import numpy as np
import array

gStyle.SetOptStat(000000)
gStyle.SetOptTitle(0)
TH1.SetDefaultSumw2()
gROOT.SetBatch(True)

start = 'std2018-like' #  'std2018-like' 'data-like' 'poisson'
# start distribution must be something that I can sample from, so gaussian or poisson, cannot take a histogram easily
target = 'data' # 'data' 'data-like'
# there is no such a limitation for target distribution
nExp=1

# RDstar profile
rdst_pu = [1.286e-05,4.360e-05,1.258e-04,2.721e-04,4.548e-04,7.077e-04,1.074e-03,1.582e-03,2.286e-03,3.264e-03,4.607e-03,6.389e-03,8.650e-03,1.139e-02,1.456e-02,1.809e-02,2.190e-02,2.589e-02,2.987e-02,3.362e-02,3.686e-02,3.938e-02,4.100e-02,4.173e-02,4.178e-02,4.183e-02,4.189e-02,4.194e-02,4.199e-02,4.205e-02,4.210e-02,4.178e-02,4.098e-02,3.960e-02,3.761e-02,3.504e-02,3.193e-02,2.840e-02,2.458e-02,2.066e-02,1.680e-02,1.320e-02,9.997e-03,7.299e-03,5.139e-03,3.496e-03,2.305e-03,1.479e-03,9.280e-04,5.729e-04,3.498e-04,2.120e-04,1.280e-04,7.702e-05,4.618e-05,2.758e-05,1.641e-05,9.741e-06,5.783e-06,3.446e-06,2.066e-06,1.248e-06,7.594e-07,4.643e-07,2.842e-07,1.734e-07,1.051e-07,6.304e-08,3.733e-08,2.179e-08,1.251e-08,7.064e-09,3.920e-09,2.137e-09,1.144e-09,6.020e-10,3.111e-10,1.579e-10,7.880e-11,3.866e-11,1.866e-11,8.864e-12,4.148e-12,1.914e-12,8.721e-13,3.928e-13,1.753e-13,7.757e-14,3.413e-14,1.496e-14,6.545e-15,2.862e-15,1.253e-15,5.493e-16,2.412e-16,1.060e-16,4.658e-17,2.045e-17,8.949e-18,3.899e-18]
bins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99]
pu_rdst_mc = TH1F('pu_rdst_mc', 'MC, like R(D*), mixing', 100, 0, 100)
for i,pu in enumerate(rdst_pu):
        pu_rdst_mc.SetBinContent(i+1,pu)

# Profiles from files
_file0 = TFile.Open("MC_PileUp_2018_Autumn18.root")
_file1 = TFile.Open("pileup_2018.root")

pu_std_mc =   _file0.Get("pileup")           # from Riccardo HNL analysis # 100 bins 1-100  # no norm
pu_bph_data = _file1.Get("pileup_2018total") # from sara's calculation, inherited through Anne-Mazarine  #200 bins 1-200 # norm
pu_std_mc.SetTitle('MC, std 2018, pre-mixing')
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
mygaus = TF1("mygaus","TMath::Gaus(x,35,11)",0,200)
pu_std_mc_gaus.FillRandom("mygaus", 100000000) # 100M   (with 500M or more, the distirbutions gets capped... boh)

# Profile close to pu_bph_data => pu_bph_data_gaus, Gaussian
pu_bph_data_gaus = TH1F('pu_bph_data_gaus', 'like Data b-park, Gaussian #mu=29,#sigma=8', 200, 0, 200)
mygaus2 = TF1("mygaus2","TMath::Gaus(x,29,8)",0,200)
pu_bph_data_gaus.FillRandom("mygaus2", 100000000) # 100M   (with 500M or more, the distirbutions gets capped... boh)


# Plot them all normalised together
c = TCanvas()

pu_bph_mc.GetXaxis().SetTitle('<#mu>')
pu_bph_mc.GetYaxis().SetTitle('Normalized to unit area')
pu_bph_mc.GetXaxis().SetRangeUser(0,100)
pu_bph_mc.SetLineColor(kRed)
pu_rdst_mc.SetLineColor(kOrange)
pu_std_mc.SetLineColor(kBlue)
pu_std_mc_gaus.SetLineColor(kCyan)
pu_bph_data.SetLineColor(kGreen+3)
pu_bph_data.SetLineWidth(1)
pu_bph_data_gaus.SetLineColor(kGreen)

for i in range(1,pu_bph_data.GetNbinsX()+1):
        print(i,pu_bph_data.GetBinContent(i),pu_bph_data_gaus.GetBinContent(i))

pu_bph_mc.DrawNormalized('hist')
pu_std_mc.DrawNormalized('samehist')
pu_std_mc_gaus.DrawNormalized('samehist')
pu_bph_data.DrawNormalized('samehist')
#pu_bph_data_gaus.DrawNormalized('samehist')
pu_rdst_mc.DrawNormalized('samehist')


c.BuildLegend(0.6,0.6,0.9,0.9)
c.SaveAs('pu_comparison.pdf')
c.SaveAs('pu_comparison.C')
c.SaveAs('pu_comparison.root')


if start=='poisson':
        pu_start = pu_bph_mc
elif start=='std2018-like':
        pu_start = pu_std_mc_gaus
#elif start=='std2018':
#       pu_start = pu_std_mc
elif start=='data-like':
        pu_start = pu_bph_data_gaus
else:
        raise RuntimeError('not supported')

if target=='data':
        pu_target = pu_bph_data
elif target=='data-like':
        pu_target = pu_bph_data_gaus
else:
        raise RuntimeError('not supported')

# Mock the PU reweighting procedure
# let's assume to have generated the sample with a starting profile pu_start and want to reweight to the PU distribution to pu_target
# will have to apply to each event a weight given by  pu_target/pu_start , where the two distributions are both normalised
norm_pu_start = pu_start.Clone()
norm_pu_start.Scale(1./norm_pu_start.Integral())
print 'integral norm_pu_start=',norm_pu_start.Integral()
norm_pu_target = pu_target.Clone()
norm_pu_target.Scale(1./norm_pu_target.Integral())
print 'integral norm_pu_target=',norm_pu_target.Integral()
ratio = norm_pu_target.Clone()
ratio.Divide(norm_pu_start)

c2 = TCanvas()
ratio.SetLineColor(kBlue)
ratio.GetXaxis().SetTitle('<#mu>')
ratio.GetYaxis().SetTitle('{} / {}'.format(norm_pu_target.GetTitle(),norm_pu_start.GetTitle()))
ratio.GetXaxis().SetRangeUser(0,100)
# large weights not fixed...
ratio.Draw("histE")
ratio.GetXaxis().SetRangeUser(0,100)
c2.SetLogy()
c2.SaveAs('ratio_histogram.pdf')

errs = []
errs_cropped = []
biases = []
biases_cropped = []

for iExp in range(nExp):

        # calculate distribution of the weights in a mock sample
        c3 = TCanvas()
        # generate a poissonian distribution with 10K events, determine the weights for those

        weights = TH1F('weights','weights', 100,0,100000)
        if   start=='std2018' or start=='std2018-like' or start=='data-like':
                bins = list(np.logspace(np.log10(0.0000001),np.log10(1000),60,base=10))
        elif start=='poisson':
                bins = list(np.logspace(np.log10(0.01),np.log10(100000),60,base=10))
        else:
                raise RuntimeError()
        weights_small = TH1F('weights_small', 'weights_small', len(bins)-1, array.array('f',bins))


        Lxy = TH1F('Lxy', 'Lxy (mock)', 25,0.,50.)
        Lxy_weighted = TH1F('Lxy_weighted', 'Lxy (mock) w/ PU weighting', 25,0.,50.)
        Lxy_weightedCropped = TH1F('Lxy_weighted', 'Lxy (mock) w/ PU weighting, w/o weights>10 ', 25,0.,50.)
        #c3.SetLogy()
        #weights.GetYaxis().SetRangeUser(0.01,100)

        if start=='poisson':
                custom_sample = list(np.random.poisson(20,10000))
        elif start=='std2018-like':
                custom_sample = list(np.random.normal(loc=35,scale=11,size=10000))
        elif start=='data-like':
                custom_sample = list(np.random.normal(loc=29,scale=8,size=10000))
        else:
                raise RuntimeError()

        exp = list(np.random.exponential(scale=20., size=10000)) # exp(-ct/ctau), ctau=20cm., assuming betagamma=1

        all_weights = []

        for i,ip in enumerate(custom_sample):
                bin = max(1, min(ratio.GetNbinsX(), ratio.GetXaxis().FindBin(ip)))
                weight = ratio.GetBinContent(bin)
                all_weights.append(weight)
                weights.Fill(weight)
                weights_small.Fill(weight)
                Lxy.Fill(exp[i])
                Lxy_weighted.Fill(exp[i],weight)
                weight_cropped = min(weight,10)
                Lxy_weightedCropped.Fill(exp[i],weight_cropped)

        #weights.Draw('hist')
        c3.SetLogx()
        c3.SetLogy()

        weights_small.DrawNormalized('hist')
        weights_small.GetXaxis().SetTitle('weight')
        #weights_small.GetXaxis().SetRangeUser(0.001,100)
        weights_small.GetYaxis().SetTitle('Entries')
        weights_small.DrawNormalized('hist')

        label1 = TLatex(0.6, 0.6, 'Mean:{:.2f}\pm{:.2f}, #sigma={:.2f}\pm{:.2f}'.format(weights_small.GetMean(),weights_small.GetMeanError(),weights_small.GetStdDev(),weights_small.GetStdDevError()))
        label1.SetNDC()
        label1.SetTextSize(0.03)
        label1.Draw()

        print 'mean={:.2f}, sigma={:.2f}'.format( weights.GetMean(),       np.sqrt(weights.GetRMS()))
        print 'mean={:.2f}, sigma={:.2f}'.format( weights_small.GetMean(), np.sqrt(weights_small.GetRMS()))
        #print all_weights

        c3.SaveAs('weights.pdf')

        # compare mock distributions before and after reweighting
        c4 = TCanvas()
        Lxy.SetLineColor(kBlue)
        Lxy_weighted.SetLineColor(kRed)
        Lxy_weightedCropped.SetLineColor(kOrange)
        Lxy_weighted.GetXaxis().SetTitle('mock L_{xy} (cm)')
        Lxy_weighted.Draw("histPE")
        Lxy.Draw("histPEsame")
        Lxy_weightedCropped.Draw("histPEsame")

        err_Lxy = Double(0)
        err_Lxy_weighted = Double(0)
        err_Lxy_weightedCropped = Double(0)
        yield_Lxy = Lxy.IntegralAndError(0, Lxy.GetNbinsX()+2, err_Lxy)
        yield_Lxy_weighted = Lxy_weighted.IntegralAndError(0, Lxy_weighted.GetNbinsX()+2, err_Lxy_weighted)
        yield_Lxy_weightedCropped = Lxy_weightedCropped.IntegralAndError(0, Lxy_weightedCropped.GetNbinsX()+2, err_Lxy_weightedCropped)
        rel_err_Lxy = err_Lxy/yield_Lxy
        rel_err_Lxy_weighted = err_Lxy_weighted /yield_Lxy_weighted
        rel_err_Lxy_weightedCropped = err_Lxy_weightedCropped / yield_Lxy_weightedCropped

        errs.append(rel_err_Lxy_weighted)
        errs_cropped.append(rel_err_Lxy_weightedCropped)
        biases.append(yield_Lxy_weighted/yield_Lxy-1)
        biases_cropped.append(yield_Lxy_weightedCropped/yield_Lxy-1)

        yield_label = 'Yield ={:.0f}, +PUrw={:.0f}, +cropping={:.0f}'.format(yield_Lxy,yield_Lxy_weighted,yield_Lxy_weightedCropped)
        err_label   = 'RelErr={:.3f}, +PUrw={:.3f}, +cropping={:.3f}'.format(rel_err_Lxy,rel_err_Lxy_weighted,rel_err_Lxy_weightedCropped)
        print yield_label
        print err_label
        label2 = TLatex(0.5, 0.5, yield_label)
        label3 = TLatex(0.5, 0.42, err_label)
        label2.SetNDC()
        label2.SetTextSize(0.03)
        label3.SetNDC()
        label3.SetTextSize(0.03)
        c4.BuildLegend(0.6,0.6,0.9,0.9)
        label3.Draw()
        label2.Draw()

        c4.SaveAs('mockLxy_weighted.pdf')


print 'Avg stat over N={} simulations, err={:.3f} bias={:.3f}'.format(nExp,sum(errs)/nExp,sum(biases)/nExp)
print 'Avg stat over N={} simulations, err={:.3f} bias={:.3f} cropped'.format(nExp, sum(errs_cropped)/nExp, sum(biases_cropped)/nExp)


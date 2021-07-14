from ROOT import gROOT, TH1F, TH1, kBlack, kRed, kOrange, kBlue, TCanvas, gStyle, TLegend, TLatex, TFile, TChain, TPad, kWhite, TMath
import ROOT
import os.path
import sys
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *

def makeHistSettings(h):
  h.GetYaxis().SetNdivisions(505);
  h.GetYaxis().SetLabelSize(0.06);
  h.GetXaxis().SetLabelSize(0.06);
  h.GetYaxis().SetTitleSize(0.06);
  h.GetXaxis().SetTitleSize(0.08);
  h.GetYaxis().SetTitleOffset(1.0);
  h.GetXaxis().SetTitleOffset(1.0);
  h.GetYaxis().SetLabelOffset(0.01);
  h.GetXaxis().SetLabelOffset(0.01);
  h.GetXaxis().SetTickLength(0.04);
  h.SetMarkerSize(4)
  h.SetMarkerStyle(1)

def makeRatioSettings(ratioMC):
  #ratioMC.GetYaxis().SetRangeUser(0.75,1.25)
  ratioMC.GetYaxis().SetRangeUser(0.0,2.0)
  ratioMC.GetYaxis().SetNdivisions(504);
  #ratioMC.GetYaxis().SetNdivisions(508);
  ratioMC.GetYaxis().SetLabelSize(0.12);
  ratioMC.GetXaxis().SetLabelSize(0.15);
  ratioMC.GetYaxis().SetTitleSize(0.12);
  ratioMC.GetXaxis().SetTitleSize(0.15);
  ratioMC.GetYaxis().SetTitleOffset(0.35);
  ratioMC.GetXaxis().SetTitleOffset(1.2);
  ratioMC.GetYaxis().SetLabelOffset(0.01);
  ratioMC.GetXaxis().SetLabelOffset(0.03);
  ratioMC.GetXaxis().SetTickLength(0.04);
  #ratioMC.SetFillColor(kBlack)
  ratioMC.SetFillStyle(3004)
 # ratioMC.SetLineColor( kBlack )
  #ratioMC.SetLineWidth(4 )
  #ratioMC.SetMarkerStyle(1)
  #ratioMC.SetMarkerColor(kBlack)
  ratioMC.SetLineWidth(2)
  ratioMC.SetMarkerSize(4)
  ratioMC.SetMarkerStyle(1)

# hnum nominal, hden var1, hden2 var2
def makeRatioPlot(hNum, hDen, hDen2="", nameNum="num", nameDen="den", nameDen2="", xtitle="pt",ytitle="Entries", ratiotitle="Ratio", norm=False, log=True, plotName="ratio", outDir='out',ratioyrange=(0,2)):
  TH1.SetDefaultSumw2()

  # prepare settings of histos
  hNum.SetLineColor(kBlack)
  hNum.SetLineWidth(2)
  hNum.SetMarkerStyle(20)
  #hNum.SetMarkerStyle(1)
  hNum.SetMarkerColor(kBlack)
  hNum.GetYaxis().SetTitle(ytitle)

  hDen.SetLineColor(kOrange+1)
  hDen.SetMarkerColor(kOrange+1)
  hDen.SetLineWidth(2)
  hDen.SetMarkerStyle(21)
  #hDen.SetMarkerStyle(1)

  if nameDen2 != "":
    hDen2.SetLineColor(kBlue)
    hDen2.SetMarkerColor(kBlue)
    hDen2.SetLineWidth(2)
    hDen2.SetMarkerStyle(22)
    #hDen2.SetMarkerStyle(1)
    makeHistSettings(hDen2)#


  makeHistSettings(hNum)
  makeHistSettings(hDen)




  # prepare canva
  canvas=TCanvas(plotName, plotName, 600, 600)
  ROOT.SetOwnership(canvas, False) # Crucial to avoid crashes due to the way python deletes the objects
  canvas.cd()
  yMinP1=0.305;
  bottomMarginP1=0.005;
  pad1 = TPad('pad1','pad1',0,yMinP1,0.99,1)
  if log: pad1.SetLogy()
  pad1.SetBottomMargin(bottomMarginP1)
  pad1.SetFillColor(kWhite)
  pad1.SetTickx()
  pad1.SetTicky()
  pad2 = TPad('pad2','pad2',0.,0.01,.99,0.300)
  pad2.SetTopMargin(0.033)
  pad2.SetBottomMargin(0.40)
  pad2.SetFillColor(kWhite)
  pad2.SetGridy()
  pad1.SetNumber(1)
  pad2.SetNumber(2)
  pad1.Draw()
  pad2.Draw()

  # prepare legend
  #leg = TLegend(0.75,0.68,0.82,0.88,'')
  #leg.SetBorderSize(0)
  #leg.SetTextSize(0.05)
  leg=defaultLegend(x1=0.30,y1=0.7,x2=0.95,y2=0.92,mult=1.5)


  # Draw
  pad1.cd()
  histo_saver = []

  opt = 'histE'
  #hNumNorm = hNum.Clone()
  if norm:
    hNumNorm = hNum.DrawNormalized('histE')
    histo_saver.append(hNumNorm)
  else:
    hNum.Draw('histE')

#  hNum.SetMaximum(hNum.GetMaximum()*4)
  #leg.Draw('same')
  ymax = max(hNum.GetMaximum(),hDen.GetMaximum())
 
  if log==True:
    hNum.SetMaximum(ymax*30)
  else:
    hNum.SetMaximum(ymax*1.3)


  #hDenNorm = hDen.Clone()
  if norm:
    hDenNorm = hDen.DrawNormalized('samehistE')
    histo_saver.append(hDenNorm)
  else:
    hDen.Draw('samehistE')

  #hDenNorm2 = hDen2.Clone()
  if nameDen2 != "":
    if norm:
      hDenNorm2 = hDen2.DrawNormalized('samehistE')
      histo_saver.append(hDenNorm2)
    else:
      hDen2.Draw('samehistE')
    leg.AddEntry(hDen2, nameDen2, 'LP')

  leg.AddEntry(hDen, nameDen, 'LP')
  leg.AddEntry(hNum, nameNum, 'LP')
  leg.Draw('same')


  #print hNumNorm.Integral(), hDenNorm.Integral(), hDenNorm2.Integral()
  #print hNum.Integral(), hDen.Integral(), hDen2.Integral()

  #######################################################
  ### RATIO PAD
  #######################################################
  pad2.cd()
  hRatio = hNum.Clone()
  if nameDen2 != "": hRatio2 = hNum.Clone()

  if norm:
    hRatio = hNumNorm.Clone()
    hRatio2 = hNumNorm.Clone()
    #print 'Ratio integral before division'
    #print hRatio.Integral(), hRatio2.Integral()

    hRatio.Divide(hDenNorm)
    if nameDen2 != "": hRatio2.Divide(hDenNorm2)

    #print 'Ratio integral after division'
    #print hRatio.Integral(), hRatio2.Integral()


  else:
    hRatio.Divide(hDen)
    if nameDen2 != "": hRatio2.Divide(hDen2)



  #hRatio.SetLineColor(kRed+2)
  if nameDen2 != "": hRatio2.SetLineColor(kBlue)
  #print hRatio.Integral()
  makeRatioSettings(hRatio)
  if nameDen2 != "": makeRatioSettings(hRatio2)
  hRatio.GetYaxis().SetRangeUser(ratioyrange[0],ratioyrange[1])

  hRatio.GetXaxis().SetTitle(xtitle)
  hRatio.GetYaxis().SetTitle(ratiotitle)
  hRatio.SetTitle('')



  #hRatio.Draw('histE')
  #if nameDen2 != "": hRatio2.Draw('histEsame')
  hRatio.Draw('PE')
  if nameDen2 != "": hRatio2.Draw('PEsame')
  #hRatio.GetXaxis().SetRangeUser(200.,2000.)


  canvas.SaveAs('{d}/{name}.pdf'.format(d=outDir, name = plotName))
  canvas.SaveAs('{d}/{name}.png'.format(d=outDir, name = plotName))

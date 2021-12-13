import ROOT


doMyWs = True

if doMyWs:
  fName = "./ws/allData.root"
else: 
  fName = "/scratch/mratti/FinalFits_tutorial/data_2016/allData.root"

f = ROOT.TFile(fName)
w = f.Get("tagsDumper").Get("cms_hgg_13TeV")
w.Print()

mass = w.var("CMS_hgg_mass")
mass.Print("v")

#IntLumi = w.var("IntLumi")
#IntLumi.Print("v")

data1 = w.data("Data_13TeV_LxyGt5_OS")
data1.Print("v")

xframe = mass.frame()
data1.plotOn(xframe)
c = ROOT.TCanvas()
xframe.GetYaxis().SetTitleOffset(1.4)
xframe.Draw()
c.SaveAs("test.png")

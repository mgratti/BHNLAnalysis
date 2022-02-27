import ROOT


#doMyWs = True
#
#if doMyWs:
#  fName = "./ws/allData.root"
#else: 
#  fName = "/scratch/mratti/FinalFits_tutorial/data_2016/allData.root"

fName = "./ws/output_M3_BHNL.root"

f = ROOT.TFile(fName)
w = f.Get("tagsDumper").Get("cms_hgg_13TeV")
w.Print()

mass = w.var("CMS_hgg_mass")
mass.Print("v")

dZ = w.var("dZ")
dZ.Print("v")

#IntLumi = w.var("IntLumi")
#IntLumi.Print("v")

#data = w.data("Data_13TeV_LxyGt5_OS")
#data.Print("v")

data = w.data("process_signal_Lxy0to1_OS")

xframe = mass.frame()
dZframe = dZ.frame()
data.plotOn(xframe)
c = ROOT.TCanvas()
xframe.GetYaxis().SetTitleOffset(1.4)
xframe.Draw()
c.SaveAs("test_mass.png")
data.plotOn(dZframe)
dZframe.Draw()
c.SaveAs("test_dZ.png")


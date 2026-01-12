import ROOT
import sys


# Script to play around with the 1D spline output from this directory

ROOT.gROOT.SetBatch(True)

f = ROOT.TFile.Open(sys.argv[1], "READ")

f.ls()


h = f.Get("h_txz_ratio_0")

hneg = h.Clone("hneg")
hneg.Reset()

for i in range(1, hneg.GetNbinsX()+1):
    if hneg.GetXaxis().GetBinCenter(i) < 0:
        c = hneg.GetXaxis().GetBinCenter(i)
        cp = hneg.GetXaxis().FindBin(abs(c))
        hneg.SetBinContent(cp, h.GetBinContent(i))


hneg.SetLineColor(4)

c = ROOT.TCanvas("c", "c", 700, 500)
c.SetGrid()
h.SetStats(0)
h.GetYaxis().SetRangeUser(0.85, 1.05)
h.GetXaxis().SetRangeUser(0, 80)
h.Draw("HISTE")
hneg.Draw("HISTE Same")
c.Update()
c.SaveAs("dummy_plot.png")

h.Divide(hneg)

c2 = ROOT.TCanvas("c2", "c2", 700, 500)
c2.SetGrid()
h.GetYaxis().SetRangeUser(0.85, 1.2)
h.GetXaxis().SetRangeUser(0, 80)
h.Draw("HISTE")
c2.Update()
c2.SaveAs("dummy_plot2.png")



#ROOT.gPad.WaitPrimitive()

f.Close()



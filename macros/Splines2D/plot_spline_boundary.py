import ROOT
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot


fbin = ROOT.TFile.Open(sys.argv[1], "READ")
h = fbin.Get("hHit0")
h2D = h.Projection(1, 0)
#h2D = fbin.Get("h0")
h2D.SetDirectory(0)
#fbin.ls()
fbin.Close()


print("Get Bin Tree with Uproot")
inFile = uproot.open(sys.argv[2])
t = inFile["bin_tree"]
df = t.arrays(t.keys(), library="pd")


x1 = df.query("idx == 0")["x1"].values
x2 = df.query("idx == 0")["x2"].values
y1 = df.query("idx == 0")["y1"].values
y2 = df.query("idx == 0")["y2"].values

x1 = np.array(x1, dtype=int)
x2 = np.array(x2, dtype=int)
y1 = np.array(y1, dtype=int)
y2 = np.array(y2, dtype=int)

def get_point(h, x1, x2, y1, y2):
    xcen = 0.5*(h.GetXaxis().GetBinCenter(x1)+h.GetXaxis().GetBinCenter(x2))
    ycen = 0.5*(h.GetYaxis().GetBinCenter(y1)+h.GetYaxis().GetBinCenter(y2))
    xerr = abs(h.GetXaxis().GetBinCenter(x1)-h.GetXaxis().GetBinCenter(x2))/2. + h2D.GetXaxis().GetBinWidth(1)/2;
    yerr = abs(h.GetYaxis().GetBinCenter(y1)-h.GetYaxis().GetBinCenter(y2))/2. + h2D.GetYaxis().GetBinWidth(1)/2;
    
    return [xcen, ycen, xerr, yerr]
    #return [xcen, ycen, 0, 0]


for num in range(6):
    x1v = int(x1[num])
    x2v = int(x2[num])
    y1v = int(y1[num])
    y2v = int(y2[num])

    x1e = h2D.GetXaxis().GetBinLowEdge(x1v)
    x2e = h2D.GetXaxis().GetBinUpEdge(x2v)
    y1e = 0
    y2e = 0
    print("y1", y1v, "y2", y2v)
    if y1v > y2v:
        y1e = h2D.GetYaxis().GetBinUpEdge(y1v)
        y2e = h2D.GetYaxis().GetBinLowEdge(y2v)
    else:
        y1e = h2D.GetYaxis().GetBinLowEdge(y1v)
        y2e = h2D.GetYaxis().GetBinUpEdge(y2v)
    print("y1 edge", y1e, "y2 edge", y2e)

    p = get_point(h2D, x1v, x2v, y1v, y2v)


    plt.plot([x1e, x2e], [y1e, y1e], c="black")
    plt.plot([x1e, x2e], [y2e, y2e], c="black")
    plt.plot([x1e, x1e], [y1e, y2e], c="black")
    plt.plot([x2e, x2e], [y1e, y2e], c="black")

    plt.errorbar([p[0]], [p[1]], xerr=[p[2]], yerr=[p[3]], c="r", fmt="o")

    plt.xlabel("X")
    plt.ylabel("TXZ")

    plt.show()




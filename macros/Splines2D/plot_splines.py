import ROOT
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot


f = ROOT.TFile.Open(sys.argv[1], "READ")
f.ls()


for num in range(6):

    g = f.Get("mc_mpv_"+str(num))
    n = g.GetN()
    xarr = g.GetX()
    yarr = g.GetY()
    zarr = g.GetZ()

    plt.scatter(xarr, yarr, c="r")
    plt.plot([0, 0], [-90, 90], c="black")
    if num < 3:
       plt.plot([-200, 0], [-90, -90], c="black")
       plt.plot([-200, 0], [90, 90], c="black")
       plt.plot([-200, -200], [-90, 90], c="black")
    else:
       plt.plot([0, 200], [-90, -90], c="black")
       plt.plot([0, 200], [90, 90], c="black")
       plt.plot([200, 200], [-90, 90], c="black")

    plt.title("MC Scatter: "+str(num))
    plt.show()

f.Close()



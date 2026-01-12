import ROOT
import sys


f = ROOT.TFile.Open(sys.argv[1], "READ")

f.ls()


for num in range(6):
	fit = f.Get("exclude_"+str(num))

	print("")
	print("Fit", num)
	for p in range(4):
		print("Param["+str(p)+"]", fit.GetParameter(p))

	print("")

f.Close()



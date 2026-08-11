import sys
import ROOT
from array import array
#import matplotlib.pyplot as plt


def ITM(h, sig_down, sig_up, tol, result):

    mean = h.GetMean()
    sd = h.GetRMS()
    if abs(result[0] - mean) < tol:
        #return result
        return array('d', [mean, sd])
    
    # Report error on the mean
    err = sd / (h.Integral()**0.5)
    result[0] = mean
    result[1] = err

    median = array('d', [0.0])
    probs = array('d', [0.5])

    h.GetQuantiles(1, median, probs)

    hnew = h.Clone()
    hnew.Reset()
    for i in range(1, h.GetNbinsX()+1):
        x = h.GetBinCenter(i)
        if h.GetBinLowEdge(i) > median[0] + sig_up * sd or h.GetBinLowEdge(i) + h.GetBinWidth(i) < median[0] + sig_down * sd:
            continue
        hnew.SetBinContent(i, h.GetBinContent(i))
        hnew.SetBinError(i, h.GetBinError(i))

    ITM(hnew, sig_down, sig_up, tol, result)




if __name__ == "__main__":

    fmc = ROOT.TFile.Open(sys.argv[1], "READ")
    fdata = ROOT.TFile.Open(sys.argv[2], "READ")

    fmc.ls()
    print("")
    print("")

    PLANE_ID = [float(x.strip()) for x in input(
    "Enter TPC and Plane separated by commas [integers]: "
    ).split(",")]
    print()
    x_lims = [float(x.strip()) for x in input(
    "Enter min and max X separated by commas [cm]: "
    ).split(",")]
    print()
    y_lims = [float(x.strip()) for x in input(
    "Enter min and max Y separated by commas [cm]: "
    ).split(",")]
    print()
    z_lims = [float(x.strip()) for x in input(
    "Enter min and max Z separated by commas [cm]: "
    ).split(",")]
    print()

    IDX = int(3*PLANE_ID[0] + PLANE_ID[1])

   
    h_mc = fmc.Get(f"hHit"+str(IDX)+";1")
    h_data = fdata.Get(f"hHit"+str(IDX)+";1")

    bin_xmin = h_mc.GetAxis(0).FindBin(x_lims[0])
    bin_xmax = h_mc.GetAxis(0).FindBin(x_lims[1])

    xmin_edge = h_mc.GetAxis(0).GetBinLowEdge(bin_xmin)
    xmax_edge = h_mc.GetAxis(0).GetBinUpEdge(bin_xmax)

    bin_ymin = h_mc.GetAxis(1).FindBin(y_lims[0])
    bin_ymax = h_mc.GetAxis(1).FindBin(y_lims[1])

    ymin_edge = h_mc.GetAxis(1).GetBinLowEdge(bin_ymin)
    ymax_edge = h_mc.GetAxis(1).GetBinUpEdge(bin_ymax)
    

    bin_zmin = h_mc.GetAxis(2).FindBin(z_lims[0])
    bin_zmax = h_mc.GetAxis(2).FindBin(z_lims[1])

    zmin_edge = h_mc.GetAxis(2).GetBinLowEdge(bin_zmin)
    zmax_edge = h_mc.GetAxis(2).GetBinUpEdge(bin_zmax)


    h_mc.GetAxis(0).SetRange(bin_xmin, bin_xmax)
    h_data.GetAxis(0).SetRange(bin_xmin, bin_xmax)
    
    h_mc.GetAxis(1).SetRange(bin_ymin, bin_ymax)
    h_data.GetAxis(1).SetRange(bin_ymin, bin_ymax)

    h_mc.GetAxis(2).SetRange(bin_zmin, bin_zmax)
    h_data.GetAxis(2).SetRange(bin_zmin, bin_zmax)

    h_mc_itm = h_mc.Projection(3)
    h_data_itm = h_data.Projection(3)

    result_mc = array('d', [0.0, 0.0])
    ITM(h_mc_itm, -2.0, 1.75, 1e-4, result_mc)
    
    result_data = array('d', [0.0, 0.0])
    ITM(h_data_itm, -2.0, 1.75, 1e-4, result_data)

    h_mc_itm.SetName("h_mc_itm_"+str(h_mc.GetAxis(3).GetTitle()))
    h_data_itm.SetName("h_data_itm_"+str(h_data.GetAxis(3).GetTitle()))


    h_mc_itm.SetDirectory(0)
    h_data_itm.SetDirectory(0)

    outf = ROOT.TFile(sys.argv[3], "RECREATE")
    outf.cd()
    

    h_mc_itm.Write()
    h_data_itm.Write()

    par_xmin = ROOT.TParameter(float)("xmin", xmin_edge)
    par_xmax = ROOT.TParameter(float)("xmax", xmax_edge)
    
    par_ymin = ROOT.TParameter(float)("ymin", ymin_edge)
    par_ymax = ROOT.TParameter(float)("ymax", ymax_edge)
 
    par_zmin = ROOT.TParameter(float)("zmin", zmin_edge)
    par_zmax = ROOT.TParameter(float)("zmax", zmax_edge)
 
    par_itm_mc =  ROOT.TParameter(float)("itm_mc", result_mc[0])
    par_itm_data =  ROOT.TParameter(float)("itm_data", result_data[0])


    par_xmin.Write()
    par_xmax.Write()

   
    par_ymin.Write()
    par_ymax.Write()

    par_zmin.Write()
    par_zmax.Write()

    par_itm_mc.Write()
    par_itm_data.Write()

    outf.Close()
    fmc.Close()
    fdata.Close()
    print("done")



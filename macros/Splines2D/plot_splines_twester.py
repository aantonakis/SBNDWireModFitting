import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys


root_file = sys.argv[1]

zlabel = sys.argv[2]

use_zlim = False

try:
    use_zlim = int(sys.argv[3])

except:
    print("Won't use z limits!")


tgraph_names = [f"splines_{i}" for i in range(6)]
tgraph_titles = ['TPC 0 plane 0', 'TPC 0 plane 1', 'TPC 0 plane 2',
                 'TPC 1 plane 0', 'TPC 1 plane 1', 'TPC 1 plane 2']


for ig, gname in enumerate(tgraph_names):
    fig, ax = plt.subplots(figsize=(5,4))
    f = ROOT.TFile.Open(root_file)
    graph = f.Get(gname).Clone()

    npoints = graph.GetN()
    x = np.array([graph.GetX()[i] for i in range(npoints)])
    y = np.array([graph.GetY()[i] for i in range(npoints)])
    z = np.array([graph.GetZ()[i] for i in range(npoints)])

    xi = np.linspace(np.min(x), np.max(x), 40)  # or nxbins
    yi = np.linspace(np.min(y), np.max(y), 40)  # or nybins
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate: 'linear' gives smooth, 'nearest' is blocky, 'cubic' is smoother
    zi = griddata((x, y), z, (xi, yi), method='linear')

    
    zmin=0.75
    zmax=1.25

    zminv = np.min(zi)
    zmaxv = np.max(zi)
    print("Zmin", zminv)
    print("Zmax", zmaxv)

    if use_zlim:
        zmin = zminv
        zmax = zmaxv
    

    pcm = ax.pcolormesh(xi, yi, zi, cmap='coolwarm', rasterized=True, vmin=zmin, vmax=zmax)

    f.Close()

    # ax.set_xlabel("z (cm)")
    # ax.set_ylabel("y (cm)")
    ax.set_xlabel("x (cm)")
    ax.set_ylabel(r"$\theta_{xw}$ (cm)")
    ax.set_title(tgraph_titles[ig])
    ax.tick_params(axis='both', direction='in', top=True, right=True)
    fig.colorbar(pcm, ax=ax, label=zlabel)
    plt.tight_layout()

    # Save with high DPI
    #plt.savefig(f"plot_{gname}.pdf", dpi=300)
    plt.show()



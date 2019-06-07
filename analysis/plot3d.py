from __future__ import print_function
import gzip
import pickle
import numpy as np
from coffea import hist
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
import scipy.stats

with gzip.open("tagtensor.pkl.gz") as fin:
    htagtensor = pickle.load(fin)

# The binning for pt and msd are very coarse, just a range and the overflows
print("Summing pt bins:", htagtensor.identifiers("AK8Puppijet0_pt", overflow='over'))
msdbin = htagtensor.identifiers("AK8Puppijet0_msd")[1]
print("Summing msd bins:", msdbin)
htagtensor = htagtensor.sum("AK8Puppijet0_pt", overflow='over').project("AK8Puppijet0_msd", msdbin)
processes = htagtensor.identifiers("process")

axes = ["AK8Puppijet0_deepdoubleb", "AK8Puppijet0_deepdoublec", "AK8Puppijet0_deepdoublecvb"]
axes = [htagtensor.axis(axname) for axname in axes]
xyz = [array.flatten() for array in np.meshgrid(*(ax.centers() for ax in axes), indexing='ij')]

sigs = ["Zcc"] #["Zcc", "Hcc"]
print("Signal processes:", sigs)
bkgs = [p for p in processes if p not in sigs]
print("Background processes:", bkgs)
bkg = htagtensor.project("process", bkgs).values(overflow='all')[()]
sig = htagtensor.project("process", sigs).values(overflow='all')[()]
def integrate(array, xdir, ydir, zdir):
    return np.cumsum(np.cumsum(np.cumsum(array[::xdir,::ydir,::zdir], axis=0), axis=1), axis=2)

# previously 1/fb normalized
lumi = 40.
# the direction to integrate controls which corner is always in the pass region
# for 1,-1,-1, we have low double-b, high double-c, high double-cvb
intdir = [1, -1, -1]
splusb = lumi*integrate(sig+bkg, *intdir)
b = np.maximum(lumi*integrate(bkg, *intdir), 1.)  # best would be based on MC sample size
#signif = scipy.stats.norm.isf(scipy.stats.poisson.sf(splusb, mu=b))
signif = (splusb-b)/(np.sqrt(b)+1)
imax = np.unravel_index(np.argmax(signif), signif.shape)
cuts = [ax.edges(overflow='all')[::d][i] for ax,i,d in zip(axes, imax, intdir)]
print("Best significance: %.4f at BvL %.2f CvL %.2f CvB %.2f" % ((signif[imax], )+tuple(cuts)))
cut = (xyz[0] < cuts[0]) & (xyz[1] >= cuts[1]) & (xyz[2] >= cuts[2])

figures = []
pdfs = {}
for proc in processes:
    array = htagtensor.values()[(proc,)]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cmap_blue = LinearSegmentedColormap.from_list(name="transparentblue", colors=["#0000ff00", "b"])
    cmap_red = LinearSegmentedColormap.from_list(name="transparentred", colors=["#ff000000", "r"])

    density = array.flatten()
    density /= density.max()  # transparency according to full scale
    pdfs[proc] = density / density.sum()
    colors = cmap_blue(density)
    colors[cut] = cmap_red(density[cut])
    print("Efficiency for %s: %.3f" % (proc, density[cut].sum()/density.sum()))
    
    # rendering all 40^3 points makes it a bit slow, cut the mostly transparent ones
    rendercut = (density > density.max()/150.)
    # print("Reduced to %.1f%%" % (100.*rendercut.sum()/rendercut.size))
    ax.scatter(xyz[0][rendercut], xyz[1][rendercut], xyz[2][rendercut], c=colors[rendercut], marker='.', s=250)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_zlim(0,1)
    ax.set_xlabel("double-b")
    ax.set_ylabel("double-c")
    ax.set_zlabel("double-cvb")
    ax.set_title(proc)
    # https://stackoverflow.com/questions/29041326/3d-plot-with-matplotlib-hide-axes-but-keep-axis-labels
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) 
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0)) 
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    figures.append(fig)

# close the windows to exit
plt.show()

if False:
    fig,ax = plt.subplots()
    spos = axes[2].edges(overflow='all')
    sprof = np.max(np.max(signif, axis=0), axis=0)[::-1]
    ax.step(x=spos[:-1], y=sprof, where='post')
    ax.set_xlim(0,1)
    ax.autoscale(axis='y')
    ax.set_xlabel("Deep double-cvb value")
    ax.set_ylabel(r"$S/\sqrt{B+1}$")
    lumi = ax.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)" % lumi, fontsize=14, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
    fig.show()

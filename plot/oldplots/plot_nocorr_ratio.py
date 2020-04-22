#!/usr/bin/env python

import pandas as pd
import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib

##### Control Area ######
# PLOTTING ratio of P?/P_nocorr
plot_mocks = True
show_plot = False
#########################

colors = ['red', 'cyan','purple','orange','gold','green','indigo','black','gray']
marker = ['s','D','^','d','*']
SMALL_SIZE = 12
MEDIUM_SIZE = 15
BIGGER_SIZE = 20
plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)
matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.major.width'] = 1.5
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['ytick.major.width'] = 1.5
matplotlib.rcParams['errorbar.capsize'] = 4

compare = pd.read_csv(inis.save_pk_path)
if "n5000" in inis.tag:
    path_nc = "../output/pk_mocks_lyb_nocorr_n5000.csv"
    nocorr = pd.read_csv(path_nc)
    print("path nocorr: ",path_nc)
else:
    path_nc = "../output/pk_mocks_lyb_nocorr.csv"
    nocorr = pd.read_csv(path_nc)
    print("path nocorr: ",path_nc)
print("path with metals: \n", inis.save_pk_path)
print("xs_ovi = 10**({0:.2f}) * xs_alpha".format(np.log10(opt.xs_ovi/opt.xs_alpha)))

fig,ax = plt.subplots(3)
fig.set_size_inches(12,11*1)
ax[-1].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")

for i,zidx in enumerate(range(4,7),start=0):
    labels = [r"$P'_{\alpha\alpha}/P_{\alpha\alpha,nocorr}$",
              r"$P'_{TT}/P_{TT,nocorr}$",
              r"$P'_{T\alpha}/P_{T\alpha,nocorr}$"]
    custom_lines = [Line2D([0], [0], color=colors[0], lw=3, marker=None),
                    Line2D([0], [0], color=colors[1], lw=3, marker=None),
                    Line2D([0], [0], color=colors[2], lw=3, marker=None)]
    if plot_mocks:
        zmask = (nocorr.z == opt.zbin_centers[zidx])
        t_nocorr = nocorr[zmask]
        t_compare = compare[zmask]
        k = t_nocorr.k
        k_x = np.log10(k)

        ax[i].text(0.85, 0.1,"z={0}".format(opt.zbin_centers[zidx]),
                    ha='center', va='center', transform=ax[i].transAxes)
        ax[i].plot(k_x,np.ones_like(k_x),color='gray',ls='dotted',alpha=0.5)
        ax[i].plot(k_x,t_compare.Paa/t_nocorr.Paa,color=colors[0])
        ax[i].plot(k_x,t_compare.Ptot/t_nocorr.Ptot,color=colors[1])
        if inis.include_metal_model:
            pab_new = opt.metal_model_pab(a=1e-3,b=1e-3,k=k_x,pab=t_compare.Pab)
            ax[i].plot(k_x,pab_new/t_nocorr.Pab,color=colors[2])
        else:
            ax[i].plot(k_x,t_compare.Pab/t_nocorr.Pab,color=colors[2])

for i in range(len(ax)):
    box = ax[i].get_position()
    ax[i].set_position([box.x0, box.y0, box.width * 0.75, box.height])
    ax[i].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))


if inis.save_nocorr_ratio_fig:
    plt.savefig(inis.save_nocorr_ratio_fig_path)
if show_plot:
    plt.show()
plt.clf()

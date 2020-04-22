#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import numpy as np
import pandas as pd
import options as opt
import inis

##### Control Area ######
plot_inis = True
show_plot = False
#########################
colors = ['red', 'green','blue','purple','orange','gold','indigo','black','gray']
markers = ['.','d','^','s','*']
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
####################################################################################################
labels = [r"$P_{\alpha \alpha}$",r"$P_{\alpha \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None)]
pkdata = pd.read_csv(inis.save_pk_with_err_path)

fig,ax = plt.subplots(1)
fig.set_size_inches(9,7)
ax.set_xlabel(r"log$_{10}(k/[km^{-1}s])$")
ax.set_ylabel(r"$k P(k) / \pi$")


if plot_inis:
    k = np.unique(pkdata.k)
    k_x = np.log10(k)
    for i,zidx in enumerate([6,4,2]):
        pkmask = (pkdata.z == opt.zbin_centers[zidx])
        t = pkdata[pkmask]
        paa,ptt,pab = t.paa,t.ptt,t.pab
        err_paa,err_pab = t.err_paa, t.err_pab
        ax.set_yscale('log')
        ax.errorbar(k_x,k*paa/np.pi,yerr=err_paa*k/np.pi,color=colors[0], fmt=markers[i])
        ax.errorbar(k_x,k*pab/np.pi,yerr=err_pab*k/np.pi,color=colors[2], fmt=markers[i])
        custom_lines.append(Line2D([0], [0], color='k', lw=0, marker=markers[i]))
        labels.append("z: {}".format(opt.zbin_centers[zidx]))
        k_x *=1.01
        #ax.text(0.40, 0.10,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
        #                                          , ha='center', va='center', transform=ax.transAxes)

    #labels.append("Wilson+19")
    #custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.75, box.height])
ax.legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))


# plt.tight_layout()
if inis.save_paper_pk_v2:
    plt.savefig(inis.save_paper_pk_v2_path)
if show_plot:
    plt.show()
else:
    plt.clf()

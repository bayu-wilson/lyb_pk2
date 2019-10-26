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


if "n5000" in inis.tag:
    pk_path_list = ["../output/pk_mocks_lyb_nocorr_n5000.csv",
                    "../output/pk_mocks_lyb_wN_n5000.csv",
                    "../output/pk_mocks_lyb_wR2_n5000.csv",
                    "../output/pk_mocks_lyb_wNR_n5000.csv"]
    # labels = [r"$P_{\alpha,wN\_n5000}/P_{\alpha,nocorr\_n5000}$",
    #           r"$P_{\alpha,wR\_n5000}/P_{\alpha,nocorr\_n5000}$",
    #           r"$P_{\alpha,wNR\_n5000}/P_{\alpha,nocorr\_n5000}$"]
else:
    pk_path_list = ["../output/pk_mocks_lyb_nocorr.csv",
                    "../output/pk_mocks_lyb_wN.csv",
                    "../output/pk_mocks_lyb_wR.csv",
                    "../output/pk_mocks_lyb_wNR.csv"]
    # labels = [r"$P_{\alpha,wN}/P_{\alpha,nocorr}$",
    #           r"$P_{\alpha,wR}/P_{\alpha,nocorr}$",
    #           r"$P_{\alpha,wNR}/P_{\alpha,nocorr}$"]

nocorr,wN,wR,wNR = [pd.read_csv(i) for i in pk_path_list]


fig,ax = plt.subplots(9)
fig.set_size_inches(12,11*4.5)
ax[-1].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")

for i,zidx in enumerate(range(4,7),start=0):
    m,n = int(i+3),int(i+6)
    labels = [r"$P_{wN}/P_{nocorr}$",
              r"$P_{wR}/P_{nocorr}$",
              r"$P_{wNR}/P_{nocorr}$"]
    custom_lines = [Line2D([0], [0], color=colors[0], lw=3, marker=None),
                    Line2D([0], [0], color=colors[1], lw=3, marker=None),
                    Line2D([0], [0], color=colors[2], lw=3, marker=None)]

    if plot_mocks:
        zmask = (nocorr.z == opt.zbin_centers[zidx])
        t_nocorr,t_wN,t_wR,t_wNR = nocorr[zmask],wN[zmask],wR[zmask],wNR[zmask]
        #print(len(t_nocorr),len(t_wN))
        k = t_nocorr.k
        k_x = np.log10(k)

        ax[i].text(0.85, 0.95,r"$P_{\alpha \alpha}$"+"(z={0})".format(opt.zbin_centers[zidx]),
                    ha='center', va='center', transform=ax[i].transAxes)
        ax[i].plot(k_x,np.ones_like(k_x),color='gray',ls='dotted',alpha=0.5)
        ax[i].plot(k_x,t_wN.Paa/t_nocorr.Paa,color=colors[0])
        ax[i].plot(k_x,t_wR.Paa/t_nocorr.Paa,color=colors[1])
        ax[i].plot(k_x,t_wNR.Paa/t_nocorr.Paa,color=colors[2])

        ax[m].text(0.85, 0.95,r"$P_{TT}$"+"(z={0})".format(opt.zbin_centers[zidx]),
                    ha='center', va='center', transform=ax[m].transAxes)
        ax[m].plot(k_x,np.ones_like(k_x),color='gray',ls='dotted',alpha=0.5)
        ax[m].plot(k_x,t_wN.Ptot/t_nocorr.Ptot,color=colors[0])
        ax[m].plot(k_x,t_wR.Ptot/t_nocorr.Ptot,color=colors[1])
        ax[m].plot(k_x,t_wNR.Ptot/t_nocorr.Ptot,color=colors[2])

        ax[n].text(0.85, 0.95,r"$P_{T\alpha}$"+"(z={0})".format(opt.zbin_centers[zidx]),
                    ha='center', va='center', transform=ax[n].transAxes)
        ax[n].plot(k_x,np.ones_like(k_x),color='gray',ls='dotted',alpha=0.5)
        ax[n].plot(k_x,t_wN.Pab/t_nocorr.Pab,color=colors[0])
        ax[n].plot(k_x,t_wR.Pab/t_nocorr.Pab,color=colors[1])
        ax[n].plot(k_x,t_wNR.Pab/t_nocorr.Pab,color=colors[2])


for i in range(len(ax)):
    box = ax[i].get_position()
    ax[i].set_position([box.x0, box.y0, box.width * 0.75, box.height])
    ax[i].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))
    #box = ax[i].get_position()
    #ax[i].set_position([box.x0, box.y0, box.width * 0.75, box.height])
    #ax[i].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

if inis.save_pk_ratio_fig:
    plt.savefig(inis.save_pk_ratio_fig_path)
if show_plot:
    plt.show()
plt.clf()
        #print(t_nocorr)
        #sys.exit()


# pkdata = pd.read_csv(inis.save_pk_path)
#
# fig,ax = plt.subplots(opt.zbinlen)
# fig.set_size_inches(12,11*opt.zbinlen/1.5)
# ax[-1].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")

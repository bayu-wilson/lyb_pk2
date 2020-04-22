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
####################################################################################################
labels = [r"$P_{\alpha \alpha}$",r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
                #,r"$P_{TT}$", r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                #Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None),
                Line2D([0], [0], color=colors[3], lw=9, marker=None)]


pkdata = pd.read_csv(inis.save_pk_with_err_path)
fig,ax = plt.subplots(3,2)
fig.delaxes(ax[-1][1])
fig.set_size_inches(9,10)
ax[-1,0].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")
#[i.set_xlabel(r"log$_{10}(k/[km^{-1}s])$") for i in ax[-1]] #xlabels on bottom row
ax[1,0].set_ylabel(r"$k P(k) / \pi$")

zidx = 2
for i in range(3):
    for j in range(2):
        if zidx<7:
            pkmask = (pkdata.z == opt.zbin_centers[zidx])
            t = pkdata[pkmask]
            k,paa,ptt,pab,pbb = t.k,t.paa,t.ptt,t.pab,t.pbb
            err_paa,err_pab,err_ptt,err_pbb= t.err_paa, t.err_pab, t.err_ptt,t.err_pbb
            k_x = np.log10(k)

            ax[i,j].set_yscale('log')
            ax[i,j].errorbar(k_x,k*paa/np.pi,yerr=err_paa*k/np.pi,color=colors[0], fmt='.')
            #ax[i,j].errorbar(k_x*0.99,k*ptt/np.pi,yerr=err_ptt*k/np.pi,color=colors[1], fmt='.')
            ax[i,j].errorbar(k_x,k*pab/np.pi,yerr=err_pab*k/np.pi,color=colors[2], fmt='.')
            ax[i,j].errorbar(k_x,k*pbb/np.pi,yerr=err_pbb*k/np.pi,color=colors[3], fmt='.')
            ax[i,j].text(0.40, 0.10,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
                                                      , ha='center', va='center', transform=ax[i,j].transAxes)
            zidx += 1
        # if (i==2)|(i==1)&(j==1):
        #     pass
        # else:
        #     ax[i,j].set_xticklabels([])
        #     ax[i,j].set_xticks([])
ax[2,0].set_ylim(8e-3,4e-1)

if plot_inis:
    labels.append("Wilson+19")
    custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
box = ax[-1,0].get_position()
ax[-1,0].set_position([box.x0, box.y0, box.width, box.height])
ax[-1,0].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()
if inis.save_paper_pk_v1:
    plt.savefig(inis.save_paper_pk_v1_path)
if show_plot:
    plt.show()
else:
    plt.clf()

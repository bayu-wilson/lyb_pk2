#!/usr/bin/env python

import pandas as pd
import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import glob

colors = ['blue', 'red','green','purple','orange','gold','indigo','black','gray']
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

##"pk_mocks_lyb_nocorr_interpolate" #switched_order_paaptt #switched_order_zidx
#pk_mocks_lyb_nocorr_dont_interpolate  #original_order_paaptt #original_order_zidx

tag1 = "pk_mocks_lyb_nocorr_original_order_kidx" #"pk_mocks_lyb_nocorr_interpolate" #switched_order_paaptt
tag2 = "pk_mocks_lyb_nocorr_switched_order_kidx" #"pk_mocks_lyb_nocorr_dont_interpolate" #switched_order_paaptt
short_tag1 = "original_order_kidx" #"interpolated" #switched_order_paaptt
short_tag2 = "switched_order_kidx" #"non-interpolated" #switched_order_paaptt

# pkdata1 = pd.read_csv("../output/pk_june26_removing_metals.csv")
# pkdata2 = pd.read_csv("../output/pk_june26_keeping_metals.csv")
pkdata1 = pd.read_csv("../output/{}.csv".format(tag1))
pkdata2 = pd.read_csv("../output/{}.csv".format(tag2))

fig,ax = plt.subplots(opt.zbinlen,2)
fig.set_size_inches(12*2,11*opt.zbinlen/1.5)
ax[-1,0].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")


for zidx in range(opt.zbinlen):#range(4,7):#(4,7):#opt.zbinlen):
    ax[zidx,0].set_ylabel(r"$k P(k) / \pi$")
    ax[zidx,0].set_yscale('log')
    labels = [r"$P_{\alpha \alpha}$",r"$P_{TT}$",r"$P_{T\alpha}$"] #[short_tag1,short_tag2]
    custom_lines = [Line2D([0], [0], color=colors[0], lw=5, marker=None),
                    Line2D([0], [0], color=colors[1], lw=5, marker=None),
                    Line2D([0], [0], color=colors[2], lw=5, marker=None)]

    pkmask = (pkdata1.z == opt.zbin_centers[zidx])
    t1 = pkdata1[pkmask]
    t2 = pkdata2[pkmask]

    k1,paa1,ptot1,pat1 = t1.k,t1.Paa,t1.Ptot,t1.Pab
    k2,paa2,ptot2,pat2 = t2.k,t2.Paa,t2.Ptot,t2.Pab

    k_x = np.log10(k1)
    ax[zidx,0].plot(k_x,k1*paa1/np.pi,color=colors[0], marker='o',lw=0)
    ax[zidx,0].plot(k_x,k1*ptot1/np.pi,color=colors[1], marker='o',lw=0)
    ax[zidx,0].plot(k_x,k1*pat1/np.pi,color=colors[2], marker='o',lw=0)
    #ax[zidx].errorbar(k_x,k*paa/np.pi,yerr=err_paa[zidx]*k/np.pi,color=colors[0], fmt='o')
    #ax[zidx].errorbar(k_x,k*ptot/np.pi,yerr=err_ptt[zidx]*k/np.pi,color=colors[1], fmt='o')
    #ax[zidx].errorbar(k_x,k*pat/np.pi,yerr=err_pab[zidx]*k/np.pi,color=colors[2], fmt='o')
    #labels.append(inis.tag)
    labels.append(short_tag1)
    custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='o'))


    ax[zidx,0].plot(k_x,k2*paa2/np.pi,color=colors[0], marker='.',lw=0)
    ax[zidx,0].plot(k_x,k2*ptot2/np.pi,color=colors[1], marker='.',lw=0)
    ax[zidx,0].plot(k_x,k2*pat2/np.pi,color=colors[2], marker='.',lw=0)
    #labels.append(inis.tag)
    labels.append(short_tag2)
    custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))

    # Ratio
    ax[zidx,1].plot(k_x,paa1/paa2,color=colors[0])#, marker='o')
    ax[zidx,1].plot(k_x,ptot1/ptot2,color=colors[1])#, marker='o')
    ax[zidx,1].plot(k_x,pat1/pat2,color=colors[2])#, marker='o')
    labels.append("{0}/{1}".format(short_tag1,short_tag2))
    custom_lines.append(Line2D([0], [0], color='k',lw=1))


    ax[zidx,0].text(0.80, 0.95,"z={0}".format(opt.zbin_centers[zidx])
                                              , ha='center', va='center', transform=ax[zidx,0].transAxes)
    box = ax[zidx,1].get_position()
    ax[zidx,1].set_position([box.x0, box.y0, box.width * 0.75, box.height])
    ax[zidx,1].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig("figures/remove.pdf")
plt.clf()
# plt.show()

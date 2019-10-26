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
plot_obs = True
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


if inis.lya_dlas_in_lybf:
    pk_path_list = ["../output/pk_REMOVING_DLAs_all.csv",
                    "../output/pk_KEEPING_DLAs.csv"]
elif inis.lyb_dlas_in_lybf:
    pk_path_list = ["../output/pk_REMOVING_DLAs_lybf.csv",
                    "../output/pk_KEEPING_DLAs.csv"]
else:
    pk_path_list = ["../output/pk_REMOVING_DLAs_lyaf.csv",
                    "../output/pk_KEEPING_DLAs.csv"]


obs_names = ["k","z","Paa","Ptt", "Pab", "errPaa", "errPtt", "errPab"]
boot_pk_err = pd.read_csv("../output/pk_errboot_obs_corrNR.txt",delim_whitespace=True, names = obs_names)

d_arr = [pd.read_csv(i) for i in pk_path_list]

fig,ax = plt.subplots(9)
fig.set_size_inches(12,11*4.5)
ax[-1].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")

for i,zidx in enumerate(range(4,7),start=0):
    m,n = int(i+3),int(i+6)
    labels = [r"$P_{removing DLAs}/P_{keeping DLAs}$",r"$N_{removing DLAs}/N_{keeping DLAs}$" ]
    custom_lines = [Line2D([0], [0], color=colors[0], lw=3, marker=None),
                    Line2D([0], [0], color=colors[0], ls='--',lw=3, marker=None)]
    if plot_obs:
        zmask = (d_arr[0].z == opt.zbin_centers[zidx])
        t_removing_dla= d_arr[0][zmask]
        t_remaining_dla= d_arr[1][zmask]
        err_subset = boot_pk_err[zmask]

        k = t_removing_dla.k
        k_x = np.log10(k)

        ax[i].text(0.85, 0.95,r"$P_{\alpha \alpha}$"+"(z={0})".format(opt.zbin_centers[zidx]),
                    ha='center', va='center', transform=ax[i].transAxes)
        ax[i].plot(k_x,t_removing_dla.Paa/t_remaining_dla.Paa,color=colors[0])
        ax[i].plot(k_x,t_removing_dla.npix_aa/t_remaining_dla.npix_aa,color=colors[0],ls='--')
        ax[i].fill_between(k_x,
                        1+err_subset.errPaa*err_subset.k/np.pi,
                        1-err_subset.errPaa*err_subset.k/np.pi, alpha = 0.5, color = 'gray')

        ax[m].text(0.85, 0.95,r"$P_{T T}$"+"(z={0})".format(opt.zbin_centers[zidx]),
                    ha='center', va='center', transform=ax[m].transAxes)
        ax[m].plot(k_x,t_removing_dla.Ptot/t_remaining_dla.Ptot,color=colors[0])
        #ax[m].plot(k_x,t_removing_dla.Ptot/t_remaining_dla.Ptot,color=colors[0])
        ax[m].fill_between(k_x,
                        1+err_subset.errPtt*err_subset.k/np.pi,
                        1-err_subset.errPtt*err_subset.k/np.pi, alpha = 0.5, color = 'gray')

        ax[n].text(0.85, 0.95,r"$P_{\alpha T}$"+"(z={0})".format(opt.zbin_centers[zidx]),
                    ha='center', va='center', transform=ax[n].transAxes)
        ax[n].plot(k_x,t_removing_dla.Pab/t_remaining_dla.Pab,color=colors[0])
        ax[n].plot(k_x,t_removing_dla.npix_ab/t_remaining_dla.npix_ab,color=colors[0],ls='--')
        ax[n].fill_between(k_x,
                        1+err_subset.errPab*err_subset.k/np.pi,
                        1-err_subset.errPab*err_subset.k/np.pi, alpha = 0.5, color = 'gray')

for i in range(len(ax)):
    box = ax[i].get_position()
    ax[i].set_position([box.x0, box.y0, box.width * 0.75, box.height])
    ax[i].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

if inis.save_ratio_dla_fig:
    plt.savefig(inis.save_ratio_dla_fig_path)

if show_plot:
    plt.show()
plt.clf()

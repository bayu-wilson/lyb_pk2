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
ylabels = [r"$P_{\alpha \alpha}$",r"$P_{TT}$",r"$P_{\alpha \beta}$"]
legend_labels = ["z=3.8","z=4.0","z=4.2"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None)]
pkdata = pd.read_csv(inis.save_pk_with_err_path)

for j,zidx in enumerate(range(4,7)):
    fig,ax = plt.subplots(nrows=3,ncols=1)
    fig.set_size_inches(12,15)
    for i in range(3):
        if i !=2:
            ax[i].set_ylim(6e-2,3e-1)
            ax[i].set_xticklabels([])
        ax[i].set_ylabel(ylabels[i])
    ax[2].set_ylim(8e-3,2e-1)
    ax[2].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")

    if plot_inis:
        for i,m in enumerate(['paa','ptt','pab']):
            k = np.unique(pkdata.k)
            k_x = np.log10(k)
            pkmask = (pkdata.z == opt.zbin_centers[zidx])
            t = pkdata[pkmask]
            p = t[m]
            err_p = t["err_"+m]
            ax[i].errorbar(k_x,k*p/np.pi,yerr=err_p*k/np.pi,color=colors[j], fmt='.')
            k_x *=1.005
    ax[1].legend([custom_lines[j]],[legend_labels[j]],fontsize=15,loc='lower center',ncol=1)
    plt.tight_layout()

    if inis.save_paper_pk_v0:
        plt.savefig(inis.save_paper_pk_path_v0[:-4] + "{}.pdf".format(opt.zbin_centers[zidx]))
    if show_plot:
        plt.show()
    else:
        plt.clf()




# fig,ax = plt.subplots(nrows=3,ncols=1)
# fig.set_size_inches(12,15)
# for i in range(3):
#     if i !=2:
#         ax[i].set_ylim(6e-2,3e-1)
#         ax[i].set_xticklabels([])
#     ax[i].set_ylabel(ylabels[i])
# ax[2].set_ylim(8e-3,2e-1)
# ax[2].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")
#
# # ax[1].set_ylabel(r"$k P(k) / \pi$")
#
#
# if plot_inis:
#     for i,m in enumerate(['paa','ptt','pab']):
#         k = np.unique(pkdata.k)
#         k_x = np.log10(k)
#         ax[i].set_yscale('log')
#         for j,zidx in enumerate(range(4,7)):
#             pkmask = (pkdata.z == opt.zbin_centers[zidx])
#             t = pkdata[pkmask]
#             p = t[m]
#             err_p = t["err_"+m]
#             ax[i].errorbar(k_x,k*p/np.pi,yerr=err_p*k/np.pi,color=colors[j], fmt='.')
#             k_x *=1.005

# box = ax[2].get_position()
# ax[2].set_position([box.x0, box.y0+ box.height * 0.1,
#                     box.width, box.height * 0.9])
# ax[1].legend(custom_lines,legend_labels,fontsize=15,loc='lower center',ncol=3)
# plt.tight_layout()
# Shrink current axis's height by 10% on the bottom
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])

# # Put a legend below current axis
# ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
#           fancybox=True, shadow=True, ncol=5)





            # paa,ptt,pab = t.paa,t.ptt,t.pab
            # err_paa,err_pab = t.err_paa, t.err_pab
            #
            # ax.errorbar(k_x,k*paa/np.pi,yerr=err_paa*k/np.pi,color=colors[0], fmt=markers[i])
            # ax.errorbar(k_x,k*pab/np.pi,yerr=err_pab*k/np.pi,color=colors[2], fmt=markers[i])
            # custom_lines.append(Line2D([0], [0], color='k', lw=0, marker=markers[i]))
            # labels.append("z: {}".format(opt.zbin_centers[zidx]))
            # k_x *=1.01

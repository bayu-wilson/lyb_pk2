#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import numpy as np
import pandas as pd
import glob
import sys
import options as opt
import inis

##### Control Area ######
# PLOTTING MF FOR LYA, total in LYB forest, and only LYB
plot_inis = True
plot_BOSS = False
plot_vid = False
show_plot = False
#########################

def fit_mf_PalDel(z):
        return np.exp(-0.0046*(1+z)**3.3)

mf_empirical = fit_mf_PalDel(opt.zbin_centers)


mf_df = pd.read_csv(inis.save_mf_with_err_path)
# print(mf_df)
#
# sys.exit()
# mf_df #test_output
mfa = mf_df.mfa # mean lyman alpha flux
mft = mf_df.mft # mean total flux in lyman beta forest (includes both lyb and lya)
mfb = mf_df.mfb # mean lyman beta flux
err_mfa = mf_df.err_mfa # mean lyman alpha flux
err_mft = mf_df.err_mft # mean total flux in lyman beta forest (includes both lyb and lya)
err_mfb = mf_df.err_mfb # mean lyman beta flux
z_mf = mf_df.z # redshift of each redshift bin

colors = ['red','blue','green']
labels = [# Measurement on each dataset
          r'$\overline{F}_{\alpha}$',r"$\overline{F}_{tot}=\overline{F}_{\alpha}(z_{low})\overline{F}_{\beta}(z)"+
            r"=\overline{F}_z \overline{F}_{\beta}$", r'$\overline{F}_{\beta}$',
          #Datasets used
          inis.tag]
lns = ["-","-.","--"]
fig,ax = plt.subplots(1)
fig.set_size_inches(15,6)
ax.set_ylabel(r"$\overline{F}$",fontsize = 15)
ax.set_xlabel("z",fontsize = 15)

if plot_inis:
    ax.errorbar(z_mf,mfa,yerr=err_mfa,color=colors[0], fmt='.',capsize=3)
    ax.errorbar(z_mf,mft,yerr=err_mft,color=colors[1], fmt='.',capsize=3)
    ax.errorbar(z_mf,mfb,yerr=err_mfb,color=colors[2], fmt='.',capsize=3)
    # ax.scatter(z_mf,mfa,s=10,color=colors[0],label=labels[0])
    # ax.scatter(z_mf,mfb,s=10,color=colors[1],label=labels[1])
    # ax.scatter(z_mf,mft,s=10,color=colors[2],label=labels[2])
if plot_vid:
    vid_meanF = pd.read_csv("../data/mocks/meanF_n5000_nocorr.txt",
                                delim_whitespace=True,names=['z','mfa','err_mfa'])
    z_vid = vid_meanF.z
    mf_vid = vid_meanF.mfa
    labels.append('lyb_nocorr_n5000 [Vid]')

    ax.plot(z_vid,mf_vid,color=colors[0],label=labels[0],ls=lns[1])
if plot_BOSS:
    z_arr = np.linspace(3.0,4.2,100)
    mfarr = fit_mf_PalDel(z_arr)
    ax.plot(z_arr,mfarr,color = colors[0],ls=lns[2])
    labels.append("BOSS '13")

custom_lines = [Line2D([0], [0], color=colors[0], lw=0, marker='o'),
                Line2D([0], [0], color=colors[1], lw=0, marker='o'),
                Line2D([0], [0], color=colors[2], lw=0, marker='o')] # labels for each dataset
                # Line2D([0], [0], color='k', lw=1,ls=lns[0],marker=None), # label for each measurement from dataset
                # Line2D([0], [0], color='k', lw=1,ls=lns[1],marker=None),
                # Line2D([0], [0], color='k', lw=1,ls=lns[2],marker=None)]

box = ax.get_position()  # Shrink current axis by 20%
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # Put a legend to the right of the current axis
ax.legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

if inis.save_mf_fig:
    plt.savefig(inis.save_mf_fig_path)
print(inis.save_mf_fig_path)
if show_plot:
    plt.show()
plt.clf()

#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import options as opt
import inis

##### Control Area ######
plot_inis = True
plot_BOSS = True
plot_becker = True
show_plot = True
fontsize = 15
#########################

mf_df = pd.read_csv(inis.save_mf_with_err_path)
tau_a = -np.log(mf_df.mfa)
tau_t = -np.log(mf_df.mft)
tau_b = -np.log(mf_df.mfb)
err_tau_a = mf_df.err_mfa/mf_df.mfa
err_tau_t = mf_df.err_mft/mf_df.mft
err_tau_b = mf_df.err_mfb/mf_df.mfb
z_mf = mf_df.z

colors = ['red','blue','green'] # 3 colors for mfa,mft,mfb
labels = [r"$\tau_{\alpha}$",
          r"$\tau_{total}$",
          r"$\tau_{\beta}$"]
# linestyles = ["-","-.","--"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9),#'o'),
                Line2D([0], [0], color=colors[1], lw=9),
                Line2D([0], [0], color=colors[2], lw=9)]
#Plotting
fig,ax = plt.subplots(1)
fig.set_size_inches(8,5)
ax.set_ylabel(r"$\tau_{eff}$",fontsize = fontsize)
ax.set_xlabel("z",fontsize = fontsize)
# ax.grid()

if plot_inis:
    ax.errorbar(z_mf,tau_b,yerr=err_tau_b,color=colors[2], fmt='.',capsize=3,markersize=4)
    ax.errorbar(z_mf,tau_t,yerr=err_tau_t,color=colors[1], fmt='.',capsize=3,markersize=4)
    ax.errorbar(z_mf,tau_a,yerr=err_tau_a,color=colors[0], fmt='.',capsize=3,markersize=4)
    #ax.errorbar(z_mf,mfb,yerr=err_mfb,color=colors[2], fmt='.',capsize=3,markersize=4)
    #ax.errorbar(z_mf,mft,yerr=err_mft,color=colors[1], fmt='.',capsize=3,markersize=4)
    #ax.errorbar(z_mf,mfa,yerr=err_mfa,color=colors[0], fmt='.',capsize=3,markersize=4)
    labels.append('XQ-100')
    custom_lines.append(Line2D([0], [0], color='k', lw=0,marker='.'))

if plot_BOSS:
    boss_t = pd.read_csv("../data/obs/BOSS/BOSS_mf_jun24.csv")
    boss_err_mfa = boss_t.mfa*0.1#np.sqrt(boss_t.err_mfa**2+boss_t.var_mfa*boss_t.mfa**2)
    boss_err_mfb = boss_t.mfb*0.1#np.sqrt(boss_t.err_mfb**2+boss_t.var_mfb*boss_t.mfb**2)
    #boss_mft = boss_t.mfb*boss_t.mfa
    ax.errorbar(boss_t.z*1.001,-np.log(boss_t.mfb),
                yerr=boss_err_mfb/boss_t.mfb,
                color=colors[2],fmt='^',capsize=3,alpha=0.5,markersize=3)
    ax.errorbar(boss_t.z*0.999,-np.log(boss_t.mfa),
                yerr=boss_err_mfa/boss_t.mfa,
                color=colors[0],fmt='^',capsize=3,alpha=0.5,markersize=3)
    labels.append('BOSS')
    custom_lines.append(Line2D([0], [0], color='k', lw=0,marker='^'))

if plot_becker:
    beck_t = pd.read_csv("../data/obs/Becker13_DR7_mf.csv")
    # ax.errorbar(beck_t.z*0.97,beck_t.mfb,yerr=boss_err_mfb,color=colors[2],
    #             fmt='s',capsize=3,alpha=0.5,markersize=3)
    beck_tau_a = -np.log(beck_t.mfa)
    ax.errorbar(beck_t.z,beck_tau_a,yerr=beck_t.err_mfa/beck_t.mfa,color='gray',
                fmt='s',capsize=3,alpha=0.5,markersize=2)
    labels.append('Becker+13')
    custom_lines.append(Line2D([0], [0], color='k', lw=0,marker='s'))



box = ax.get_position()  # Shrink current axis by 20%
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # Put a legend to the right of the current axis
ax.legend(custom_lines,labels,fontsize=fontsize,loc='center left', bbox_to_anchor=(1, 0.5))
#
if inis.save_paper_tau:
    plt.savefig(inis.save_paper_tau_path)
if show_plot:
    plt.show()
plt.clf()

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
show_plot = False
fontsize = 18
#########################

# opt.zbinlen = 5 #sept12

# mf_df = pd.read_csv(inis.save_mf_with_err_path)
mf_df = pd.read_csv(inis.save_mf_path)
mf_boot = np.loadtxt(inis.save_boot_mf_path)
errorbars = np.reshape(np.cov(mf_boot).diagonal()**0.5,(2,opt.zbinlen))

# Continuum correction
C_a_ratio = opt.continuum_correction(opt.zbin_centers)
C_t_ratio = opt.continuum_correction(opt.zbin_centers-0.2)
C_b_ratio = np.abs(np.abs(C_a_ratio)-np.abs(C_t_ratio))

mfa = mf_df.mf_a # mean lyman alpha flux #continuum corrected
mfa_uncorr = mfa / (1-C_a_ratio)

nans_mft = np.ones_like(opt.zbin_centers)
nans_mft[:2] = np.nan
mft = mf_df.mf_tot*nans_mft # mean total flux in lyman beta forest (includes both lyb and lya)
mft_uncorr = mft / (1-C_t_ratio)*nans_mft
# one data point is bad


mfb = mf_df.mf_b # mean lyman beta flux
mfb_uncorr = mfb / (1-C_b_ratio)

err_mfa = errorbars[0]#mf_df.err_mfa
err_mft = errorbars[1]#mf_df.err_mft
err_mfb = opt.find_err_mf_beta(mfb,mfa,err_mfa,mft,err_mft)

z_mf = mf_df.z

# LYMAN BETA MEAN FLUX #sept12
# mf_beta = mf_total / mf_alpha
# dFb = Fb * sqrt((dFt/Ft)^2+(dFa/Fa)^2)

# zab_centers = opt.find_za(zbin_centers) #converting lyb zbins to the equivalent, lower, lya zbins
# len_zab = len(zab_centers)
#
# bin_zab=np.ones(len_zab)*np.nan
# for i in range(len_zab):
#     for j in range(len_zab):
#         if (zab_centers[i]>zbin_edges[j])&(zab_centers[i]<zbin_edges[j+1]):
#             bin_zab[i] = (zbin_centers[j])



colors = ['red','blue','green'] # 3 colors for mfa,mft,mfb
labels = [r"$\overline{F}_{\beta}$",
          r"$\overline{F}_{total}$",
          r"$\overline{F}_{\alpha}$"]
# linestyles = ["-","-.","--"]
custom_lines = [Line2D([0], [0], color=colors[2], lw=9),#'o'),
                Line2D([0], [0], color=colors[1], lw=9),
                Line2D([0], [0], color=colors[0], lw=9)]
#Plotting
fig,ax = plt.subplots(1)
fig.set_size_inches(9,6)#(8,5)
ax.set_ylabel(r"$\overline{F}$",fontsize = 20)
ax.set_xlabel("z",fontsize = 20)
# ax.grid()

if plot_inis:
    #ax.errorbar(z_mf,mfb,yerr=err_mfb,color=colors[2], fmt='.',capsize=4,markersize=0,lw=2,capthick=1.5)
    ax.errorbar(z_mf,mfa,yerr=err_mfa,color=colors[0],capsize=4,markersize=0,lw=1.5,capthick=2)
    ax.errorbar(z_mf,mft,yerr=err_mft,color=colors[1],capsize=4,markersize=0,lw=1.5,capthick=2)
    ax.errorbar(z_mf,mfb,yerr=err_mfb,color=colors[2],capsize=4,markersize=0,lw=1.5,capthick=2)
    ax.errorbar(z_mf,mfa_uncorr,color=colors[0],markersize=0,lw=1.5,ls='dashed')
    ax.errorbar(z_mf,mft_uncorr,color=colors[1],markersize=0,lw=1.5,ls='dashed')
    ax.errorbar(z_mf,mfb_uncorr,color=colors[2],markersize=0,lw=1.5,ls='dashed')

    #labels.append('XQ-100')
    #custom_lines.append(Line2D([0], [0], color='k', lw=0,marker='.'))

if plot_BOSS:
    boss_t = pd.read_csv("../data/obs/BOSS/BOSS_mf_jun24.csv")
    boss_err_mfa = boss_t.mfa*0.1#np.sqrt(boss_t.err_mfa**2+boss_t.var_mfa*boss_t.mfa**2)
    boss_err_mfb = boss_t.mfb*0.1#np.sqrt(boss_t.err_mfb**2+boss_t.var_mfb*boss_t.mfb**2)
    #boss_mft = boss_t.mfb*boss_t.mfa
    ax.errorbar(boss_t.z*1.001,boss_t.mfb,yerr=boss_err_mfb,color=colors[2],
                fmt='^',capsize=3,alpha=1.0,markersize=4,lw=1,capthick=1)
    ax.errorbar(boss_t.z*0.999,boss_t.mfa,yerr=boss_err_mfa,color=colors[0],
                fmt='^',capsize=3,alpha=1.0,markersize=4,lw=1,capthick=1)
    #labels.append('BOSS')
    #custom_lines.append(Line2D([0], [0], color='k', lw=0,marker='^'))

if plot_becker:
    beck_t = pd.read_csv("../data/obs/Becker13_DR7_mf.csv")
    # ax.errorbar(beck_t.z*0.97,beck_t.mfb,yerr=boss_err_mfb,color=colors[2],
    #             fmt='s',capsize=3,alpha=0.5,markersize=3)
    ax.errorbar(beck_t.z,beck_t.mfa,yerr=beck_t.err_mfa,color=colors[0],
                fmt='s',capsize=3,markersize=0,lw=1,capthick=1)
    ##labels.append('Becker+13')
    #custom_lines.append(Line2D([0], [0], color='k', lw=0,marker='s'))



ax.tick_params(axis='both', which='major', labelsize=15)
# ax.tick_params(axis='both', which='minor', labelsize=8)


# box = ax.get_position()  # Shrink current axis by 20%
#ax.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # Put a legend to the right of the current axis
ax.legend(custom_lines,labels,fontsize=15,loc='upper right', ncol=3)# bbox_to_anchor=(1, 0.5))

if inis.save_paper_mf:
    plt.savefig(inis.save_paper_mf_path)
print(inis.save_paper_mf_path)
if show_plot:
    plt.show()
plt.clf()

#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import options as opt
import inis

##### Control Area ######
plot_inis = True
plot_BOSS = False
plot_becker = False
plot_models = True
show_plot = False

cont_corrected = False
fontsize = 18
#########################

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

colors = ['red','blue','green'] # 3 colors for mfa,mft,mfb
labels = [r"$\overline{F}_{\alpha}$",
          r"$\overline{F}_{total}$",
          r"$\overline{F}_{\beta}$"]
# linestyles = ["-","-.","--"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9),#'o'),
                Line2D([0], [0], color=colors[1], lw=9),
                Line2D([0], [0], color=colors[2], lw=9)]
#Plotting
fig,ax = plt.subplots(1)
fig.set_size_inches(9,6)#(8,5)
ax.set_ylabel(r"$\overline{F}$",fontsize = 20)
ax.set_xlabel("z",fontsize = 20)
# ax.grid()

if plot_inis:
    if cont_corrected:
        ax.set_title("Continuum Corrected mf with models",fontsize=15)
        ax.errorbar(z_mf,mfa,yerr=err_mfa,color=colors[0],capsize=4,markersize=0,lw=1.5,capthick=2)
        ax.errorbar(z_mf,mft,yerr=err_mft,color=colors[1],capsize=4,markersize=0,lw=1.5,capthick=2)
        ax.errorbar(z_mf,mfb,yerr=err_mfb,color=colors[2],capsize=4,markersize=0,lw=1.5,capthick=2)
    else:
        ax.set_title("Continuum not corrected mf with models",fontsize=15)
        ax.errorbar(z_mf,mfa_uncorr,yerr=err_mfa,color=colors[0],capsize=4,markersize=0,lw=1.5,capthick=2)
        ax.errorbar(z_mf,mft_uncorr,yerr=err_mft,color=colors[1],capsize=4,markersize=0,lw=1.5,capthick=2)
        ax.errorbar(z_mf,mfb_uncorr,yerr=err_mfb,color=colors[2],capsize=4,markersize=0,lw=1.5,capthick=2)

if plot_models:
    model_names = ["redshift","model_name", "F_Lya", "F_Lyb", "F_tot"]
    model_types = ["LCDM","LCDM_HOT","LCDM_COLD","LCDM_g1.3","LCDM_COLDg1.3",
                   "LCDM_HOTg1.3","LCDM_g1.0","LCDM_COLDg1.0","LCDM_HOTg1.0"][0:1]
    model_ls = ['dashed','-.','dotted']#,'green','blue','purple','fuchsia', 'cyan','deeppink']
    if cont_corrected:
        model_table = pd.read_csv("../data/models/fbar_models_corrected.txt", delim_whitespace=True,names=model_names)
    else:
        model_table = pd.read_csv("../data/models/fbar_models_uncorrected.txt", delim_whitespace=True,names=model_names)

    for i,model_type in enumerate(model_types):
        mask = model_table.model_name == model_type
        table_subset = model_table[mask]
        ax.plot(table_subset.redshift,table_subset.F_Lya,ls=model_ls[i],color='k',lw=1.0)
        ax.plot(table_subset.redshift,table_subset.F_Lyb,ls=model_ls[i],color='k',lw=1.0)
        ax.plot(table_subset.redshift,table_subset.F_tot,ls=model_ls[i],color='k',lw=1.0)
        custom_lines.append(Line2D([0], [0], ls=model_ls[i],color='k',lw=1.0))
        labels.append(model_types[i])

if cont_corrected:
    ax.set_ylim(0.30,0.90)
else:
    ax.set_ylim(0.39,0.90)






# if plot_BOSS:
#     boss_t = pd.read_csv("../data/obs/BOSS/BOSS_mf_jun24.csv")
#     boss_err_mfa = boss_t.mfa*0.1
#     boss_err_mfb = boss_t.mfb*0.1
#     ax.errorbar(boss_t.z*1.001,boss_t.mfb,yerr=boss_err_mfb,color=colors[2],
#                 fmt='^',capsize=3,alpha=1.0,markersize=4,lw=1,capthick=1)
#     ax.errorbar(boss_t.z*0.999,boss_t.mfa,yerr=boss_err_mfa,color=colors[0],
#                 fmt='^',capsize=3,alpha=1.0,markersize=4,lw=1,capthick=1)
#
# if plot_becker:
#     beck_t = pd.read_csv("../data/obs/Becker13_DR7_mf.csv")
#     ax.errorbar(beck_t.z,beck_t.mfa,yerr=beck_t.err_mfa,color=colors[0],
#                 fmt='s',capsize=3,markersize=0,lw=1,capthick=1)


ax.tick_params(axis='both', which='major', labelsize=15)
ax.legend(custom_lines,labels,fontsize=15,loc='upper left',ncol=2)# bbox_to_anchor=(1, 0.5))


if inis.save_paper_mf_v2:
    plt.savefig(inis.save_paper_mf_v2_path)
print(inis.save_paper_mf_v2_path)
if show_plot:
    plt.show()
plt.clf()

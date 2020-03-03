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
plot_models = True
show_plot = False

# cont_corrected = False
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

colors = ['red','blue','green','black','gray'] # 3 colors for mfa,mft,mfb
labels = [r"$\overline{F}_{\beta}$",r"$\overline{F}_{T}$",r"$\overline{F}_{\alpha}$", " ",
          r"$\gamma = 1.5$",r"$\gamma = 1.0$", #" ",
          "Continuum Uncorrected", "Continuum Corrected"]
linestyles = ['-','dashed','-.','dotted']
# linestyles = ["-","-.","--"]
custom_lines = [Line2D([0], [0], color=colors[2], lw=9),
                Line2D([0], [0], color=colors[1], lw=9),
                Line2D([0], [0], color=colors[0], lw=9),
                Line2D([0], [0], lw=0),
                Line2D([0], [0], color=colors[3], lw=9),
                Line2D([0], [0], color=colors[4], lw=9),
                #Line2D([0], [0], lw=0),
                Line2D([0], [0], color='k',ls=linestyles[1]),
                Line2D([0], [0], color='k',ls=linestyles[0])]
#Plotting
fig,ax = plt.subplots(1)
fig.set_size_inches(12,8)#(8,5)
ax.set_ylabel(r"$\overline{F}$",fontsize = 30)
ax.set_xlabel("z",fontsize = 30)
# ax.grid()

if plot_inis:
    ax.errorbar(z_mf,mfa,yerr=err_mfa,color=colors[0],capsize=5,markersize=0,lw=2,capthick=2.5)
    ax.errorbar(z_mf*1.001,mfa_uncorr,yerr=err_mfa,color=colors[0],capsize=5,markersize=0,
                      lw=1.5,capthick=2.5,ls=linestyles[1])#[-1][0].set_linestyle('--')

    ax.errorbar(z_mf,mft,yerr=err_mft,color=colors[1],capsize=5,markersize=0,lw=2,capthick=2.5)
    ax.errorbar(z_mf*1.000,mft_uncorr,yerr=err_mft,color=colors[1],capsize=5,markersize=0,
                      lw=1.5,capthick=2.5,ls=linestyles[1])

    # ax.errorbar(z_mf,mfb,yerr=err_mfb,color=colors[2],capsize=4,markersize=0,lw=1.5,capthick=2)
    ax.errorbar(z_mf,mfb,yerr=err_mfb,color=colors[2],capsize=5,markersize=0,lw=2,capthick=2.5,ls=linestyles[0])
    ax.errorbar(z_mf*0.999,mfb_uncorr,yerr=err_mfb,color=colors[2],capsize=5,markersize=0,
                      lw=1.5,capthick=2.5,ls=linestyles[1])
    #labels.append(r"$\overline{F}_{\beta}$ Continuum Corrected")
    # custom_lines.append(Line2D([0], [0], color=colors[0], lw=9))

    # error_plot = ax.errorbar(z_mf*0.999,mfb_uncorr,yerr=err_mfb,color=colors[2],capsize=4,markersize=0,
    #                   lw=1.5,capthick=2,ls=linestyles[1])
    # error_plot[-1][0].set_linestyle('--') #error_plot[-1][0] is the LineCollection objects of the errorbar lines
    #[-1][0].set_linestyle('--')
    #[-1][0].set_linestyle('--')


    # labels.append(r"$\overline{F}_{\beta}$")
    # custom_lines.append(Line2D([0], [0], color=colors[0], lw=9))

if plot_models:
    model_names = ["redshift","model_name", "F_Lya", "F_Lyb", "F_tot"]
    model_types = ["LCDM","LCDM_HOT","LCDM_COLD","LCDM_g1.3","LCDM_COLDg1.3",
                   "LCDM_HOTg1.3","LCDM_g1.0","LCDM_COLDg1.0","LCDM_HOTg1.0"]
    for n,i in enumerate([0,6]):
        model_type = model_types[i] #6

        model_table_corr = pd.read_csv("../data/models/fbar_models_corrected.txt", delim_whitespace=True,names=model_names)
        mask_corr = model_table_corr.model_name == model_type
        table_subset_corr = model_table_corr[mask_corr]
        ax.plot(table_subset_corr.redshift,table_subset_corr.F_Lyb,ls=linestyles[0],color=colors[3+n],lw=1.0)
        # custom_lines.append(Line2D([0], [0], ls=linestyles[0],color='k',lw=1.0))
        # labels.append(model_type + " Continuum Corrected")

        model_table_uncorr = pd.read_csv("../data/models/fbar_models_uncorrected.txt", delim_whitespace=True,names=model_names)
        mask_uncorr = model_table_uncorr.model_name == model_type
        table_subset_uncorr = model_table_uncorr[mask_uncorr]
        ax.plot(table_subset_uncorr.redshift,table_subset_uncorr.F_Lyb,ls=linestyles[1],color=colors[3+n],lw=1.0)
        # custom_lines.append(Line2D([0], [0], ls=linestyles[1],color='k',lw=1.0))
        #labels.append(model_type)


if plot_BOSS:
    boss_t = pd.read_csv("../data/obs/BOSS/BOSS_mf_jun24.csv")
    boss_err_mfa = boss_t.mfa*0.1
    boss_err_mfb = boss_t.mfb*0.1
    ax.errorbar(boss_t.z*1.001,boss_t.mfb,yerr=boss_err_mfb,color=colors[2],
                fmt='^',capsize=3,alpha=1.0,markersize=7,lw=1,capthick=1)
                #fmt='^',capsize=3,alpha=1.0,markersize=4,lw=1,capthick=1)
    ax.errorbar(boss_t.z*0.999,boss_t.mfa,yerr=boss_err_mfa,color=colors[0],
                fmt='^',capsize=3,alpha=1.0,markersize=7,lw=1,capthick=1)
                #fmt='^',capsize=3,alpha=1.0,markersize=4,lw=1,capthick=1)

if plot_becker:
    beck_t = pd.read_csv("../data/obs/Becker13_DR7_mf.csv")
    ax.errorbar(beck_t.z,beck_t.mfa,yerr=beck_t.err_mfa,color=colors[0],
                fmt='s',capsize=3,markersize=0,lw=1,capthick=1,alpha=0.8)


ax.set_ylim(0.25,1.05)#0.31,0.85)

ax.tick_params(axis='both', which='major', labelsize=20)
ax.legend(custom_lines,labels,fontsize=20,loc='lower left',ncol=2,frameon=False)# bbox_to_anchor=(1, 0.5))

if inis.save_paper_mf_v3:
    plt.savefig(inis.save_paper_mf_v3_path)
print(inis.save_paper_mf_v3_path)
if show_plot:
    plt.show()
plt.clf()



    #ax.plot(table_subset.redshift,table_subset.F_Lya,ls=model_ls[i],color='k',lw=1.0)


    #ax.plot(table_subset.redshift,table_subset.F_tot,ls=model_ls[i],color='k',lw=1.0)


    # if cont_corrected:
    #     model_table = pd.read_csv("../data/models/fbar_models_corrected.txt", delim_whitespace=True,names=model_names)
    # else:
    #     model_table = pd.read_csv("../data/models/fbar_models_uncorrected.txt", delim_whitespace=True,names=model_names)
    #
    # for i,model_type in enumerate(model_types):
    #     mask = model_table.model_name == model_type
    #     table_subset = model_table[mask]
    #     #ax.plot(table_subset.redshift,table_subset.F_Lya,ls=model_ls[i],color='k',lw=1.0)
    #     ax.plot(table_subset.redshift,table_subset.F_Lyb,ls=model_ls[i],color='k',lw=1.0)
    #     #ax.plot(table_subset.redshift,table_subset.F_tot,ls=model_ls[i],color='k',lw=1.0)
    #     custom_lines.append(Line2D([0], [0], ls=model_ls[i],color='k',lw=1.0))
    #     labels.append(model_types[i])
#
# if cont_corrected:
#     ax.set_ylim(0.30,0.90)
# else:
#     ax.set_ylim(0.39,0.90)








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

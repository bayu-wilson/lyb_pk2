#!/usr/bin/env python

import pandas as pd
import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import glob

##### Control Area ######
# PLOTTING PK
plot_standard = True
plot_bb = False
show_plot = False
plot_sims_thermal = False
plot_sims_symp = False
plot_sims_march = True
plot_sims_feb = True
#########################
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

# if inis.lya_dlas_in_lybf:
#     pkdata_dla_remaining = pd.read_csv("../output/pk_KEEPING_DLAs.csv")
#     pkdata_dla_removed = pd.read_csv("../output/pk_REMOVING_DLAs_all.csv")
# elif inis.lyb_dlas_in_lybf:
#     pkdata_dla_remaining = pd.read_csv("../output/pk_KEEPING_DLAs.csv")
#     pkdata_dla_removed = pd.read_csv("../output/pk_REMOVING_DLAs_lybf.csv")
# else:
#     pkdata_dla_remaining = pd.read_csv("../output/pk_KEEPING_DLAs.csv")
#     pkdata_dla_removed = pd.read_csv("../output/pk_REMOVING_DLAs_lyaf.csv")
pkdata_dla_remaining = pd.read_csv("../output/pk_KEEPING_DLAs.csv")
pkdata_dla_removed = pd.read_csv(inis.save_dla_path)

if inis.mock_or_obs == "obs":
    p = np.loadtxt("../output/pk_boot_obs_corrNR.csv")
    pk_bootstrap = np.reshape(p,(8,opt.kbinlen*opt.zbinlen,inis.M))
    err_paa = np.reshape(np.cov(pk_bootstrap[1]).diagonal()**0.5,
             (opt.zbinlen,opt.kbinlen))
    err_ptt = np.reshape(np.cov(pk_bootstrap[2]).diagonal()**0.5,
             (opt.zbinlen,opt.kbinlen))
    err_pab = np.reshape(np.cov(pk_bootstrap[3]).diagonal()**0.5,
             (opt.zbinlen,opt.kbinlen))
    # if inis.save_pk_with_err:
    #     pk_everything = np.column_stack((pkdata.k,pkdata.z,
    #                     pkdata.Paa, pkdata.Ptot, pkdata.Pab,
    #                     np.concatenate(err_paa),
    #                     np.concatenate(err_ptt),
    #                     np.concatenate(err_pab)))
    #     np.savetxt("../output/pk_errboot_obs_corrNR.txt",pk_everything)

fig,ax = plt.subplots(opt.zbinlen)
fig.set_size_inches(12,11*opt.zbinlen/1.5)
ax[-1].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")

for zidx in range(opt.zbinlen):#range(4,7):#(4,7):#opt.zbinlen):
    ax[zidx].set_ylabel(r"$k P(k) / \pi$")
    ax[zidx].set_yscale('log')
    labels = [r"$P_{\alpha \alpha}$",r"$P_{T\alpha}$","N(k,z)"]
    #[r"$P_{\alpha \alpha}$",r"$P_{TT}$",r"$P_{T\alpha}$"]
    custom_lines = [Line2D([0], [0], color=colors[0], lw=3, marker=None),
                    #Line2D([0], [0], color=colors[1], lw=3, marker=None),
                    Line2D([0], [0], color=colors[1], lw=3, marker=None),
                    Line2D([0], [0], color="gray", lw=3, marker=None)]
    if plot_standard:
        pkmask = (pkdata_dla_remaining.z == opt.zbin_centers[zidx])
        t = pkdata_dla_remaining[pkmask]
        k,paa,ptot,pat = t.k,t.Paa,t.Ptot,t.Pab
        k_x = np.log10(k)
        ax[zidx].plot(k_x,k*paa/np.pi,color=colors[0],ls = '-')
        ax[zidx].plot(k_x,k*pat/np.pi,color=colors[1],ls = '-')
        ax2 = ax[zidx].twinx()
        ax2.set_yscale('log')
        ax2.plot(k_x,t.npix_aa,color=colors[0],ls='-',alpha=0.3)
        ax2.plot(k_x,t.npix_ab,color=colors[1],ls='-',alpha=0.3)
        #ax[zidx].errorbar(k_x,k*paa/np.pi,yerr=0*err_paa[zidx]*k/np.pi,color=colors[0], fmt='o')
        #ax[zidx].errorbar(k_x,k*ptot/np.pi,yerr=err_ptt[zidx]*k/np.pi,color=colors[1], fmt='o')
        #ax[zidx].errorbar(k_x,k*pat/np.pi,yerr=0*err_pab[zidx]*k/np.pi,color=colors[1], fmt='o')
        #labels.append(inis.tag)
        labels.append("Keeping DLAs")
        custom_lines.append(Line2D([0], [0], color='k',lw=1,ls = '-'))

        pkmask = (pkdata_dla_removed.z == opt.zbin_centers[zidx])
        t = pkdata_dla_removed[pkmask]
        k,paa,ptot,pat = t.k,t.Paa,t.Ptot,t.Pab
        k_x = np.log10(k)
        ax[zidx].plot(k_x,k*paa/np.pi,color=colors[0],ls = 'dotted')
        ax[zidx].plot(k_x,k*pat/np.pi,color=colors[1],ls = 'dotted')
        #ax2 = ax[zidx].twinx()
        ax2.plot(k_x,t.npix_aa,color=colors[0],ls='dotted',alpha=0.3)
        ax2.plot(k_x,t.npix_ab,color=colors[1],ls='dotted',alpha=0.3)
        #ax[zidx].errorbar(k_x,k*paa/np.pi,yerr=0*err_paa[zidx]*k/np.pi,color=colors[0], fmt='^')
        #ax[zidx].errorbar(k_x,k*ptot/np.pi,yerr=err_ptt[zidx]*k/np.pi,color=colors[1], fmt='o')
        #ax[zidx].errorbar(k_x,k*pat/np.pi,yerr=0*err_pab[zidx]*k/np.pi,color=colors[1], fmt='^')
        #labels.append(inis.tag)
        labels.append("Removing DLAs")
        custom_lines.append(Line2D([0], [0], color='k',lw=1,ls = 'dotted'))
    if (plot_sims_march) & (zidx>3):
        new_idx = zidx - 4
        ref_path = glob.glob("../data/sims/pkT_LCDM_f1.0_z*.txt")
        ref_path.sort()
        t=pd.read_csv(ref_path[new_idx],skiprows=4, delim_whitespace=True,
                       names=["k","paa","pbb","pab","pzz","ptt","pta","pza","pzb","pzpb"])
        k_sim = t.k#,t.pk
        k_mask = (k_sim<=opt.kmax)#&(k_sim>=opt.kmin)
        k_sim = k_sim[k_mask]
        k_sim_x = np.log10(k_sim)

        ax[zidx].plot(k_sim_x,k_sim*t.paa[k_mask]/np.pi,marker='None', color=colors[0],ls='--')
        #ax[zidx].plot(k_sim_x,k_sim*(t.ptt)[k_mask]/np.pi,marker='None', color=colors[1],ls='--')
        ax[zidx].plot(k_sim_x,k_sim*t.pta[k_mask]/np.pi,marker='None', color=colors[1],ls='--')
        labels.append("Sherwood Sims (3/19)")
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='dashed'))


    ax[zidx].text(0.80, 0.95,"z={0}".format(opt.zbin_centers[zidx])
                                              , ha='center', va='center', transform=ax[zidx].transAxes)
    box = ax[zidx].get_position()
    ax[zidx].set_position([box.x0, box.y0, box.width * 0.75, box.height])
    ax[zidx].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0, box.width * 0.75, box.height])
    #ax[zidx].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

# plt.tight_layout()
if inis.save_dla_fig:
    plt.savefig(inis.save_dla_fig_path)
plt.clf()

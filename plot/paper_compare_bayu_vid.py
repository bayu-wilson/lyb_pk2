#!/usr/bin/env python

# FRACTIONAL DIFFERENCE BETWEEN VID'S AND MY POWER SPECTRA. (X,Y) AXES: (P_VID/P_BAYU,K)
# OUTPUT: "figures/paper_compare_vid.pdf"

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
tick_length = 7.5
tick_width = 1.5
matplotlib.rcParams['xtick.major.size'] = tick_length
matplotlib.rcParams['xtick.major.width'] = tick_width
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.major.size'] = tick_length
matplotlib.rcParams['ytick.major.width'] = tick_width
matplotlib.rcParams['errorbar.capsize'] = 4
####################################################################################################
labels = [r"$\widehat{P}_{\alpha \alpha, Vid}$",r"$\widehat{P}_{\alpha \alpha, Bayu}$",]#r"$\widehat{\cal{P}}_{\alpha \beta}$"]
                #,, r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[0], lw=9, alpha=0.2,marker=None)]
                # Line2D([0], [0], color=colors[2], lw=9, marker=None),
                # Line2D([0], [0], lw=0)]
                # Line2D([0], [0], color=colors[3], lw=9, marker=None)]


names=['z','k','paa','err_paa','idk1','idk2','idk3']
## OBS ###
bayu_data = pd.read_csv("../output/pk_obs_corrNR.csv")#pk_errboot_obs_corrNR.txt")
vid_data = pd.read_csv("../data/obs/vid_2020/Pk_xq_wRb.txt",delim_whitespace=True,names=names)

# ### MOCKS ###
# bayu_data = pd.read_csv("../output/pk_mocks_lyb_wR_n5000.csv")#pk_errboot_obs_corrNR.txt")
# vid_data = pd.read_csv("../data/obs/vid_2020/Pk_n5000_wR.txt",delim_whitespace=True,names=names)


# plt.style.use('classic')
fig,ax = plt.subplots(2,3,sharex=True, sharey=True,gridspec_kw={'hspace': 0,'wspace': 0}) #ROW COLUMN

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# fig.text(0.5, 0.04, r"$log_{10}(k/[km^{-1}s])$", ha='center', va='center',fontsize=20)
fig.text(0.02, 0.5, r"$P_{\alpha \alpha, Vid}/P_{\alpha \alpha, Bayu}$", ha='center', va='center', rotation='vertical',fontsize=25)
#r"$log_{10}(k P(k) / \pi)$", ha='center', va='center', rotation='vertical',fontsize=25)
plt.grid(False)
# plt.xlabel(r"log$_{10}(k/[km^{-1}s])$")
# plt.ylabel(r"$log_{10}(k P(k) / \pi)$")

# fig.delaxes(ax[-1][2])
fig.set_size_inches(15,9)
ax[-1,1].set_xlabel(r"$log_{10}(k/[km^{-1}s])$",fontsize=25)
#[i.set_xlabel(r"log$_{10}(k/[km^{-1}s])$") for i in ax[-1]] #xlabels on bottom row
# ax[1,1].set_ylabel(r"$k P(k) / \pi$")

zidx = 2
for i in range(2):
    for j in range(3):
        if zidx<7:
            bayu_pkmask = (bayu_data.z == opt.zbin_centers[zidx])
            t = bayu_data[bayu_pkmask]
            k,paa,ptt,pab,pbb = t.k.values,t.Paa.values,t.Ptot.values,t.Pab.values,t.Pbb.values
            #err_paa,err_pab,err_ptt,err_pbb= t.err_paa.values, t.err_pab.values, t.err_ptt.values,t.err_pbb.values

            # ### INCLUDING RESOLUTION UNCERTAINTY
            # sigma_res_aa = opt.find_res_uncertainty(k,t.z.values,paa)
            # err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
            # sigma_res_tt = opt.find_res_uncertainty(k,t.z.values,ptt)
            # err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
            # sigma_res_ab = opt.find_res_uncertainty(k,t.z.values,pab)
            # err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
            # sigma_res_bb = opt.find_res_uncertainty(k,t.z.values,pbb)
            # err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)


            # log_err_paa = [np.abs(np.log10((paa-err_paa))-np.log10(paa)),
            #                np.abs(np.log10((paa+err_paa))-np.log10(paa))] #0.434*err_paa/paa
            # log_err_ptt = [np.abs(np.log10((ptt-err_ptt))-np.log10(ptt)),
            #                np.abs(np.log10((ptt+err_ptt))-np.log10(ptt))] #0.434*err_ptt/ptt
            # log_err_pab = [np.abs(np.log10(np.abs(pab-err_pab))-np.log10(pab)),
            #                np.abs(np.log10((pab+err_pab))-np.log10(pab))] #0.434*err_pab/pab
            k_x = np.log10(k)
            if (i == 0)&(j==2):
                pass
            else:
                #ax[i,j].set_yscale('log')
                vid_pkmask = (vid_data.z == opt.zbin_centers[zidx])
                t2 = vid_data[vid_pkmask]
                #k,paa,err_paa = t2.k.values,t2.paa.values,t2.err_paa.values
                ax[i,j].plot(k_x,t2.paa.values/paa,color='k')
                ax[i,j].plot(k_x,np.ones_like(k_x),color='gray',ls='--')

                # ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.')
                ax[i,j].text(0.82, 0.10,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
                                                          , ha='center', va='center',
                                                          transform=ax[i,j].transAxes,fontsize=25)
                ax[i,j].xaxis.set_ticks_position('both')
                ax[i,j].yaxis.set_ticks_position('both')
                ax[i,j].xaxis.set_tick_params(direction='in')#, which='top')
                ax[i,j].yaxis.set_tick_params(direction='in')#, which='top')
                #ax[i,j].set_yticks(ax[i,j].get_yticks()[::1])
                #print(zidx,i,j)


                ################################################################
                ################################################################
                ################################################################


                # ### INCLUDING RESOLUTION UNCERTAINTY
                # sigma_res_aa = opt.find_res_uncertainty(k,t2.z.values,paa)
                # err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
                #
                # log_err_paa = [np.abs(np.log10((paa-err_paa))-np.log10(paa)),
                #                np.abs(np.log10((paa+err_paa))-np.log10(paa))] #0.434*err_paa/paa
                # k_x = np.log10(k)*1.01

                #ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.',alpha=0.2)
                #ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.',alpha=0.2)
                # ax[i,j].errorbar(k_x*0.99,np.log10(k*ptt/np.pi),yerr=log_err_ptt,color=colors[1], fmt='.',alpha=0.2)
                # ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi),yerr=log_err_pab,color=colors[2], fmt='.',alpha=0.2)
                ################################################################
                ################################################################
                ################################################################
                zidx += 1

# ax[1,2].set_ylim(np.log10(8e-3),np.log10(4e-1))
ax[0,0].set_xlim(-2.7,-0.9)
# ax[0,0].set_ylim(0.5,1.5)
#ax[i,j].set_ylim(1-0.01,1+0.01)
# ax[0,0].set_ylim(-2.24,-0.35)
fig.delaxes(ax[0][2])

# if plot_inis:
#     #labels.append("Wilson+19")
#     custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
# box = ax[-1,1].get_position()
# ax[-1,1].set_position([box.x0, box.y0, box.width, box.height])
# ax[-1,1].legend(custom_lines,labels,fontsize=25,loc='center left', bbox_to_anchor=(1.1, 0.5),frameon=False)
# ax[1,0].legend(custom_lines[:2],labels[:2],loc='lower left',ncol=1,frameon=False)
# ax[1,1].legend(custom_lines[2:4],labels[2:4],loc='lower left',ncol=1,frameon=False)


# for i in ax[0,2].get_xticklabels():
#     i.set_visible(True)
ax[0,2].xaxis.set_tick_params(labelbottom=True)

# ax[0,1].set_xticklabels(['-2.5','-2.0','-1.5','-1.0'])

plt.tight_layout()
if inis.save_paper_compare_vid:
    plt.savefig(inis.save_paper_compare_vid_path)
if show_plot:
    plt.show()
else:
    plt.clf()

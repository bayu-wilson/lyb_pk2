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
labels = [r"$\widehat{P}_{\alpha \alpha}/ \widehat{P}_{\alpha \alpha+OVI}$",
            r"$\widehat{P}_{TT}/\widehat{P}_{TT+OVI}$",
            r"$\widehat{P}_{\alpha \beta}/\widehat{P}_{\alpha \beta+OVI}$",
            r"$\widehat{\cal{P}}_{\beta \beta} /\widehat{\cal{P}}_{\beta \beta+OVI}$"]
            #,r"$P_{\beta \beta}$"]
                #,, r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None),
                Line2D([0], [0], color=colors[3], lw=9, marker=None),
                Line2D([0], [0], lw=0)]
                # Line2D([0], [0], color=colors[3], lw=9, marker=None)]




# pkdata = pd.read_csv("")
# pkdata = pd.read_csv(inis.save_pk_with_err_path)
# pkdata = pd.read_csv("../output/OVI_final/pk_mocks_noB_addedB.csv")
pkdata = pd.read_csv("../output/OVI_final/pk_mocks_lyb_nocorr_noOVI.csv")
pkdata_addedOVI = pd.read_csv("../output/OVI_final/pk_mocks_lyb_nocorr_addedOVI.csv")
# pkdata_addedOVI = pd.read_csv("../output/OVI_final/pk_mocks_noB_addedB_addedOVI.csv")



fig,ax = plt.subplots(1,3,sharex=True, sharey=True,gridspec_kw={'hspace': 0,'wspace': 0}) #ROW COLUMN

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# fig.text(0.5, 0.04, r"$log_{10}(k/[km^{-1}s])$", ha='center', va='center',fontsize=20)
fig.text(0.02, 0.5, r"$P/P_{OVI}$", ha='center', va='center', rotation='vertical',fontsize=25)
plt.grid(False)
# plt.xlabel(r"log$_{10}(k/[km^{-1}s])$")
# plt.ylabel(r"$log_{10}(k P(k) / \pi)$")

# fig.delaxes(ax[-1][2])
fig.set_size_inches(15,5)
ax[1].set_xlabel(r"$log_{10}(k/[km^{-1}s])$",fontsize=25)
#[i.set_xlabel(r"log$_{10}(k/[km^{-1}s])$") for i in ax[-1]] #xlabels on bottom row
# ax[1,1].set_ylabel(r"$k P(k) / \pi$")

zidx = 4
for j in range(3):
    if zidx<7:
        pkmask = (pkdata.z == opt.zbin_centers[zidx])
        t = pkdata[pkmask]
        t_ovi = pkdata_addedOVI[pkmask]
        k = t.k.values
        k_x = np.log10(k)
        #pkmask = (pkdata.z == opt.zbin_centers[zidx])
        #k_x = np.log10(k)

        ax[j].plot(k_x,t.Paa.values/t_ovi.Paa.values,color=colors[0],marker='.')
        ax[j].plot(k_x,t.Ptot.values/t_ovi.Ptot.values,color=colors[1],marker='.')
        ax[j].plot(k_x,t.Pab.values/t_ovi.Pab.values,color=colors[2],marker='.')
        ax[j].plot(k_x,t.Pbb.values/t_ovi.Pbb.values,color=colors[3],marker='.')

        fit_k = np.linspace(opt.kmin,opt.kmax,1000)
        phi = -3*np.pi/5
        A = 0.09
        fit_ovi = 1+A*np.cos(1600*fit_k+phi)
        ax[j].plot(np.log10(fit_k),fit_ovi,color='gray',ls='--')
        # t = pkdata_keepingDLA[pkmask]
        # k,paa,ptt,pab,pbb = t.k,t.paa,t.ptt,t.pab,t.pbb
        # err_paa,err_pab,err_ptt,err_pbb= t.err_paa, t.err_pab, t.err_ptt,t.err_pbb

        # k,paa_keepDLA,ptt_keepDLA,pab_keepDLA,pbb_keepDLA = (t_keepDLA.k.values,
        #             t_keepDLA.paa.values,t_keepDLA.ptt.values,
        #             t_keepDLA.pab.values,t_keepDLA.pbb.values)
        # k,paa_removeDLA,ptt_removeDLA,pab_removeDLA,pbb_removeDLA = (t_removeDLA.k.values,
        #             t_removeDLA.paa.values,t_removeDLA.ptt.values,
        #             t_removeDLA.pab.values,t_removeDLA.pbb.values)
        # ax[j].plot(k_x,paa_metals/paa,color=colors[0], marker='.')
        # ax[j].plot(k_x,ptt_metals/ptt,color=colors[1], marker='.')
        # ax[j].plot(k_x,pab_metals/pab,color=colors[2], marker='.')
        # ax[j].plot(k_x,pbb_metals/pbb,color=colors[3], marker='.',ls='--')

        ### INCLUDING RESOLUTION UNCERTAINTY
        # sigma_res_aa = opt.find_res_uncertainty(k.values,t.z.values,paa.values)
        # err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
        # sigma_res_tt = opt.find_res_uncertainty(k.values,t.z.values,ptt.values)
        # err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
        # sigma_res_ab = opt.find_res_uncertainty(k.values,t.z.values,pab.values)
        # err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
        # sigma_res_bb = opt.find_res_uncertainty(k.values,t.z.values,pbb.values)
        # err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)

        # log_err_paa = [np.abs(np.log10((paa-err_paa))-np.log10(paa)),
        #                np.abs(np.log10((paa+err_paa))-np.log10(paa))] #0.434*err_paa/paa
        # log_err_ptt = [np.abs(np.log10((ptt-err_ptt))-np.log10(ptt)),
        #                np.abs(np.log10((ptt+err_ptt))-np.log10(ptt))] #0.434*err_ptt/ptt
        # log_err_pab = [np.abs(np.log10(np.abs(pab-err_pab))-np.log10(pab)),
        #                np.abs(np.log10((pab+err_pab))-np.log10(pab))] #0.434*err_pab/pab
        # log_err_pbb = [np.abs(np.log10(np.abs(pbb-err_pbb))-np.log10(pbb)),
        #                np.abs(np.log10((pbb+err_pbb))-np.log10(pbb))] #0.434*err_pab/pab
        #k_x = np.log10(k)

        #ax[i,j].set_yscale('log')
        # ax[j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.')
        # ax[j].errorbar(k_x*0.99,np.log10(k*ptt/np.pi),yerr=log_err_ptt,color=colors[1], fmt='.')
        # ax[j].errorbar(k_x     ,np.log10(k*pab/np.pi),yerr=log_err_pab,color=colors[2], fmt='.')
        # ax[j].errorbar(k_x*1.01,np.log10(k*pbb/np.pi),yerr=log_err_pbb,color=colors[3], fmt='.')
        # ax[j].text(0.50, 0.10,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
        #                                           , ha='center', va='center',
        #                                           transform=ax[j].transAxes,fontsize=20)
        ax[j].xaxis.set_ticks_position('both')
        ax[j].yaxis.set_ticks_position('both')
        ax[j].xaxis.set_tick_params(direction='in')#, which='top')
        ax[j].yaxis.set_tick_params(direction='in')#, which='top')
        #ax[i,j].set_yticks(ax[i,j].get_yticks()[::1])
        #print(zidx,i,j)
        zidx += 1
    # if (i==2)|(i==1)&(j==1):
    #     pass
    # else:
    #     ax[i,j].set_xticklabels([])
    #     ax[i,j].set_xticks([])
# ax[1,2].set_ylim(np.log10(8e-3),np.log10(4e-1))
# ax[0].set_xlim(-2.7,-0.9)
ax[0].set_ylim(0.7,1.3)

if plot_inis:
    #labels.append("Wilson+19")
    custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
# box = ax[1].get_position()
# ax[1].set_position([box.x0, box.y0, box.width, box.height*0.9])
# ax[1].legend(custom_lines,labels,fontsize=25,loc='center left', bbox_to_anchor=(1.1, 0.5),frameon=False)
# ax[1].legend(custom_lines,labels,ncol=4,bbox_to_anchor= (0.3, 1.01), frameon=False)
# fig.legend(custom_lines,labels,loc="lower center",borderaxespad=0.7,ncol=4,frameon=False)
ax[0].legend(custom_lines[:2],labels[:2],loc='upper left',ncol=1,frameon=False)
ax[1].legend(custom_lines[2:4],labels[2:4],loc='upper left',ncol=1,frameon=False)

plt.tight_layout()
if inis.save_paper_OVI:
    plt.savefig(inis.save_paper_OVI_path)
print(inis.save_paper_OVI_path)
if show_plot:
    plt.show()
else:
    plt.clf()



### plotting two rows
# # plt.style.use('classic')
# fig,ax = plt.subplots(2,3,sharex=True, sharey=True,gridspec_kw={'hspace': 0,'wspace': 0}) #ROW COLUMN
#
# fig.add_subplot(111, frameon=False)
# # hide tick and tick label of the big axes
# plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# # fig.text(0.5, 0.04, r"$log_{10}(k/[km^{-1}s])$", ha='center', va='center',fontsize=20)
# fig.text(0.02, 0.5, r"$P/P_{OVI}$", ha='center', va='center', rotation='vertical',fontsize=25)
# plt.grid(False)
# # plt.xlabel(r"log$_{10}(k/[km^{-1}s])$")
# # plt.ylabel(r"$log_{10}(k P(k) / \pi)$")
#
# # fig.delaxes(ax[-1][2])
# fig.set_size_inches(15,9)
# ax[-1,1].set_xlabel(r"$log_{10}(k/[km^{-1}s])$",fontsize=25)
# #[i.set_xlabel(r"log$_{10}(k/[km^{-1}s])$") for i in ax[-1]] #xlabels on bottom row
# # ax[1,1].set_ylabel(r"$k P(k) / \pi$")
#
# zidx = 2
# for i in range(2):
#     for j in range(3):
#         if zidx<7:
#             pkmask = (pkdata.z == opt.zbin_centers[zidx])
#             t = pkdata[pkmask]
#             t_addedOVI = pkdata_addedOVI[pkmask]
#             k = t.k.values
#             k_x = np.log10(k)
#             if (i == 0)&(j==2):
#                 pass
#             else:
#                 ax[i,j].plot(k_x,t.Paa.values/t_addedOVI.Paa.values,color=colors[0],marker='.')
#                 ax[i,j].plot(k_x,t.Ptot.values/t_addedOVI.Ptot.values,color=colors[1],marker='.')
#                 ax[i,j].plot(k_x,t.Pab.values/t_addedOVI.Pab.values,color=colors[2],marker='.')
#                 ax[i,j].plot(k_x,t.Pbb.values/t_addedOVI.Pbb.values,color=colors[3],marker='.')
#
#
#                 #ax[i,j].set_yscale('log')
#                 # ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.')
#                 # ax[i,j].errorbar(k_x*0.99,np.log10(k*ptt/np.pi),yerr=log_err_ptt,color=colors[1], fmt='.')
#                 # ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi),yerr=log_err_pab,color=colors[2], fmt='.')
#                 #ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi))
#                 #ax[i,j].errorbar(k_x,k*pbb/np.pi,yerr=err_pbb*k/np.pi,color=colors[3], fmt='.')
#                 ax[i,j].text(0.50, 0.9,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
#                                                           , ha='center', va='center',
#                                                           transform=ax[i,j].transAxes,fontsize=20)
#                 ax[i,j].xaxis.set_ticks_position('both')
#                 ax[i,j].yaxis.set_ticks_position('both')
#                 ax[i,j].xaxis.set_tick_params(direction='in')#, which='top')
#                 ax[i,j].yaxis.set_tick_params(direction='in')#, which='top')
#                 #ax[i,j].set_yticks(ax[i,j].get_yticks()[::1])
#                 #print(zidx,i,j)
#                 zidx += 1
# ax[1,2].set_ylim(0.5,1.5)
# ax[0,0].set_xlim(-2.7,-0.9)
# # ax[0,0].set_ylim(-2.8,-0.30)
# fig.delaxes(ax[0][2])
#
# # if plot_inis:
# #     #labels.append("Wilson+19")
# #     custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
# # box = ax[-1,1].get_position()
# # ax[-1,1].set_position([box.x0, box.y0, box.width, box.height])
# # ax[-1,1].legend(custom_lines,labels,fontsize=25,loc='center left', bbox_to_anchor=(1.1, 0.5),frameon=False)
# ax[1,0].legend(custom_lines[:2],labels[:2],loc='upper left',ncol=1,frameon=False)
# ax[1,1].legend(custom_lines[2:4],labels[2:4],loc='upper left',ncol=1,frameon=False)
#
#
# # for i in ax[0,2].get_xticklabels():
# #     i.set_visible(True)
# ax[0,2].xaxis.set_tick_params(labelbottom=True)
#
# # ax[0,1].set_xticklabels(['-2.5','-2.0','-1.5','-1.0'])
#
# plt.tight_layout()
# if inis.save_paper_OVI:
#     plt.savefig(inis.save_paper_OVI_path)
# print(inis.save_paper_OVI_path)
# if show_plot:
#     plt.show()
# else:
#     plt.clf()

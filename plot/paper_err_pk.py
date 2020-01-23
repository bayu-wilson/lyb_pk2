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
labels = [r"$\sigma_{P_{\alpha \alpha}}/P_{\alpha \alpha}$",r"$\sigma_{P_{TT}}/P_{TT}$",
          r"$\sigma_{P_{\alpha \beta}}/P_{\alpha \beta}$","data","mocks/mocks_n5000",r"$\sigma_{real}/P_{real}$",r"$\widebar{ \sigma}_{boot}/P_{real}$"]
#r"${\sigma_{\hat P_{\beta \beta}}}}$","data","mocks"]#,r"$P_{\beta \beta}$"]
                #,, r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None),
                #Line2D([0], [0], color=colors[3], lw=9, marker=None),
                Line2D([0], [0], color='k', ls='-'),
                Line2D([0], [0], color='k', ls='--'),
                Line2D([0], [0], color='k', ls='-.'),
                Line2D([0], [0], color='k', ls='dotted')]


pkdata = pd.read_csv("../output/pk_errboot_obs_corrNR.txt")
# pkdata = pd.read_csv("../output/pk_errboot_obs_uncorr.txt")

#inis.save_pk_with_err_path)
# pkdata_mocks = pd.read_csv("../output/pk_errboot_mocks_lyb_nocorr.txt")
pkdata_mocks = pd.read_csv("../output/pk_errboot_mocks_lyb_nocorr.txt")
pkdata_mocks_big = pd.read_csv("../output/bootstrap_final/pk_errboot_mocks_lyb_nocorr_n5000.txt")

mean_pk_table = np.reshape(np.loadtxt('../output/realizations/mean_pk.txt').T,(50,21,13))
var_pk_table = np.reshape(np.loadtxt('../output/realizations/var_pk.txt').T,(50,21,13))
p_xx_z = np.reshape(np.mean(mean_pk_table,axis=0),(3,7,13)) #paa,ptt,pab for 7 zbins and 13 kbins
sigma_xx_z = np.reshape(np.std(mean_pk_table,axis=0),(3,7,13)) #sigma_aa,sigma_ptt,sigma_pab for 7 zbins and 13 kbins
boot_sigma_xx_z = np.reshape(np.mean(np.sqrt(var_pk_table),axis=0),(3,7,13)) #mean of the std's
sigma_boot_sigma_xx_z = np.reshape(np.std(np.sqrt(var_pk_table),axis=0),(3,7,13)) #std of the std's


# pk_errboot_mocks_lyb_nocorr_n5000


# plt.style.use('classic')
fig,ax = plt.subplots(2,3,sharex=True, sharey=True,gridspec_kw={'hspace': 0,'wspace': 0}) #ROW COLUMN

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# fig.text(0.5, 0.04, r"$log_{10}(k/[km^{-1}s])$", ha='center', va='center',fontsize=20)
fig.text(0.02, 0.5, r"$\sigma_{P_k}/P_{k}$", ha='center', va='center', rotation='vertical',fontsize=25)
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
            pkmask = (pkdata.z == opt.zbin_centers[zidx])
            t = pkdata[pkmask]
            k,paa,ptt,pab,pbb = t.k.values,t.paa.values,t.ptt.values,t.pab.values,t.pbb.values
            err_paa,err_pab,err_ptt,err_pbb= t.err_paa.values, t.err_pab.values, t.err_ptt.values,t.err_pbb.values

            ### INCLUDING RESOLUTION UNCERTAINTY
            sigma_res_aa = opt.find_res_uncertainty(k,t.z.values,paa)
            err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
            sigma_res_tt = opt.find_res_uncertainty(k,t.z.values,ptt)
            err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
            sigma_res_ab = opt.find_res_uncertainty(k,t.z.values,pab)
            err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
            sigma_res_bb = opt.find_res_uncertainty(k,t.z.values,pbb)
            err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)
            # log_err_paa = [np.abs(np.log10((paa-err_paa))-np.log10(paa)),
            #                np.abs(np.log10((paa+err_paa))-np.log10(paa))] #0.434*err_paa/paa
            # log_err_ptt = [np.abs(np.log10((ptt-err_ptt))-np.log10(ptt)),
            #                np.abs(np.log10((ptt+err_ptt))-np.log10(ptt))] #0.434*err_ptt/ptt
            # log_err_pab = [np.abs(np.log10(np.abs(pab-err_pab))-np.log10(pab)),
            #                np.abs(np.log10((pab+err_pab))-np.log10(pab))] #0.434*err_pab/pab
            k_x = np.log10(k)
            # for test_idx in range(len(log_err_pab)):
            #     #if np.log10(k[test_idx]*pab[test_idx]/np.pi)-(log_err_pab[0])[test_idx]>0:
            #     if pab[test_idx]-err_pab[test_idx]<0:
            #         (log_err_pab[0])[test_idx] = 0
            #         print("!")

            ### MOCKS ###
            pkmask_mocks = (pkdata_mocks.z == opt.zbin_centers[zidx])
            t_mocks = pkdata_mocks[pkmask_mocks]
            k_mocks,paa_mocks,ptt_mocks = t_mocks.k.values,t_mocks.paa.values,t_mocks.ptt.values
            pab_mocks,pbb_mocks = t_mocks.pab.values,t_mocks.pbb.values
            err_paa_mocks,err_pab_mocks = t_mocks.err_paa.values, t_mocks.err_pab.values
            err_ptt_mocks,err_pbb_mocks = t_mocks.err_ptt.values,t_mocks.err_pbb.values


            ### MOCKS ###
            pkmask_mocks_big = (pkdata_mocks_big.z == opt.zbin_centers[zidx])
            t_mocks_big = pkdata_mocks_big[pkmask_mocks_big]
            k_mocks_big,paa_mocks_big,ptt_mocks_big = t_mocks_big.k.values,t_mocks_big.paa.values,t_mocks_big.ptt.values
            pab_mocks_big,pbb_mocks_big = t_mocks_big.pab.values,t_mocks_big.pbb.values
            err_paa_mocks_big,err_pab_mocks_big = t_mocks_big.err_paa.values, t_mocks_big.err_pab.values
            err_ptt_mocks_big,err_pbb_mocks_big = t_mocks_big.err_ptt.values,t_mocks_big.err_pbb.values

            if (i == 0)&(j==2):
                pass
            else:
                ax[i,j].plot(k_x,err_paa/paa,color=colors[0])
                ax[i,j].plot(k_x,err_ptt/ptt,color=colors[1])
                ax[i,j].plot(k_x,err_pab/pab,color=colors[2])
                #ax[i,j].plot(k_x,np.log10(err_pbb),color=colors[3])

                ### mocks ###
                ax[i,j].plot(k_x,err_paa_mocks/paa_mocks,color=colors[0],ls='--',alpha=0.3)
                ax[i,j].plot(k_x,err_ptt_mocks/ptt_mocks,color=colors[1],ls='--',alpha=0.3)
                ax[i,j].plot(k_x,err_pab_mocks/pab_mocks,color=colors[2],ls='--',alpha=0.3)


                ### mocks n5000 ###
                ax[i,j].plot(k_x,err_paa_mocks_big/paa_mocks_big*np.sqrt(50),color=colors[0],ls='--')
                ax[i,j].plot(k_x,err_ptt_mocks_big/ptt_mocks_big*np.sqrt(50),color=colors[1],ls='--')
                ax[i,j].plot(k_x,err_pab_mocks_big/pab_mocks_big*np.sqrt(50),color=colors[2],ls='--')



                ### Realization stuff ###
                ax[i,j].plot(k_x,boot_sigma_xx_z[2][zidx]/p_xx_z[2][zidx],color=colors[2],ls='dotted')
                ax[i,j].plot(k_x,sigma_xx_z[2][zidx]/p_xx_z[2][zidx],color=colors[2],ls='-.')
                ax[i,j].fill_between(x=k_x, y1=(boot_sigma_xx_z[2][zidx]-sigma_boot_sigma_xx_z[2][zidx])/p_xx_z[2][zidx],
                                            y2=(boot_sigma_xx_z[2][zidx]+sigma_boot_sigma_xx_z[2][zidx])/p_xx_z[2][zidx],
                                            color=colors[2],alpha=0.5)

                ax[i,j].plot(k_x,boot_sigma_xx_z[1][zidx]/p_xx_z[1][zidx],color=colors[1],ls='dotted')
                ax[i,j].plot(k_x,sigma_xx_z[1][zidx]/p_xx_z[1][zidx],color=colors[1],ls='-.')
                ax[i,j].fill_between(x=k_x, y1=(boot_sigma_xx_z[1][zidx]-sigma_boot_sigma_xx_z[1][zidx])/p_xx_z[1][zidx],
                                            y2=(boot_sigma_xx_z[1][zidx]+sigma_boot_sigma_xx_z[1][zidx])/p_xx_z[1][zidx],
                                            color=colors[1],alpha=0.5)

                ax[i,j].plot(k_x,boot_sigma_xx_z[0][zidx]/p_xx_z[0][zidx],color=colors[0],ls='dotted')
                ax[i,j].plot(k_x,sigma_xx_z[0][zidx]/p_xx_z[0][zidx],color=colors[0],ls='-.')
                ax[i,j].fill_between(x=k_x, y1=(boot_sigma_xx_z[0][zidx]-sigma_boot_sigma_xx_z[0][zidx])/p_xx_z[0][zidx],
                                            y2=(boot_sigma_xx_z[0][zidx]+sigma_boot_sigma_xx_z[0][zidx])/p_xx_z[0][zidx],
                                            color=colors[0],alpha=0.5)



                ax[i,j].text(0.80, 0.90,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
                                                          , ha='center', va='center',
                                                          transform=ax[i,j].transAxes,fontsize=20)
                ax[i,j].xaxis.set_ticks_position('both')
                ax[i,j].yaxis.set_ticks_position('both')
                ax[i,j].xaxis.set_tick_params(direction='in')#, which='top')
                ax[i,j].yaxis.set_tick_params(direction='in')#, which='top')
                #ax[i,j].set_yticks(ax[i,j].get_yticks()[::1])
                #print(zidx,i,j)
                zidx += 1
        # if (i==2)|(i==1)&(j==1):
        #     pass
        # else:
        #     ax[i,j].set_xticklabels([])
        #     ax[i,j].set_xticks([])
# ax[1,2].set_ylim(np.log10(8e-3),np.log10(4e-1))
# ax[0,0].set_xlim(-2.7,-0.9)
# ax[0,0].set_ylim(-2.8,-0.30)
fig.delaxes(ax[0][2])

# if plot_inis:
#     #labels.append("Wilson+19")
#     custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
# box = ax[-1,1].get_position()
# ax[-1,1].set_position([box.x0, box.y0, box.width, box.height])
# ax[-1,1].legend(custom_lines,labels,fontsize=25,loc='center left', bbox_to_anchor=(1.1, 0.5),frameon=False)
ax[0,0].legend(custom_lines[:2],labels[:2],loc='upper left',ncol=1,frameon=False)
ax[0,1].legend(custom_lines[2:3],labels[2:3],loc='upper left',ncol=1,frameon=False)
ax[1,0].legend(custom_lines[3:5],labels[3:5],loc='upper left',ncol=1,frameon=False)
ax[1,1].legend(custom_lines[5:],labels[5:],loc='upper left',ncol=1,frameon=False)


# for i in ax[0,2].get_xticklabels():
#     i.set_visible(True)
ax[0,2].xaxis.set_tick_params(labelbottom=True)

ax[0,0].set_ylim(0,0.65)

# ax[0,1].set_xticklabels(['-2.5','-2.0','-1.5','-1.0'])

plt.tight_layout()
if inis.save_paper_err_pk:
    plt.savefig(inis.save_paper_err_pk_path)
if show_plot:
    plt.show()
else:
    plt.clf()

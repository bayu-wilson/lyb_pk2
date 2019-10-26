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
          r"$\sigma_{P_{\alpha \beta}}/P_{\alpha \beta}$","realization"]
#r"${\sigma_{\hat P_{\beta \beta}}}}$","data","mocks"]#,r"$P_{\beta \beta}$"]
                #,, r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None),
                #Line2D([0], [0], color=colors[3], lw=9, marker=None),
                Line2D([0], [0], color='k', ls='-'),
                Line2D([0], [0], color='k', ls='--')]


# pkdata = pd.read_csv(inis.save_pk_with_err_path)
# pkdata_mocks = pd.read_csv("../output/pk_errboot_mocks_lyb_nocorr.txt")
# np.savetxt('../output/realizations/mean_pk.txt',np.column_stack(np.reshape(realization_table_pk,(50,21*13))))
# np.savetxt('../output/realizations/var_pk.txt',np.column_stack(np.reshape(realization_table_var,(50,21*13))))

mean_pk_table = np.reshape(np.loadtxt('../output/realizations/mean_pk.txt').T,(50,21,13))
var_pk_table = np.reshape(np.loadtxt('../output/realizations/var_pk.txt').T,(50,21,13))
p_xx_z = np.reshape(np.mean(mean_pk_table,axis=0),(3,7,13)) #paa,ptt,pab for 7 zbins and 13 kbins
sigma_xx_z = np.reshape(np.std(mean_pk_table,axis=0),(3,7,13)) #sigma_aa,sigma_ptt,sigma_pab for 7 zbins and 13 kbins
boot_sigma_xx_z = np.reshape(np.mean(np.sqrt(var_pk_table),axis=0),(3,7,13))
sigma_boot_sigma_xx_z = np.reshape(np.std(np.sqrt(var_pk_table),axis=0),(3,7,13))


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
            k = opt.kbin_centers
            k_x = np.log10(k)

            if (i == 0)&(j==2):
                pass
            else:
                #ax[i,j].plot(k_x,sigma_xx_z[0][zidx]/p_xx_z[0][zidx],color=colors[0])
                #ax[i,j].plot(k_x,sigma_xx_z[0][zidx]/p_xx_z[1][zidx],color=colors[1])
                ax[i,j].plot(k_x,sigma_xx_z[2][zidx]/p_xx_z[2][zidx],color=colors[2])
                ax[i,j].plot(k_x,boot_sigma_xx_z[2][zidx]/p_xx_z[2][zidx],color=colors[2],ls='--')
                ax[i,j].fill_between(x=k_x, y1=(boot_sigma_xx_z[2][zidx]-sigma_boot_sigma_xx_z[2][zidx])/p_xx_z[2][zidx],
                                            y2=(boot_sigma_xx_z[2][zidx]+sigma_boot_sigma_xx_z[2][zidx])/p_xx_z[2][zidx],
                                            color=colors[2],alpha=0.5)





            # (mean_bootstrap_ii - std deviation bootstra_ii) / mean_i and (mean_bootstrap_ii + std deviation bootstra_ii) / mean_i)



                #ax[i,j].plot(k_x,sigma_boot_sigma_xx_z[0][zidx]/boot_sigma_xx_z[2][zidx],color=colors[2],ls='--')
                #ax[i,j].plot(k_x,np.log10(err_pbb),color=colors[3])

                # ### mocks ###
                # ax[i,j].plot(k_x,err_paa_mocks/paa_mocks,color=colors[0],ls='--')
                # ax[i,j].plot(k_x,err_ptt_mocks/ptt_mocks,color=colors[1],ls='--')
                # ax[i,j].plot(k_x,err_pab_mocks/pab_mocks,color=colors[2],ls='--')
                #ax[i,j].plot(k_x,np.log10(err_pbb_mocks),color=colors[3],ls='--')
                #ax[i,j].set_yscale('log')
                # ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.')
                # ax[i,j].errorbar(k_x*0.99,np.log10(k*ptt/np.pi),yerr=log_err_ptt,color=colors[1], fmt='.')
                # ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi),yerr=log_err_pab,color=colors[2], fmt='.')
                #ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi))
                #ax[i,j].errorbar(k_x,k*pbb/np.pi,yerr=err_pbb*k/np.pi,color=colors[3], fmt='.')
                ax[i,j].text(0.80, 0.90,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
                                                          , ha='center', va='center',
                                                          transform=ax[i,j].transAxes,fontsize=20)
                ax[i,j].xaxis.set_ticks_position('both')
                ax[i,j].yaxis.set_ticks_position('both')
                ax[i,j].xaxis.set_tick_params(direction='in')#, which='top')
                ax[i,j].yaxis.set_tick_params(direction='in')#, which='top')

                zidx += 1

fig.delaxes(ax[0][2])

# if plot_inis:
#     #labels.append("Wilson+19")
#     custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
# box = ax[-1,1].get_position()
# ax[-1,1].set_position([box.x0, box.y0, box.width, box.height])
# ax[-1,1].legend(custom_lines,labels,fontsize=25,loc='center left', bbox_to_anchor=(1.1, 0.5),frameon=False)
ax[0,0].legend(custom_lines[:2],labels[:2],loc='upper left',ncol=1,frameon=False)
ax[0,1].legend(custom_lines[2:3],labels[2:3],loc='upper left',ncol=1,frameon=False)
ax[1,0].legend(custom_lines[3:4],labels[3:4],loc='upper left',ncol=1,frameon=False)


# for i in ax[0,2].get_xticklabels():
#     i.set_visible(True)
ax[0,2].xaxis.set_tick_params(labelbottom=True)

ax[0,0].set_ylim(0,0.65)

# ax[0,1].set_xticklabels(['-2.5','-2.0','-1.5','-1.0'])

plt.tight_layout()
if inis.save_paper_realization_err_pk:
    plt.savefig(inis.save_paper_realization_err_pk_path)
if show_plot:
    plt.show()
else:
    plt.clf()

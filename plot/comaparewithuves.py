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
labels = [r"$P_{\alpha \alpha}$",r"$P_{TT}$",r"$P_{\alpha \beta}$"]#,r"$P_{\beta \beta}$"]
                #,, r"$P_{\alpha \beta}$",r"$P_{\beta \beta}$"]
custom_lines = [Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[2], lw=9, marker=None)]
                # Line2D([0], [0], color=colors[3], lw=9, marker=None)]


pkdata = pd.read_csv("../output/pk_errboot_obs_lyb_wNR_DLATrue_metalFalse_res0_nb100.csv") #../output/dla_comparison/pk_errboot_obs_corrNR_keepingDLAs.txt")#inis.save_pk_with_err_path)
pkdata2 = pd.read_csv("../output_uves/pk_errboot_obs_lyb_nocorr_DLAFalse_metalFalse_SNR20_nb1000.csv")   #"../output/dla_comparison/pk_errboot_obs_corrNR_removingDLAs.txt")#inis.save_pk_with_err_path)


# plt.style.use('classic')
fig,ax = plt.subplots(2,3,sharex=True, sharey=True,gridspec_kw={'hspace': 0,'wspace': 0}) #ROW COLUMN
# ax2 = ax.twinx()

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# fig.text(0.5, 0.04, r"$log_{10}(k/[km^{-1}s])$", ha='center', va='center',fontsize=20)
fig.text(0.02, 0.5, r"$log_{10}(k P(k) / \pi)$", ha='center', va='center', rotation='vertical',fontsize=25)
plt.grid(False)
# plt.xlabel(r"log$_{10}(k/[km^{-1}s])$")
# plt.ylabel(r"$log_{10}(k P(k) / \pi)$")

# fig.delaxes(ax[-1][2])
fig.set_size_inches(15,9)
ax[-1,1].set_xlabel(r"$log_{10}(k/[km^{-1}s])$",fontsize=25)
#[i.set_xlabel(r"log$_{10}(k/[km^{-1}s])$") for i in ax[-1]] #xlabels on bottom row
# ax[1,1].set_ylabel(r"$k P(k) / \pi$")

zidx = 0
for i in range(2):
    for j in range(3):
        #plot Walther data
        if opt.zbin_centers[zidx] <= 3.4:
            walther = np.loadtxt("waltherall.dat")

            def find_nearest(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return array[idx]
            
            znearest = find_nearest(walther[:,1], opt.zbin_centers[zidx])
            print("znearest = ", znearest, np.log10(walther[:,2]),np.log10(walther[:,3]))
            zwmask = (znearest == walther[:,1])
            ax[i,j].errorbar(np.log10(walther[zwmask,2]),np.log10(walther[zwmask,3]),yerr=np.log10((walther[zwmask,4]+walther[zwmask,3])/walther[zwmask,3]),color='orange', fmt='.-',alpha = 0.3)
            
        if zidx<7:
            pkmask = (pkdata.z == opt.zbin_centers[zidx])
            t = pkdata[pkmask]
            t2 = pkdata2[pkmask]
            k,paa,ptt,pab,pbb = t.k,t.paa,t.ptt,t.pab,t.pbb
            k2,paa2,ptt2,pab2,pbb2 = t2.k,t2.paa,t2.ptt,t2.pab,t2.pbb

            err_paa,err_ptt,err_pab,err_pbb= t.err_paa, t.err_ptt, t.err_pab, t.err_pbb
            err_paa2,err_ptt2,err_pab2,err_pbb2= t2.err_paa, t2.err_ptt, t2.err_pab, t2.err_pbb

            ### INCLUDING RESOLUTION UNCERTAINTY
            # sigma_res_aa = opt.find_res_uncertainty(k.values,t.z.values,paa.values)
            # err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
            # sigma_res_tt = opt.find_res_uncertainty(k.values,t.z.values,ptt.values)
            # err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
            # sigma_res_ab = opt.find_res_uncertainty(k.values,t.z.values,pab.values)
            # err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
            # sigma_res_bb = opt.find_res_uncertainty(k.values,t.z.values,pbb.values)
            # err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)


            N_aa, N_tt, N_ab= t.N_aa, t.N_tt, t.N_ab
            N_aa2, N_tt2, N_ab2= t2.N_aa, t2.N_tt, t2.N_ab

            log_err_paa = 0.434*err_paa/paa
            log_err_ptt = 0.434*err_ptt/ptt
            log_err_pab = 0.434*err_pab/pab

            log_err_paa2 = 0.434*err_paa2/paa2
            log_err_ptt2 = 0.434*err_ptt2/ptt2
            log_err_pab2 = 0.434*err_pab2/pab2
            k_x = np.log10(k)

            if (i == 0)&(j==2):
                pass
            else:
                #ax[i,j].set_yscale('log')
                print("paa ratio = ", paa2/paa)
                ax[i,j].errorbar(k_x,np.log10(k*paa/np.pi),yerr=log_err_paa,color=colors[0], fmt='.',alpha = 0.5)
                ax[i,j].errorbar(k_x,np.log10(k*ptt/np.pi),yerr=log_err_ptt,color=colors[1], fmt='.',alpha = 0.5)
                ax[i,j].errorbar(k_x,np.log10(k*pab/np.pi),yerr=log_err_pab,color=colors[2], fmt='.',alpha = 0.5)

                ax[i,j].errorbar(k_x,np.log10(k*paa2/np.pi),yerr=log_err_paa2,color="gray", fmt='-',alpha = 0.5)
                ax[i,j].errorbar(k_x,np.log10(k*ptt2/np.pi),yerr=log_err_ptt2,color="gray", fmt='--',alpha = 0.5)
                ax[i,j].errorbar(k_x,np.log10(k*pab2/np.pi),yerr=log_err_pab2,color="gray", fmt='.-',alpha = 0.5)
                #ax[i,j].errorbar(k_x,k*pbb/np.pi,yerr=err_pbb*k/np.pi,color=colors[3], fmt='.')
                ax[i,j].text(0.50, 0.10,"z={0}".format(opt.zbin_centers[zidx]) #0.80, 0.95
                                                          , ha='center', va='center',
                                                          transform=ax[i,j].transAxes,fontsize=20)
                #ax[i,j].xaxis.set_ticks_position('both')
                #ax[i,j].yaxis.set_ticks_position('both')
                if j!=0:
                    ax[i,j].xaxis.set_tick_params(length=0)#, which='top')
                    ax[i,j].yaxis.set_tick_params(length=0)#, which='top')
                else:
                    ax[i,j].xaxis.set_tick_params(direction='in')#, which='top')
                    ax[i,j].yaxis.set_tick_params(direction='in')#, which='top')

                ax2 = ax[i,j].twinx()
                #ax2.set_yscale('log')
                #ax2.plot(k_x,(N_aa2-N_aa)/N_aa2,color=colors[0])
                #ax2.plot(k_x,(N_tt2-N_tt)/N_tt2,color=colors[1])
                #ax2.plot(k_x,(N_ab2-N_ab)/N_ab2,color=colors[2])

                #ax2.plot(k_x,N_aa2,color="gray")
                #ax2.plot(k_x,N_tt2,color="gray")
                #ax2.plot(k_x,N_ab2,color="gray")
                #ax2.set_ylim(1e1,2200)
                ax2.set_ylim(-.02,1)


                #ax[i,j].set_yticks(ax[i,j].get_yticks()[::1])
                #print(zidx,i,j)
                zidx += 1
        # if (i==2)|(i==1)&(j==1):
        #     pass
        # else:
        #     ax[i,j].set_xticklabels([])
        #     ax[i,j].set_xticks([])
# ax[1,2].set_ylim(np.log10(8e-3),np.log10(4e-1))
ax[0,0].set_xlim(-2.7,-0.9)
ax[0,0].set_ylim(-2.8,-0.30)
ax[0,1].text(1.2, 0.5,r"frac. pixels lost", ha='center', va='center', transform=ax[0,1].transAxes,fontsize=20,rotation=-90)

fig.delaxes(ax[0][2])

# if plot_inis:
#     #labels.append("Wilson+19")
#     custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='.'))
# box = ax[-1,1].get_position()
# ax[-1,1].set_position([box.x0, box.y0, box.width, box.height])
# ax[-1,1].legend(custom_lines,labels,fontsize=25,loc='center left', bbox_to_anchor=(1.1, 0.5),frameon=False)
ax[1,0].legend(custom_lines[:2],labels[:2],loc='lower left',ncol=1,frameon=False)
ax[1,1].legend(custom_lines[2:4],labels[2:4],loc='lower left',ncol=1,frameon=False)


# for i in ax[0,2].get_xticklabels():
#     i.set_visible(True)
ax[0,2].xaxis.set_tick_params(labelbottom=True)

# ax[0,1].set_xticklabels(['-2.5','-2.0','-1.5','-1.0'])

plt.tight_layout()

figname = "figures/uvesXShooterComp2.pdf"
print("printing to ", figname)
plt.savefig(figname)

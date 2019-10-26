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

pkdata = pd.read_csv(inis.save_pk_path)
if inis.mock_or_obs == "obs":
    # mf = np.loadtxt(inis.save_boot_mf_path)
    # mf_bootstrap = np.reshape(y,(12,opt.zbinlen,inis.M))
    # err_mfa = np.cov(mf_bootstrap[0]).diagonal()**0.5
    # err_mft = np.cov(mf_bootstrap[0]).diagonal()**0.5
    # err_mfb = np.cov(mf_bootstrap[0]).diagonal()**0.5
    # if inis.save_mf_with_err:
    #     mf_everything = np.column_stack((pkdata.z,
    #                     pkdata.Paa, np.concatenate(err_paa),
    #                     pkdata.Ptot, np.concatenate(err_ptt),
    #                     pkdata.Pab,np.concatenate(err_pab)))
    #
    #     np.savetxt(inis.save_pk_with_err_path,pk_everything)


    p = np.loadtxt(inis.save_boot_pk_path)
    pk_bootstrap = np.reshape(p,(8,opt.kbinlen*opt.zbinlen,inis.M))
    err_paa = np.reshape(np.cov(pk_bootstrap[1]).diagonal()**0.5,
             (opt.zbinlen,opt.kbinlen))
    err_ptt = np.reshape(np.cov(pk_bootstrap[2]).diagonal()**0.5,
             (opt.zbinlen,opt.kbinlen))
    err_pab = np.reshape(np.cov(pk_bootstrap[3]).diagonal()**0.5,
             (opt.zbinlen,opt.kbinlen))
    # if inis.save_pk_with_err:
    #     # pk_everything = np.column_stack((pkdata.k,pkdata.z,
    #     #                 pkdata.Paa, pkdata.Ptot, pkdata.Pab,
    #     #                 np.concatenate(err_paa),
    #     #                 np.concatenate(err_ptt),
    #     #                 np.concatenate(err_pab)))
    #     pk_everything = np.column_stack((pkdata.k,pkdata.z,
    #                     pkdata.Paa, np.concatenate(err_paa),
    #                     pkdata.Ptot, np.concatenate(err_ptt),
    #                     pkdata.Pab,np.concatenate(err_pab)))
    #
    #     np.savetxt(inis.save_pk_with_err_path,pk_everything)


fig,ax = plt.subplots(opt.zbinlen)
fig.set_size_inches(12,11*opt.zbinlen/1.5)
ax[-1].set_xlabel(r"log$_{10}(k/[km^{-1}s])$")

for zidx in range(opt.zbinlen):#range(4,7):#(4,7):#opt.zbinlen):
    ax[zidx].set_ylabel(r"$k P(k) / \pi$")
    ax[zidx].set_yscale('log')
    labels = ["Paa","Ptt","Pab"]
    #[r"$P_{\alpha \alpha}$",r"$P_{TT}$",r"$P_{T\alpha}$"]
    custom_lines = [Line2D([0], [0], color=colors[0], lw=3, marker=None),
                    Line2D([0], [0], color=colors[1], lw=3, marker=None),
                    Line2D([0], [0], color=colors[2], lw=3, marker=None)]
    if plot_standard:
        if inis.mock_or_obs == "mocks":
            pkmask = (pkdata.z == opt.zbin_centers[zidx])
            t = pkdata[pkmask]
            k,paa,ptot,pat = t.k,t.Paa,t.Ptot,t.Pab
            k_x = np.log10(k)
            ax[zidx].plot(k_x,k*paa/np.pi,color=colors[0])
            ax[zidx].plot(k_x,k*ptot/np.pi,color=colors[1])

            pos_mask = pat>0
            ax[zidx].scatter(k_x[pos_mask],k[pos_mask]*pat[pos_mask]/np.pi,color=colors[2])
            ax[zidx].scatter(k_x[~pos_mask],-1*k[~pos_mask]*pat[~pos_mask]/np.pi,color=colors[2],marker='^')

            # mask= y>0
            # plt.scatter(x[mask],y[mask])
            # plt.scatter(x[~mask],-1*y[~mask],marker='^')
            labels.append(inis.tag)
            custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='-'))
        if inis.mock_or_obs == "obs":
            pkmask = (pkdata.z == opt.zbin_centers[zidx])
            t = pkdata[pkmask]
            k,paa,ptot,pat = t.k,t.Paa,t.Ptot,t.Pab
            k_x = np.log10(k)
            ax[zidx].errorbar(k_x,k*paa/np.pi,yerr=err_paa[zidx]*k/np.pi,color=colors[0], fmt='o')
            #ax[zidx].errorbar(k_x,k*ptot/np.pi,yerr=err_ptt[zidx]*k/np.pi,color=colors[1], fmt='o')
            ax[zidx].errorbar(k_x,k*pat/np.pi,yerr=err_pab[zidx]*k/np.pi,color=colors[2], fmt='o')
            #labels.append(inis.tag)
            labels.append("Wilson+19")
            custom_lines.append(Line2D([0], [0], color='k',lw=0,marker='o'))
    if plot_bb:
        pkmask = (pkdata.z == opt.zbin_centers[zidx])
        t = pkdata[pkmask]
        k,pbb = t.k,t.Pbb
        k_x = np.log10(k)
        ax[zidx].plot(k_x,k*pbb/np.pi,color=colors[4])
        labels.append(inis.tag)
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='-'))


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
        ax[zidx].plot(k_sim_x,k_sim*t.pta[k_mask]/np.pi,marker='None', color=colors[2],ls='--')
        labels.append("Sherwood Sims (3/19)")
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='dashed'))
    if (plot_sims_symp) & (zidx==4):
        T0,gamma,pi = 7197.95, 0.988483, np.pi
        path = "../data/sims/symp_model/pkT_LCDM_COLDg1.0_f1.0_z3.800.txt"
        pk_cold_g10_z38 = np.loadtxt(path,skiprows=4)
        ## defines the window function
        def wink(k,T1,g,line='lya'):
            if (line=='lya'):
                dT = T1*3**(g-1.)
            elif (line=='lyb'):
                dT = T1*12**(g-1.)
            sigma = 14.0 * np.sqrt(dT/1e4) # in km/s
            print(np.sqrt(T0*T0 + dT*dT), np.sqrt(T0*T0 + T1*T1))
            return np.exp(-k*k*sigma*sigma)
        kp = pk_cold_g10_z38[:,0]
        ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1],color='blue',linestyle='-.',alpha=1)
        #ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1]*wink(kp,5000,1.0,'lya'),color='blue',linestyle=':',alpha=1)
        ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1]*wink(kp,5000,1.6,'lya'),color='blue',linestyle='-',alpha=1)
        ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1]*wink(kp,10000,1.6,'lya'),color='blue',linestyle='--',alpha=1)

        ## T-alpha
        ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3],color='green',linestyle='-.',alpha=1)
        #ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3]*wink(kp,5000,1.0,'lyb'),color='green',linestyle=':',alpha=1)
        ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3]*wink(kp,5000,1.6,'lyb'),color='green',linestyle='-',alpha=1)
        ax[zidx].plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3]*wink(kp,10000,1.6,'lyb'),color='green',linestyle='--',alpha=1)
        ax[zidx].set_xlim(-2.5,-1.0)
        ax[zidx].set_ylim(1e-2,3e-1)
        ax[zidx].set_yticks([])
        ax[zidx].set_xticks([])
        ax[zidx].set_ylabel('')

        labels.append(r"T=7,000 K")
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='-.'))
        #labels.append(r"T=9,000 K")
        #custom_lines.append(Line2D([0], [0], color='k', lw=1,ls=':'))
        labels.append(r"T=9,000 K")
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='-'))
        labels.append(r"T=12,000 K")
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='--'))


    if (plot_sims_thermal) & (zidx>3):
        new_idx = zidx - 4
        ref_path_hot = glob.glob("../data/sims/thermal_models/pkT_LCDM_HOT_f1.0_z*.txt")
        ref_path_cold = glob.glob("../data/sims/thermal_models/pkT_LCDM_COLD_f1.0_z*.txt")
        ref_path_neutral = glob.glob("../data/sims/thermal_models/pkT_LCDM_f1.0_z*.txt")

        ref_path_hot.sort()
        ref_path_cold.sort()
        ref_path_neutral.sort()

        t_hot=pd.read_csv(ref_path_hot[new_idx],skiprows=4, delim_whitespace=True,
                       names=["k","paa","pbb","pab","pzz","ptt","pta","pza","pzb","pzpb"])
        t_cold=pd.read_csv(ref_path_cold[new_idx],skiprows=4, delim_whitespace=True,
                       names=["k","paa","pbb","pab","pzz","ptt","pta","pza","pzb","pzpb"])
        t_neu=pd.read_csv(ref_path_neutral[new_idx],skiprows=4, delim_whitespace=True,
                       names=["k","paa","pbb","pab","pzz","ptt","pta","pza","pzb","pzpb"])


        k_sim = t_neu.k#,t.pk
        k_mask = (k_sim<=opt.kmax)#&(k_sim>=opt.kmin)
        k_sim = k_sim[k_mask]
        k_sim_x = np.log10(k_sim)

        ax[zidx].plot(k_sim_x,k_sim*t_hot.paa[k_mask]/np.pi,marker='None', color=colors[0],ls='--')
        ax[zidx].plot(k_sim_x,k_sim*(t_hot.ptt)[k_mask]/np.pi,marker='None', color=colors[1],ls='--')
        ax[zidx].plot(k_sim_x,k_sim*t_hot.pta[k_mask]/np.pi,marker='None', color=colors[2],ls='--')

        ax[zidx].plot(k_sim_x,k_sim*t_cold.paa[k_mask]/np.pi,marker='None', color=colors[0],ls='--')
        ax[zidx].plot(k_sim_x,k_sim*(t_cold.ptt)[k_mask]/np.pi,marker='None', color=colors[1],ls='--')
        ax[zidx].plot(k_sim_x,k_sim*t_cold.pta[k_mask]/np.pi,marker='None', color=colors[2],ls='--')

        ax[zidx].plot(k_sim_x,k_sim*t_neu.paa[k_mask]/np.pi,marker='None', color=colors[0],ls='--')
        ax[zidx].plot(k_sim_x,k_sim*(t_neu.ptt)[k_mask]/np.pi,marker='None', color=colors[1],ls='--')
        ax[zidx].plot(k_sim_x,k_sim*t_neu.pta[k_mask]/np.pi,marker='None', color=colors[2],ls='--')



        labels.append("Simulations")
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='dashed'))

    if (plot_sims_feb) & (zidx<4):
        ref_path = glob.glob("../data/sims/pk_LCDM_REF_Feb5/pkb_LCDM_f1.0_z*.txt")
        ref_path.sort()
        t=pd.read_csv(ref_path[zidx],skiprows=4, delim_whitespace=True,
                        names=["k","paa",'pbb','pab','qab'])
        k_sim = t.k#,t.pk
        k_mask = (k_sim<=opt.kmax)#&(k_sim>=opt.kmin)
        k_sim = k_sim[k_mask]
        k_sim_x = np.log10(k_sim)

        ax[zidx].plot(k_sim_x,k_sim*t.paa[k_mask]/np.pi,marker='None', color=colors[0],ls='dotted')
        #ax.plot(k_sim_x,k_sim*t.pbb[k_mask]/np.pi,marker='None', color=colors[1],ls='dotted')
        #ax.plot(k_sim_x,k_sim*t.pab[k_mask]/np.pi,marker='None', color=colors[2],ls='dotted')
        #ax.plot(k_sim,k_sim*t.qab[k_mask]/np.pi,marker='None', color=colors[3],ls='dotted')
        labels.append("Sherwood Sims (2/19)")
        custom_lines.append(Line2D([0], [0], color='k', lw=1,ls='dotted'))




    ax[zidx].text(0.80, 0.95,"z={0}".format(opt.zbin_centers[zidx])
                                              , ha='center', va='center', transform=ax[zidx].transAxes)
    box = ax[zidx].get_position()
    ax[zidx].set_position([box.x0, box.y0, box.width * 0.75, box.height])
    ax[zidx].legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))

if inis.save_pk_fig:
    plt.savefig(inis.save_pk_fig_path)
if show_plot:
    plt.show()
plt.clf()
# print(pkdata[pkmask])

# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D
# import matplotlib
# import glob
# plot_sim_feb = False
# plot_sim_march = True
# plot_dataset = True
# plot_boss = False
# save_fig=False
#
#
# # plot_vid_obs = False
# # plot_vid_mocks = False
# # plot_BOSS = False
#
#
# colors = ['red', 'cyan','purple','orange','gold','green','indigo','black','gray']
# marker = ['s','D','^','d','*']
# #user_input_fnames = ['lyb_nocorr_n5000 [Bayu]']
# #len_files = len(user_input_fnames)
#
# SMALL_SIZE = 12
# MEDIUM_SIZE = 15
# BIGGER_SIZE = 20
#
# plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)
#
# matplotlib.rcParams['xtick.major.size'] = 10
# matplotlib.rcParams['xtick.major.width'] = 1.5
# matplotlib.rcParams['xtick.minor.size'] = 5
# matplotlib.rcParams['xtick.minor.width'] = 1.5
# matplotlib.rcParams['ytick.major.size'] = 10
# matplotlib.rcParams['ytick.major.width'] = 1.5
# matplotlib.rcParams['errorbar.capsize'] = 4
#
#
#
# for i in range(4,7):#(4,7):#opt.zbinlen):
#     fig,ax = plt.subplots(1)
#     fig.set_size_inches(12,8)
#     ax.set_ylabel(r"$k P(k) / \pi$")
#     ax.set_xlabel(r"log$_{10}(k/[km^{-1}s])$")
#     #
#     ax.set_yscale('log')
#
#     if plot_dataset:
#         power_table = pd.DataFrame(test_mat[i].T,columns=['k','Paa','Ptot','Pab','Qab','Pbb','num'])
#         k = power_table.k
#         k_x = np.log10(k)
#         paa = power_table.Paa
#         pat = power_table.Pab
#         qat = power_table.Qab
#         ptt = power_table.Ptot
#         pbb = power_table.Pbb
#         ax.plot(k_x,k*paa/np.pi,color=colors[0])
#         ax.plot(k_x,k*ptt/np.pi,color=colors[1])
#         ax.plot(k_x,k*pat/np.pi,color=colors[2])
#         #ax.plot(k,k*pbb/np.pi,color=colors[2])
#         #ax.plot(k,k*qab/np.pi,color=colors[3])
#
#
# #     if plot_BOSS:  #USING BOSS MF
# #         t = pd.read_csv("ptotal_BOSS_lybforest_z{}.csv".format(opt.zbin_centers[i]))
# #         k_BOSS = t.k.values
# #         paa_BOSS = t.Paa.values
# #         pbb_BOSS = t.Pbb.values
# #         ax.plot(k_BOSS,k_BOSS*t.Paa.values/np.pi, color=colors[0],ls='dashdot')
# #         ax.plot(k_BOSS,k_BOSS*t.Pab.values/np.pi, color=colors[1],ls='dashdot')
# #         ax.plot(k_BOSS,k_BOSS*t.Pbb.values/np.pi, color=colors[2],ls='dashdot')
# #         ax.plot(k_BOSS,k_BOSS*t.Ptot.values/np.pi, color=colors[4],ls='dashdot')
#
#     if plot_sim_feb:
#         ref_path = glob.glob("../data/pk_LCDM_REF_Feb5/pkb_LCDM_f1.0_z*.txt")
#         ref_path.sort()
#         t=pd.read_csv(ref_path[i],skiprows=4, delim_whitespace=True,
#                         names=["k","paa",'pbb','pab','qab'])
#         k_sim = t.k#,t.pk
#         k_mask = (k_sim<=opt.kmax)&(k_sim>=opt.kmin)
#         k_sim = k_sim[k_mask]
#         k_sim_x = np.log10(k_sim)
#
#         ax.plot(k_sim_x,k_sim*t.paa[k_mask]/np.pi,marker='None', color=colors[0],ls='dotted')
#         #ax.plot(k_sim_x,k_sim*t.pbb[k_mask]/np.pi,marker='None', color=colors[1],ls='dotted')
#         ax.plot(k_sim_x,k_sim*t.pab[k_mask]/np.pi,marker='None', color=colors[2],ls='dotted')
#         #ax.plot(k_sim,k_sim*t.qab[k_mask]/np.pi,marker='None', color=colors[3],ls='dotted')
#
#     if (plot_sim_march) & (i>3):
#             new_idx = i - 4
#             ref_path = glob.glob("../data/sims/pkT_LCDM_f1.0_z*.txt")
#             ref_path.sort()
#             t=pd.read_csv(ref_path[new_idx],skiprows=4, delim_whitespace=True,
#                            names=["k","paa","pbb","pab","pzz","ptt","pta","pza","pzb","pzpb"])
#             k_sim = t.k#,t.pk
#             k_mask = (k_sim<=opt.kmax)#&(k_sim>=opt.kmin)
#             k_sim = k_sim[k_mask]
#             k_sim_x = np.log10(k_sim)
#
#             ax.plot(k_sim_x,k_sim*t.paa[k_mask]/np.pi,marker='None', color=colors[0],ls='--')
#             ax.plot(k_sim_x,k_sim*(t.ptt)[k_mask]/np.pi,marker='None', color=colors[1],ls='--')
#             ax.plot(k_sim_x,k_sim*t.pta[k_mask]/np.pi,marker='None', color=colors[2],ls='--')
#
#             #p_test = t.pzpb[k_mask]#-t.pab[k_mask] #t.pzz[k_mask] + t.pbb[k_mask] + t.pzpb[k_mask] +
#             #ax.plot(k_sim,k_sim*p_test/np.pi,marker='None', color=colors[3],ls='--')
#
#
#
#             #ax.plot(k_sim,k_sim*t.pzz[k_mask]/np.pi,marker='None', color=colors[3],ls='--')
#
#             #test = t.pbb + t.pzz #+ t.pzpb
#             #ax.plot(k_sim,k_sim*test[k_mask]/np.pi,marker='None', color=colors[4],ls='--')
#         #ax.plot(k_sim,k_sim*t.qab[k_mask]/np.pi,marker='None', color=colors[3],ls='dotted')
#
#
#     custom_lines = [Line2D([0], [0], color=colors[0], lw=3, marker=None),
#                     Line2D([0], [0], color=colors[1], lw=3, marker=None),
#                     Line2D([0], [0], color=colors[2], lw=3, marker=None),
#                     Line2D([0], [0], color='k', lw=1,ls='-'),
#                     Line2D([0], [0], color='k', lw=1,ls='dashed')]
#     labels = [r"$P_{\alpha \alpha}$",r"$P_{TT}$",r"$P_{T\alpha}$",
#               "log kbin "+tag,"Sherwood Sims (3/14)"]
# #     # Shrink current axis by 20%
#     box = ax.get_position()
# #     # Put a legend to the right of the current axis
#     ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
#     ax.legend(custom_lines,labels,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))
#     ax.text(0.85, 0.95,"z={0}".format(opt.zbin_centers[i]
#                                           ), ha='center', va='center', transform=ax.transAxes)
#
#     print("test_output/pk_logbins_{0}_z{1:.1f}.pdf".format(tag,opt.zbin_centers[i]))
#     if save_fig:
#         plt.savefig("test_output/pk_logbins_{0}_z{1:.1f}.pdf".format(tag,opt.zbin_centers[i]))
#     #print(ptot.values)
#     #np.savetxt()
#     #print("ptotal_lybforest_z{}.csv".format(opt.zbin_centers[i]))
#     #save_this = power_table[['k','Ptot']]
#     #power_table.to_csv("ptotal_BOSS_lybforest_z{}.csv".format(opt.zbin_centers[i]), index=False)
#     #ax.vlines(opt.kbin_edges,0,0.1)

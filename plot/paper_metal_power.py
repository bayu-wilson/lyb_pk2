#!/usr/bin/env python

import pandas as pd
import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
import glob

colors = ["red","green","blue"]
fontsize = 15
fontsize_big = 20

custom_lines = [Line2D([0], [0], color=colors[1], lw=9, marker=None),
                Line2D([0], [0], color=colors[0], lw=9, marker=None),
                Line2D([0], [0], color='k', lw=9, marker=None),
                Line2D([0], [0], color='grey',lw=9,ls='--',marker=None)]
legend_labels = ["BOSS+13","Iršič+17","Avg. Iršič+17","0.0015cos(500k)+0.008"]

names_boss = ["k","paa","pm","err_paa","err2"]
data_boss = glob.glob("../data/BOSS_metals/pk_z*.txt")
data_boss.sort()
data_boss = data_boss[4:-1]
# del data_boss[-1]

names_sdss = ["k","paa","pm","err_paa"]
data_sdss = glob.glob("../data/BOSS_metals/pk-pat_z*.txt")
data_sdss.sort()
# print(data_sdss)
# del data_sdss[-2]
# del data_sdss[1]

data_xshoot = pd.read_csv("../data/obs/pk_xs_final.csv")

data_xshoot_avg = np.loadtxt('../data/obs/pk_xs_avg.txt')

fig,ax = plt.subplots(1)
fig.set_size_inches(9,7)

ax.set_xlabel(r"log$_{10}(k/[km^{-1}s])$",fontsize=fontsize_big)
ax.set_ylabel(r"$log_{10}(k P(k) / \pi)$",fontsize=fontsize_big)

kmin = opt.kmin #np.min(opt.kbin_centers)
kmax = 10**-1.92

for b in data_boss:
    b_t = pd.read_csv(b,delim_whitespace=True,names=names_boss)
    kvals = b_t.k.values
    kmask = (kvals>kmin)&(kvals<kmax)
    kvals = kvals[kmask]
    #s_t = pd.read_csv(s,delim_whitespace=True,names=names_sdss)
    k_x = np.log10(kvals)
    ax.plot(k_x,np.log10(b_t.pm.values[kmask]*kvals), color = colors[1], alpha=0.5,ls='-.')
    #plt.plot(s_t.k,s_t.pm, color = colors[2], alpha=0.5)

for zidx in range(4,7):
    pkmask = (data_xshoot.z == opt.zbin_centers[zidx])
    t = data_xshoot[pkmask]
    k_x = np.log10(t.k)
    ax.plot(k_x,np.log10(t.pm*t.k), color = colors[0])

# ax.set_ylim(0,0.016)

################################################################################
################################################################################
#### Adding standard deviation... probably remove this ###
from scipy import interp
# data_xshoot = pd.read_csv("../data/obs/pk_xs_final.csv")
# np.shape(data_xshoot[data_xshoot.z>=3.4])
# 95/19
# 19*5*7
x = np.reshape(data_xshoot[data_xshoot.z>=3.4].values.T,(7,5,19))
mylist = []
for i in [-1,-2,-3]:
    mylist.append(interp(opt.kbin_centers,x[1][i], x[5][i]))
    #plt.plot(interp(opt.kbin_centers,x[1][i], x[5][i])*opt.kbin_centers)
# np.mean(mylist,axis=0)
data_xshoot_std = np.std(mylist,axis=0)
log_err_xshoot = [np.abs(np.log10((data_xshoot_avg-data_xshoot_std))-np.log10(data_xshoot_std)),
                            np.abs(np.log10((data_xshoot_avg+data_xshoot_std))-np.log10(data_xshoot_avg))]
################################################################################
################################################################################

ax.errorbar(np.log10(opt.kbin_centers),np.log10(data_xshoot_avg*opt.kbin_centers),
            yerr=log_err_xshoot,color='k',lw=5,elinewidth=1,capsize=3)
# ax.plot(np.log10(opt.kbin_centers), np.log10(data_xshoot_avg*opt.kbin_centers),color='k',lw=5)
fit_x = np.linspace(opt.kmin,opt.kmax,1000)
fit_y = np.cos(500*fit_x-np.pi*0)
ax.plot(np.log10(fit_x),np.log10(0.0015*fit_y+0.008),color='grey',lw=3,ls='--')

ax.legend(custom_lines,legend_labels,fontsize=fontsize,loc='lower left',ncol=2,frameon=False)
ax.tick_params(axis='both', which='major', labelsize=fontsize)
# plt.tight_layout()

if inis.save_paper_metal_power:
    plt.savefig(inis.save_paper_metal_power_path)

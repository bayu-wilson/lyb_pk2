#!/usr/bin/env python

import pandas as pd
import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib
from QuasarSpectrum import QuasarSpectrum

##### Control Area ######
# This program is only meant to plot resolution vs wavelength for wR and WR2
plot_mocks = True
show_plot = False
#########################

colors = ['red', 'orange','purple','cyan','gold','green','indigo','black','gray']
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



cat_name = "mocks/XQ-100_catalogue_n5000"
tags = ["lyb_wR_n5000","lyb_wR2_n5000"]#inis.tag
rescale_flux = 1.0 #0.1 #0.1#0.1 #0.1 #0.0 #0.1

QuasarSpectrum.load_cat(cat_name)
nqso = QuasarSpectrum.nqso
# spec_one = QuasarSpectrum.load_qso_data(0,tag=tag,rescale_flux=1.0)
# print("Loading Data")
# qso_arr = []
q = QuasarSpectrum.load_qso_data(0,tag=tags[0],rescale_flux=rescale_flux)
q2 = QuasarSpectrum.load_qso_data(0,tag=tags[1],rescale_flux=rescale_flux)
res =  q.resolution
res2 = q2.resolution
wav = q.wavelength
# print(min(res2))
# print(res)

fig,ax = plt.subplots(1)
fig.set_size_inches(12,6)

ax.plot(wav,res,color=colors[0])
ax.plot(wav,res2,color=colors[1])
# ax.set_ylim(10.5,11)
ax.set_xlabel("Wavelength [Angstroms]")
ax.set_ylabel(r"$\sigma_R$ [km/s]")


custom_lines = [Line2D([0], [0], color=colors[0], lw=3, marker=None),
                Line2D([0], [0], color=colors[1], lw=3, marker=None)]

plt.tight_layout()

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
ax.legend(custom_lines,tags,fontsize=15,loc='center left', bbox_to_anchor=(1, 0.5))



if inis.save_res_fig:
    plt.savefig(inis.save_res_fig_path)
if show_plot:
    plt.show()



plt.clf()

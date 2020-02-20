#!/usr/bin/env python

from QuasarSpectrum import QuasarSpectrum
import inis
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

import matplotlib.pyplot as plt

cat_name = inis.cat_name
QuasarSpectrum.load_cat(cat_name)
nqso = QuasarSpectrum.nqso
qso_arr = []

for i in range(nqso):
    released_path = "../data/obs/XQ-100/released/{0}_uvb-vis.txt".format(QuasarSpectrum.all_names[i])
    qso_table = pd.read_csv(released_path, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
                                    names = ("wavelength", "flux", "ferr", "res","dloglam"))
    qso_arr.append(qso_table)


zbin_centers = np.arange(3.6,4.8,0.2)
zbin_edges = np.arange(3.5,4.9,0.2)
zbinlen = len(zbin_centers)
wavelength_arr = 10**np.arange(np.log10(978.02938096),np.log10(1600),3e-5)
quasar_sums = np.zeros((zbinlen,len(wavelength_arr)))
quasar_nums = np.zeros(zbinlen)

redside_avg = np.zeros(zbinlen)

for qidx in range(nqso):
    for zidx in range(zbinlen):
        zqso = QuasarSpectrum.all_redshifts[qidx]
        if (zqso>zbin_edges[zidx])&(zqso<zbin_edges[zidx+1]):
            if len(qso_arr[qidx].wavelength)>2982:
                flux = qso_arr[qidx].flux[:-1]
                wave = qso_arr[qidx].wavelength[:-1]/(1+zqso)
                quasar_sums[zidx]+=griddata(points=wave,values=flux,xi=wavelength_arr,method='linear')
                quasar_nums[zidx]+=1
            else: #sometimes the length is 2983??
                flux = qso_arr[qidx].flux
                wave = qso_arr[qidx].wavelength/(1+zqso)
                quasar_sums[zidx]+=griddata(points=wave,values=flux,xi=wavelength_arr,method='linear')
                quasar_nums[zidx]+=1

new_arr = (quasar_sums.T/quasar_nums).T*1e15
mask_wave = (wavelength_arr>1425)&(wavelength_arr<1450)
for zidx in range(zbinlen):
    redside_avg[zidx] = np.mean(new_arr[zidx][mask_wave])

fig,ax = plt.subplots(1)
fig.set_size_inches(10,6)
colors = ['red','orange','green','blue','indigo','purple']
labels = ["{0:.1f}<z<{1:.1f}".format(zbin_edges[i],zbin_edges[i+1]) for i in range(zbinlen)]
for zidx in range(zbinlen-1):
    ax.plot(wavelength_arr,new_arr[zidx],lw=0.75, color=colors[zidx],label=labels[zidx])
ax.yaxis.set_ticks_position('both')
ax.yaxis.set_tick_params(direction='in')
# ax.text(0.82, 0.80,"z={0:.1f}, N={1:d}".format(zbin_centers[zidx],int(quasar_nums[zidx])),ha='center', va='center',
#                                         transform=ax.transAxes,fontsize=18)
# for zidx in range(zbinlen-1):
#     ax[zidx].plot(wavelength_arr,new_arr[zidx],lw=0.75, color='k')
#     ax[zidx].yaxis.set_ticks_position('both')
#     ax[zidx].yaxis.set_tick_params(direction='in')
#     ax[zidx].text(0.82, 0.80,"z={0:.1f}, N={1:d}".format(zbin_centers[zidx],int(quasar_nums[zidx])),ha='center', va='center',
#                                             transform=ax[zidx].transAxes,fontsize=18)
ax.set_ylim(0,4.4)
ax.xaxis.set_ticks_position('bottom')
ax.xaxis.set_tick_params(direction='in')#, which='top')
ax.set_xlabel(r"Wavelength [$\AA$]",fontsize=18)
fig.text(0.07, 0.5, r"Flux [$10^{-16}$ ergs/s/cm${^2}$/$\AA$]", ha='center', va='center', rotation='vertical',fontsize=18)
ax.legend()

if inis.save_composite_spectra:
    plt.savefig(inis.save_composite_spectra_path)
else:
    plt.clf()
if inis.save_redside_avg:
    np.savetxt(inis.save_redside_avg_path,redside_avg)




        # q = QuasarSpectrum.load_qso_data(i,tag=tag,rescale_flux=rescale_flux)
        # if ("noB" in tag)&(inis.add_beta): # adding Lyb
        #     q.get_new_forest(rescale_flux=rescale_flux,wrange =
        #         (opt.lyb_min,opt.lyb_max,opt.lyb_rest, opt.xs_beta))
        # if inis.cat_name.startswith('mocks')&(inis.add_ovi): # adding OVI
        #     q.get_new_forest(rescale_flux=rescale_flux,wrange =
        #         (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d1, opt.xs_ovi))
        #     q.get_new_forest(rescale_flux=rescale_flux,wrange =
        #         (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d2, opt.xs_ovi))
        # if inis.cat_name.startswith('mocks')&(inis.add_sithree): #adding SiIII
        #     q.get_new_forest(rescale_flux=rescale_flux,wrange =
        #         (opt.sithree_min,opt.sithree_max,opt.sithree_rest_d1, opt.xs_sithree))
        #     # q.get_new_forest(rescale_flux=rescale_flux,wrange =
        #     #     (opt.sithree_min,opt.sithree_max,opt.sithree_rest_d2, opt.xs_sithree))

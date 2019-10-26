#!/usr/bin/env python

import matplotlib.pyplot as plt
from QuasarSpectrum import QuasarSpectrum
import inis
import pandas as pd
import glob
import os
import options as opt
import numpy as np


cat_name = "mocks/XQ-100_catalogue_n100.mock" #inis.cat_name
tag = "lyb_nocorr" #inis.tag
rescale_flux = 1.0
nqso = 100

QuasarSpectrum.load_cat(cat_name)

# for zidx in range(opt.zbinlen):
zidx = 2
import sys
for qidx in range(nqso):
    q = QuasarSpectrum.load_qso_data(qidx,tag=tag,rescale_flux=rescale_flux)
    name = q.name
    zpix_a = q.get_zpix(opt.lya_rest)
    zmask_a = q.get_zmask(forest=(opt.lya_min,opt.lya_max,opt.lya_rest),
                                    zpix=zpix_a,zidx=zidx,zedges=opt.zbin_edges,name=name)
    zpix_tot = q.get_zpix(opt.lyb_rest)
    zmask_tot = q.get_zmask(forest=(opt.lyb_min,opt.lyb_max,opt.lyb_rest),
                                    zpix=zpix_tot,zidx=zidx+1,zedges=opt.zbin_edges,name=name)
    nchunks_a = opt.how_many_chunks(zmask_a)
    nchunks_tot = opt.how_many_chunks(zmask_tot)
    fig,ax = plt.subplots(1)
    fig.set_size_inches(11,8.5)
    #plt.figure(figsize=(11,8.5))
    ax.set_title("QUASAR: {0}".format(qidx))
    ax.plot(q.wavelength,q.flux)
    ax.plot(q.wavelength[zmask_a],q.flux[zmask_a], color = 'red', label = 'lya forest (z=3.8)')
    ax.plot(q.wavelength[zmask_tot],q.flux[zmask_tot], color = 'green', label = 'lyb forest (z=4.0)')
    ax.axvspan(opt.lya_rest*(1+opt.zbin_edges[2]),opt.lya_rest*(1+opt.zbin_edges[3]),color='red',alpha = 0.5)
    ax.axvspan(opt.lyb_rest*(1+opt.zbin_edges[3]),opt.lyb_rest*(1+opt.zbin_edges[4]),color='green',alpha = 0.5)
    #plt.axvline(opt.lya_rest*(1+opt.zbin_edges[2]))
    #plt.axvline(opt.lya_rest*(1+opt.zbin_edges[3]))
    ax.set_xlim(4500,6000)
    ax.legend(loc=1)
    # plt.show()
    # sys.exit()
    plt.savefig("sightlines/q{:02}.png".format(qidx))
    plt.clf()

# QuasarSpectrum.load_cat(cat_name)
# catalog_table = pd.read_csv(cat_name,delim_whitespace=True,usecols=(3,6),
#                     names = ("qso_name","redshifts"))
#
# all_paths = glob.glob("../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-*.txt.gz".format(tag))
# all_paths.sort()




# for qidx in range(100):
#     path = all_paths[qidx]
#     qdata = pd.read_csv(path, delim_whitespace=True,compression='gzip',skiprows = 1,usecols=(0,1,2,3,4),
#                             names = ("wav", "flx", "ferr", "res","dloglam"))

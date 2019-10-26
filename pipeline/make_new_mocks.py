#!/usr/bin/env python

from QuasarSpectrum import QuasarSpectrum
import inis
import pandas as pd
import glob
import os
import options as opt
import numpy as np

#########################################################################################################
# When you modify the mocks, you cut out all the pixels between lyman_alpha*(1+z1)
# and lyman_alpha*(1+z2) with z1 and z2  being 3.3 and 3.5 (effectively cutting out redshift bin z=3.4).
#########################################################################################################

cat_name = "../data/mocks/XQ-100_catalogue_n5000.mock" #inis.cat_name
tag = "lyb_nocorr_n5000" #inis.tag
new_tag = "gaussian_n5000"
nqso = 5000
# cat_name = "../data/mocks/XQ-100_catalogue_n100.mock"#inis.cat_name
# tag = "lyb_nocorr" #inis.tag
# new_tag = "gaussian"
# nqso = 100
remove_new_directory = True

# QuasarSpectrum.load_cat(cat_name)
catalog_table = pd.read_csv(cat_name,delim_whitespace=True,usecols=(3,6),
                    names = ("qso_name","redshifts"))

all_paths = glob.glob("../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-*.txt.gz".format(tag))
all_paths.sort()

### MAKING A NEW DIRECTORY FOR NEW MOCKS
if remove_new_directory:
    try:
        os.system("rm -r ../data/mocks/XQ-100_{0}".format(new_tag))
    except:
        pass
os.system("mkdir ../data/mocks/XQ-100_{0}".format(new_tag))
os.system("mkdir ../data/mocks/XQ-100_{0}/lyb".format(new_tag))
os.system("mkdir ../data/mocks/XQ-100_{0}/lyb/ref/".format(new_tag))
os.system("mkdir ../data/mocks/XQ-100_{0}/lyb/ref/spectra".format(new_tag))

if new_tag == "zsubset3.8":
    for qidx in range(100):
        path = all_paths[qidx]
        qdata = pd.read_csv(path, delim_whitespace=True,compression='gzip',skiprows = 1,usecols=(0,1,2,3,4),
                                names = ("wav", "flx", "ferr", "res","dloglam"))
        wave_min = opt.lya_rest*(1+3.7)
        wave_max = opt.lya_rest*(1+3.9)
        mask = (qdata.wav.values>wave_min)&(qdata.wav.values<wave_max)
        #mask = (qdata.wav.values<wave_min)|(qdata.wav.values>wave_max)
        qdata.flx.values[mask] = 1.0#0.0
        new_qdata = np.column_stack((np.zeros(5),qdata.values.T)).T # the zeros are because of the row skip

        four_digit = all_paths[qidx][-28:-24]
        new_path = "../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-{1}_xq100_{0}.txt".format(new_tag,four_digit)
        np.savetxt(new_path,new_qdata)
if "gaussian" in new_tag:
    mu = 0
    sigma=0.1
    np.random.seed(1)
    # np.random.normal(mu, sigma, 100000)
    for qidx in range(nqso):
        path = all_paths[qidx]
        qdata = pd.read_csv(path, delim_whitespace=True,compression='gzip',skiprows = 1,usecols=(0,1,2,3,4),
                                names = ("wav", "flx", "ferr", "res","dloglam"))
        # wave_min = opt.lya_rest*(1+3.7)
        # wave_max = opt.lya_rest*(1+3.9)
        #mask = (qdata.wav.values>wave_min)&(qdata.wav.values<wave_max)
        #mask = (qdata.wav.values<wave_min)|(qdata.wav.values>wave_max)
        #qdata.flx.values[mask] = 1.0#0.0

        qdata.flx = 1-np.random.normal(mu, sigma, len(qdata.flx))
        new_qdata = np.column_stack((np.zeros(5),qdata.values.T)).T # the zeros are because of the row skip

        if "n5000" in new_tag:
            four_digit = all_paths[qidx][-29:-25]
        else:
            four_digit = all_paths[qidx][-28:-24]
        new_path = "../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-{1}_xq100_{0}.txt".format(new_tag,four_digit)
        np.savetxt(new_path,new_qdata)
        os.system("gzip {0}".format(new_path))


# if "n5000" in new_tag:
#     os.system("gzip ../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-0*.txt".format(new_tag,four_digit))
#     os.system("gzip ../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-1*.txt".format(new_tag,four_digit))
#     os.system("gzip ../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-2*.txt".format(new_tag,four_digit))
#     os.system("gzip ../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-3*.txt".format(new_tag,four_digit))
#     os.system("gzip ../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock-4*.txt".format(new_tag,four_digit))
#
# else:
#     os.system("gzip ../data/mocks/XQ-100_{0}/lyb/ref/spectra/mock*.txt".format(new_tag,four_digit))

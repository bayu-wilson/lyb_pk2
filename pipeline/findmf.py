#!/usr/bin/env python

"""
Bayu Wilson, 2018-2019, University of Washington Astronomy

Usage
    ./findmf.py [data origin] [mock correction type] [saving data?]
    ./findmf.py obs uncorr 1
    ./findmf.py obs uncorr 0
    ./findmf.py obs corrNR 1
    ./findmf.py mocks lyb_nocorr 1
    ./findmf.py mocks nocorr 1
    ./findmf.py mocks nocorr_n5000 1
    ./findmf.py mocks lyb_nocorr_n5000 1
    ./findmf.py mocks noB 1
    ./findmf.py mocks noB_onlyA 1

Purpose
   To find mean flux, flux variance, and noise power in each redshift bin.
   Additionally, measures quasar distribution. Note: I only include quasar
   contribution if >100 pixels are involved in analysis.

 NOTES for Matt & Vid:
     1) When saving the tables, my code expects there to exist a directory called,
     `../output/`
     2) My code thinks data is located here: `../data/mocks/*` or `../data/obs/*`
     3) Many of the functions are located in a file called options.py
     4) This code is in Python 3 (not Python 2).
"""

import numpy as np
import pandas as pd
import options as opt
from QuasarSpectrum import QuasarSpectrum

num_imsrquant = 10
#'mf_a', 'nvar_a','dloglambda_a', 'npow_a','mf_tot', 'nvar_tot','dloglambda_tot', 'npow_tot','z','mf_b'
testing_matrix = np.zeros((opt.zbinlen,num_imsrquant))
FLUX_PDF=[[] for tmp in range(opt.zbinlen)]
pdf_bins = np.arange(-0.025,1.05,.05)
# pdf_bins =np.arange(0-0.1/2,1.+0.1/2,.1)

remove_me = np.zeros(opt.zbinlen)#DELETE THIS

for zidx in range(opt.zbinlen):
    #zbin = []

    testing_msrmnt =  np.zeros(num_imsrquant)
    count_testing = np.zeros(num_imsrquant)

    zbin_msrmnt = [[] for idx in range(num_imsrquant)]


    for i in range(nqso):
    #for i in range():
        #qso_arr[i].get_lybforest()
        zpix = qso_arr[i].get_zpix(opt.lya_rest)

        mask = qso_arr[i].get_zmask(forest=(opt.lya_min,opt.lya_max,opt.lya_rest),
                                     zpix=zpix,zidx=zidx,zedges=opt.zbin_edges)

        remove_me[zidx]+=np.sum(mask)
        #print('mock-{0}, {1}:{2},nomask:{3}'.format(i,opt.zbin_centers[zidx],np.sum(mask),len(qso_arr[i].wavelength)))
        #break
        ###
        #FLUX PDF
        FLUX_PDF[zidx].append(np.histogram(qso_arr[i].flux[mask],bins=pdf_bins)[0])#/np.nansum(mask))
        ###

        testing_msrmnt[0] += np.sum(qso_arr[i].flux[mask]) # mf lya
        count_testing[0]+= np.sum(mask) # len lya

        testing_msrmnt[1] += np.sum(qso_arr[i].err_flux[mask]**2) # var lya
        count_testing[1]+= np.sum(mask) # len lya
        testing_msrmnt[2] += np.sum(qso_arr[i].dloglambda[mask]) # dloglam lya
        count_testing[2]+= np.sum(mask) # len lya

        zpix = qso_arr[i].get_zpix(opt.lyb_rest) # Here is where I want to change optical depth of lyb pixels
        mask_b = qso_arr[i].get_zmask(forest=(opt.lyb_min,opt.lyb_max,opt.lyb_rest),
                                     zpix=zpix,zidx=zidx,zedges=opt.zbin_edges)

        #if 'noB' in tag.split('_'):
        #    testing_msrmnt[4] += np.nansum(qso_arr[i].get_lybforest(zidx=zidx,zedges=opt.zbin_edges,
        #                                                            rescale_lyb = rescale_lyb)) #flux[mask])
        #    #get_lybforest(self,zidx,zedges, rescale_lyb = 1/0.1):
        #else:
        testing_msrmnt[4] += np.sum(qso_arr[i].flux[mask_b]) #mf tot
        #print(qso_arr[i].get_lybforest())
        count_testing[4] += np.sum(mask_b) # len tot

        testing_msrmnt[5]+=np.sum(qso_arr[i].err_flux[mask_b]**2) # var tot
        count_testing[5] += np.sum(mask_b) # len tot
        testing_msrmnt[6]+=np.sum(qso_arr[i].dloglambda[mask_b]) # dloglam tot
        count_testing[6] += np.sum(mask_b) # len tot


    testing_matrix[zidx] = testing_msrmnt/count_testing
    #print(count_testing[0],opt.zbin_centers[zidx])
testing_matrix.T[8] = opt.zbin_centers
testing_matrix.T[3] = list(QuasarSpectrum.get_npow(mf=testing_matrix.T[0],
                                        nvar=testing_matrix.T[1],
                                        dloglambda=testing_matrix.T[2]))
testing_matrix.T[7] = list(QuasarSpectrum.get_npow(mf=testing_matrix.T[4],
                                        nvar=testing_matrix.T[5],
                                        dloglambda=testing_matrix.T[6]))
initial_results = pd.DataFrame(testing_matrix)
initial_results.columns = ['mf_a', 'nvar_a','dloglambda_a', 'npow_a',
                           'mf_tot', 'nvar_tot','dloglambda_tot', 'npow_tot',
                            'z','mf_b']


zab_centers = opt.find_za(opt.zbin_centers) #converting lyb zbins to the equivalent, lower, lya zbins
len_zab = len(zab_centers)

#Gives corresponding lya bin for each lyb bin. organized by increasing z.
bin_zab=np.ones(len_zab)*np.nan
for i in range(len_zab):
    for j in range(len_zab):
        if (zab_centers[i]>opt.zbin_edges[j])&(zab_centers[i]<opt.zbin_edges[j+1]):
            bin_zab[i] = (opt.zbin_centers[j])

mf_lyb = np.ones(len_zab)*np.nan #nan until proven otherwise
for i in range(len_zab):
    if bin_zab[i] in initial_results.z.values:
        za_idx = initial_results.z == bin_zab[i]
        ztot_idx = i

        mf_lyb[i] = initial_results.mf_tot[ztot_idx]/initial_results.mf_a[za_idx]
initial_results['mf_b'] = mf_lyb

# # Importing Modules
# import numpy as np
# import pandas as pd
# import options as opt
# import sys
# import glob
# from scipy.interpolate import griddata
#
# # User inputted arguements
# tag = str(sys.argv[1]) #obs,mocks
# corr_tag = str(sys.argv[2]) #lyb_nocorr,wN,wNR,wR,wN_n5000 etc.
# save_data = bool(float(sys.argv[3])) #0,1
# verbose = True
#
# saving_mf_here = "../output/mf_{0}_{1}.csv".format(tag,corr_tag)
# saving_qso_dist_here = "../output/qso_dist_{0}_{1}.csv".format(tag,corr_tag)
#
# if tag=="mocks": #these are to be printed out for UI
#     data_path = "../data/mocks/XQ-100_{0}/*".format(corr_tag)
# if tag=="obs":
#     data_path = "../data/obs/XQ-100/released/*"
#
# print("------------------ Mean Flux ------------------")
# print("Data path:            ", data_path)
# print("mf path:              ", saving_mf_here)
# print("QSO distribution path:", saving_qso_dist_here)
# print("saving data:          ", save_data)
# print("Verbose:              ", verbose)
#
#
# # Loading in related constants
# LYA_REST = opt.lya_rest
# LYA_MIN = opt.lya_min
# LYA_MAX = opt.lya_max
# LYB_REST = opt.lyb_rest
# LYB_MIN = opt.lyb_min
# LYB_MAX = opt.lyb_max
# zbin_centers = opt.zbin_centers
# zbin_edges = opt.zbin_edges
# zbinlen = len(zbin_centers)
#
# # Loading table with quasar names and redshifts
# ### OBS ###
# if tag == "obs":
#     path_to_cat = "../data/obs/XQ-100_catalogue.txt"
#     catalog_table = pd.read_csv(path_to_cat,delim_whitespace=True,usecols=(3,6),
#                                 names = ("qso_name","redshift"))
#     qname_array,z_array = catalog_table['qso_name'].values,catalog_table['redshift'].values
#     nqso = 100 #it will always be 100 quasars fro XQ-100 dataset
#
# ### MOCKS ###
# if tag == "mocks":
#     which_catalog = int(corr_tag.endswith("n5000")) # index for `path_to_cat`
#     path_to_cat = ["../data/mocks/XQ-100_catalogue_n100.mock",
#                    "../data/mocks/XQ-100_catalogue_n5000.mock"][which_catalog]
#     catalog_table = pd.read_csv(path_to_cat,delim_whitespace=True,usecols=(3,6),
#                                 names = ("qso_name","redshift"))
#     qname_array,z_array = catalog_table['qso_name'].values,catalog_table['redshift'].values
#     nqso = int(path_to_cat.split("_n")[1][:-5]) # cutting of _n and .mock gives 100 or 5000
#     print(nqso)
#
# print("Catalog Path:         ", path_to_cat)
# print("\nDoing Calculations ...")
#
# # Grouping flux, flux variance into redshift bins
# # Counting how many quasars contribute >100 pixels into each redshift bin
# bin_aflx = [[] for bin in range(zbinlen)]
# bin_avar = [[] for bin in range(zbinlen)]
# bin_alam = [[] for bin in range(zbinlen)]
# bin_bflx = [[] for bin in range(zbinlen)]
# bin_bvar = [[] for bin in range(zbinlen)]
# bin_blam = [[] for bin in range(zbinlen)]
# bin_aqsonum = np.zeros(7)
# bin_bqsonum = np.zeros(7)
# bin_aqsoname = [[] for bin in range(zbinlen)]
# bin_bqsoname = [[] for bin in range(zbinlen)]
#
# if tag == "mocks":
#     #print(corr_tag)
#     #print(tag)
#     #print("../data/{0}/XQ-100_{1}/*/ref/spectra/*.txt.gz".format(tag,corr_tag))
#     r_filepaths = glob.glob("../data/{0}/XQ-100_{1}/*/ref/spectra/*.txt.gz".format(tag,corr_tag))
#     #print(len(r_filepaths))
#     r_filepaths.sort() #THIS IS VERY IMPORTANT. It must be sorted to match the .mock catalog names and z's.
#
#     print(len(r_filepaths))
#
# for i in range(nqso):
#     if tag == "obs":
#         r_file = "../data/obs/XQ-100/released/{0}_uvb-vis.txt".format(qname_array[i])
#         c_file = "../data/obs/XQ-100/continuum/{0}_cont.txt".format(qname_array[i])
#         r_data = pd.read_csv(r_file, delim_whitespace=True, skiprows = 1, usecols=(0,1,2,4),
#                                         names = ("wav", "flx", "ferr", "dloglam"))
#         c_data = pd.read_csv(c_file, delim_whitespace=True, skiprows = 1,usecols=[0,1], names = ("wav",'flx'))
#         wave = r_data.wav.values
#         flux = r_data.flx.values/c_data.flx.values
#         ferr = r_data.ferr.values/c_data.flx.values
#
#     elif tag == "mocks":
#         r_file = r_filepaths[i]
#         r_data = pd.read_csv(r_file, delim_whitespace=True,compression='gzip',
#                             skiprows = 1,usecols=(0,1,2,4),
#                             names = ("wav", "flx", "ferr", "dloglam"))
#         wave = r_data.wav.values
#         flux = r_data.flx.values
#         ferr = r_data.ferr.values
#
#     oamin = LYA_MIN*(1+z_array[i])
#     oamax = LYA_MAX*(1+z_array[i])
#     obmin = LYB_MIN*(1+z_array[i]) #observed-beta-min
#     obmax = LYB_MAX*(1+z_array[i]) #observed-beta-max etc.
#     mask_af = (wave>oamin)&(wave<oamax)&(flux>opt.min_trans)
#     mask_bf = (wave>obmin)&(wave<obmax)&(flux>opt.min_trans)
#     wave_af = wave[mask_af] #alpha-forest (af) pixels
#     flux_af = flux[mask_af]
#
#     if "noB" in corr_tag.split('_'):
#         #noB = no beta. I create my own lyb forest.
#         oabmax = obmax*LYA_REST/LYB_REST # lyb max converted to lya wavelength
#         mask_abf = (wave>obmin)&(wave<oabmax)&(flux>opt.min_trans) # this gives larger lya forest since oabmax>oamax
#         wave_abf = wave[mask_abf] # alpha wavelength (large mask)
#         tmp_wave_abf = wave_abf*LYB_REST/LYA_REST # beta wavelength but incorrect grid
#         flux_abf = flux[mask_abf]  #alpha fluxes (larger mask)
#         tmp_flux_abf = flux_abf**(opt.f_osc_lyb/opt.f_osc_lya*LYB_REST/LYA_REST) # alpha flux (incorrect grid)
#         wave_bf = wave[mask_bf] # beta wavelength (correct grid)
#         flux_bf = griddata(points=tmp_wave_abf, values=tmp_flux_abf, xi=wave_bf, method='linear') #now fluxes are on correct grid
#         flux_bf = flux_bf*flux[mask_bf] # F_tot = F_a * F_b = e^(-tau_a-tau_b)
#
#     else:
#         wave_bf = wave[mask_bf] #beta-forest (bf) pixels
#         flux_bf = flux[mask_bf]
#     ferr_af = ferr[mask_af]
#     ferr_bf = ferr[mask_bf]
#     dloglam_af = r_data.dloglam.values[mask_af]
#     dloglam_bf = r_data.dloglam.values[mask_bf]
#     z_af = wave_af/LYA_REST - 1 #spatial grid to redshift grid from
#     z_bf = wave_bf/LYB_REST - 1 #lyman alpha and lyman beta transisions
#
#     for j in range(zbinlen): #cutting up data into redshift bins
#         mask_za = (z_af>zbin_edges[j])&(z_af<=zbin_edges[j+1])
#         mask_zb = (z_bf>zbin_edges[j])&(z_bf<=zbin_edges[j+1])
#         flux_za = flux_af[mask_za]
#         flux_zb = flux_bf[mask_zb]
#         var_za = ferr_af[mask_za]**2
#         var_zb = ferr_bf[mask_zb]**2
#         dloglam_za = dloglam_af[mask_za]
#         dloglam_zb = dloglam_bf[mask_zb]
#
#         if len(flux_za)>opt.min_pix: #we only use data if >100 pixels are in a zbin
#             bin_aqsonum[j]=bin_aqsonum[j]+1 #counting quasars
#             bin_aqsoname[j].append(qname_array[i])
#             bin_aflx[j].append(flux_za)
#             bin_avar[j].append(var_za)
#             bin_alam[j].append(dloglam_za)
#
#         if len(flux_zb)>opt.min_pix:
#             bin_bqsonum[j]=bin_bqsonum[j]+1
#             bin_bqsoname[j].append(qname_array[i])
#             bin_bflx[j].append(flux_zb)
#             bin_bvar[j].append(var_zb)
#             bin_blam[j].append(dloglam_zb)
#
#     opt.updt(nqso, i)    # Loading Bar. See options.py for more details.
# opt.updt(nqso, nqso)
# print("Complete!")
#
# # Calculating mean flux, flux variance, and noise power in each redshift bin
# mf_list_lya = []
# mf_list_tot = []
# var_list_lya = []
# var_list_tot = []
# Pn_list_lya = []
# Pn_list_tot = []
# mlam_list_lya = []
# mlam_list_tot = []
# FLUX_PDF=[]
#
# pdf_bins = np.arange(-0.025,1.05,.05)
#
# #First find dloglambda in each redshift bin.
# for i in range(zbinlen):
#     dlam_mean = np.median(np.concatenate(bin_alam[i]))
#     mlam_list_lya.append(dlam_mean)
#     try:
#         dlam_mean = np.median(np.concatenate(bin_blam[i]))
#         mlam_list_tot.append(dlam_mean)
#     except:
#         mlam_list_tot.append(np.nan)
#
# #Then find velocity separation, k_nyq,mf,nvar,npow in each z bin
# #for lya and lyb forest.
# for i in range(zbinlen):
#     dv_lya = opt.c_kms * mlam_list_lya[i] * np.log(10)
#     dv_tot = opt.c_kms * mlam_list_tot[i] * np.log(10)
#
#     k_nyquist_lya = np.pi/dv_lya
#     k_nyquist_tot = np.pi/dv_tot
#
#     ### FLUX PDF
#     #print(bin_aflx[i])
#     tmp = np.concatenate(bin_aflx[i])
#     #np.histogram(tmp)
#     FLUX_PDF.append(np.histogram(tmp,bins=pdf_bins)[0]/len(tmp))
#     ###
#
#
#     mf = np.mean(np.concatenate(bin_aflx[i]))
#     nvar = np.mean(np.concatenate(bin_avar[i])) # THIS IS <n^2>
#     npow = mf**(-2) * nvar * np.pi / k_nyquist_lya
#     mf_list_lya.append(mf)
#     var_list_lya.append(nvar)
#     Pn_list_lya.append(npow)
#
#     try:
#         mf = np.nanmean(np.concatenate(bin_bflx[i])) # process fails if uses mean
#         nvar = np.mean(np.concatenate(bin_bvar[i]))
#         npow = mf**(-2) * nvar * np.pi / k_nyquist_tot
#         mf_list_tot.append(mf)
#         var_list_tot.append(nvar)
#         Pn_list_tot.append(npow)
#     except:
#         mf_list_tot.append(np.nan)
#         var_list_tot.append(np.nan)
#         Pn_list_tot.append(np.nan)
#
# # Routine for writing results to data files for further use
# mean_flux_table = pd.DataFrame({"z":zbin_centers,
#                                 "mf_lya":mf_list_lya,  "mf_tot":mf_list_tot,
#                                 "vf_lya":var_list_lya, "vf_tot":var_list_tot,
#                                 "pn_lya":Pn_list_lya,  "pn_tot":Pn_list_tot})
#
# #Finding mean flux of just lyman beta. Requires statistical removal of lyman alpha.
# zab_centers = opt.find_za(zbin_centers) #converting lyb zbins to the equivalent, lower, lya zbins
# len_zab = len(zab_centers)
#
# #Gives corresponding lya bin for each lyb bin. organized by increasing z.
# bin_zab=np.ones(len_zab)*np.nan
# for i in range(len_zab):
#     for j in range(len_zab):
#         if (zab_centers[i]>zbin_edges[j])&(zab_centers[i]<zbin_edges[j+1]):
#             bin_zab[i] = (zbin_centers[j])
#
# mf_lyb = np.ones(len_zab)*np.nan #nan until proven otherwise
# var_lyb = np.ones(len_zab)*np.nan
# pow_lyb = np.ones(len_zab)*np.nan
# for i in range(len_zab):
#     if bin_zab[i] in mean_flux_table.z.values:
#         mask_mfab = mean_flux_table.z == mean_flux_table.z[i]
#         mask_mfa = mean_flux_table.z == bin_zab[i]
#
#         mf_tot = mean_flux_table.mf_tot[mask_mfab].values[0]
#         mf_lya = mean_flux_table.mf_lya[mask_mfa].values[0]
#
#         err_tot = mean_flux_table.vf_tot[mask_mfab].values[0]
#         err_lya = mean_flux_table.vf_lya[mask_mfa].values[0]
#
#         dv_lyb = opt.c_kms * mlam_list_tot[i] * np.log(10)
#         k_nyquist_lyb = np.pi/dv_lyb
#
#         mf_lyb[i] = mf_tot/mf_lya
#         var_lyb[i] = mf_lyb[i]**2*(err_tot**2+err_lya**2) #error propagation??
#         pow_lyb[i] = mf_lyb[i]**(-2) * var_lyb[i] * np.pi / k_nyquist_lyb # are my errors correct??
#
# # Adding results to our already existing table
# mean_flux_table["mf_lyb"] = mf_lyb
# mean_flux_table["vf_lyb"] = var_lyb
# mean_flux_table["pn_lyb"] = pow_lyb
#
#
# if verbose:
#     print("\nResults(Mean Flux Table):")
#     print(mean_flux_table)
#
# qsos_zbin_table = pd.DataFrame({"zbin":zbin_centers,"hist_lya":bin_aqsonum,"hist_lyb":bin_bqsonum})
#
# pdf_table = pd.DataFrame(np.column_stack(FLUX_PDF),columns=zbin_centers)
# pdf_table.to_csv("../output/initial_pdf_n5000.txt",index=False)
#
#
# # Saving data here
# if save_data:
#     mean_flux_table.to_csv(saving_mf_here,index=False)
#     qsos_zbin_table.to_csv(saving_qso_dist_here)
#     #pdf_table.to_csv("../output/initial_pdf_n5000.txt",index=False)

#!/usr/bin/env python

"""
Bayu Wilson, 2018-2019, University of Washington Astronomy

Usage
    ./findpk.py [data origin] [data type] [saving data?]
    ./findpk.py [obs,mocks] [corr_NR,uncorr,corrN,corrR] [0,1]

    ./findpk.py obs uncorr 1
    ./findpk.py obs corrNR 1
    ./findpk.py mocks noB_onlyA 1
    ./findpk.py mocks lyb_nocorr 1
    ./findpk.py mocks lyb_nocorr_n5000 0
Purpose
    Calculate the lya and lyb power spectrum as well as lya-lyb cross power spectrum in k,z bins

NOTES for Matt & Vid:
     1) When saving the tables, my code expects there to exist a directory called,
     `../output/`
     2) My code thinks data is located here: `../data/mocks/*` or `../data/obs/*`
     3) Many of the functions are located in a file called options.py
     4) This code is in Python 3 (not Python 2).
     5) This code depends on the output of `findmf.py`. For examply if you are using
     `lyb_nocorr` tag then you must have already run it and saved it with `findmf.py`.
     The outputted file would be `../output/mf_mocks_lyb_nocorr.csv`. After running this
     program, the power spectrum will be here: `../output/pk_mocks_lyb_nocorr.csv`.
"""

# Importing Modules
import numpy as np
import pandas as pd
import options as opt
import sys
import glob
from scipy.interpolate import griddata

# User inputted arguements
tag = str(sys.argv[1]) #obs,mocks
corr_tag = str(sys.argv[2]) #lyb_nocorr,wN,wNR,wR,wN_n5000 etc.
save_data = bool(float(sys.argv[3])) #0,1
# split_tag = corr_tag.split('_')

verbose = True

#UI Assistance
if tag=="mocks": #these are to be printed out for UI
    data_path = "../data/mocks/XQ-100_{0}/*".format(corr_tag)
if tag=="obs":
    data_path = "../data/obs/XQ-100/released/*"

saving_here = "../output/pk_{0}_{1}.csv".format(tag,corr_tag)
mf_table_path = "../output/mf_{0}_{1}.csv".format(tag,corr_tag)

print("------------------ Power Spectrum ------------------")
print("Data path:            ", data_path)
print("pk path:              ", saving_here)
print("mf path               ", mf_table_path)
#print("QSO distribution path:", saving_qso_dist_here)
print("saving data:          ", save_data)
print("Verbose:              ", verbose)

# Loading in related constants, outputs, contraints etc.
LYA_REST = opt.lya_rest
LYA_MIN = opt.lya_min
LYA_MAX = opt.lya_max
LYB_REST = opt.lyb_rest
LYB_MIN = opt.lyb_min
LYB_MAX = opt.lyb_max
zbin_centers = opt.zbin_centers
zbin_edges = opt.zbin_edges
zbinlen = len(zbin_centers)
kbin_centers = opt.kbin_centers
kbin_edges = opt.kbin_edges
kbinlen = len(kbin_centers)

# Loading table with quasar names and redshifts
### OBS ###
if tag == "obs":
    path_to_cat = "../data/obs/XQ-100_catalogue.txt"
    catalog_table = pd.read_csv(path_to_cat,delim_whitespace=True,usecols=(3,6),
                                names = ("qso_name","redshift"))
    qname_array,z_array = catalog_table['qso_name'].values,catalog_table['redshift'].values
    nqso = 100 #it will always be 100 quasars fro XQ-100 dataset

### MOCKS ###
if tag == "mocks":
    which_catalog = int(corr_tag.endswith("n5000")) # index for `path_to_cat`
    path_to_cat = ["../data/mocks/XQ-100_catalogue_n100.mock",
                   "../data/mocks/XQ-100_catalogue_n5000.mock"][which_catalog]
    catalog_table = pd.read_csv(path_to_cat,delim_whitespace=True,usecols=(3,6),
                                names = ("qso_name","redshift"))
    qname_array,z_array = catalog_table['qso_name'].values,catalog_table['redshift'].values
    nqso = int(path_to_cat.split("_n")[1][:-5]) # cutting of _n and .mock gives 100 or 5000

# print(opt.load_catalog(path_to_cat))
# if tag == "mocks":
#     r_filepaths = glob.glob("../data/{0}/XQ-100_{1}/*/ref/spectra/*.txt".format(tag,corr_tag))
#     r_filepaths.sort()

def find_configuration_data(name,zqso):
    """
    Input quasar name and quasar redshift
    Output relative flux fluctuation, pixel redshifts, normalized flux, and dloglambda
    in LYA FOREST ONLY & LYB FOREST
    """
    if tag == "obs":
        r_file = "../data/obs/XQ-100/released/%s_uvb-vis.txt"%(name)
        c_file = "../data/obs/XQ-100/continuum/%s_cont.txt"%(name)
        r_data = pd.read_csv(r_file, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
                                        names = ("wav", "flx", "ferr", "res","dloglam"))
        c_data = pd.read_csv(c_file, delim_whitespace=True, skiprows = 1,usecols=[0,1], names = ("wav",'flx'))
        wave = r_data.wav.values
        flux = r_data.flx.values/c_data.flx.values
        ferr = r_data.ferr.values/c_data.flx.values

    elif tag == "mocks":
        r_file = glob.glob("../data/mocks/XQ-100_{0}/*/ref/spectra/{1}*".format(corr_tag,name))[0]
        #"../data/mocks/XQ-100_{0}/lyb/ref/spectra/{1}_xq{2}_{0}.txt".format(corr_tag,name,nqso)
        r_data = pd.read_csv(r_file, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
                            compression='gzip',names = ("wav", "flx", "ferr", "res","dloglam"))
        wave = r_data.wav.values
        flux = r_data.flx.values
        ferr = r_data.ferr.values

    # res = r_data.res.values

    oamin = LYA_MIN*(1+zqso)
    oamax = LYA_MAX*(1+zqso)
    obmin = LYB_MIN*(1+zqso) #observed-beta-min
    obmax = LYB_MAX*(1+zqso) #observed-beta-max etc.
    mask_af = (wave>oamin)&(wave<oamax)&(flux>opt.min_trans)
    mask_bf = (wave>obmin)&(wave<obmax)&(flux>opt.min_trans)
    wave_af = wave[mask_af] #alpha-forest (af) pixels
    flux_af = flux[mask_af]

    if "noB" in corr_tag.split('_'):
        # print("!")
        #noB = no beta. I create my own lyb forest.
        oabmax = obmax*LYA_REST/LYB_REST*1.1 # lyb max converted to lya wavelength
        mask_abf = (wave>obmin)&(wave<oabmax)&(flux>opt.min_trans) # this gives larger lya forest since oabmax>oamax
        wave_abf = wave[mask_abf] # alpha wavelength (large mask)
        tmp_wave_abf = wave_abf*LYB_REST/LYA_REST # beta wavelength but incorrect grid
        flux_abf = flux[mask_abf]  #alpha fluxes (larger mask)
        tmp_flux_abf = flux_abf**(opt.f_osc_lyb/opt.f_osc_lya*LYB_REST/LYA_REST) # alpha flux (incorrect grid)
        wave_bf = wave[mask_bf] # beta wavelength (correct grid)
        flux_bf = griddata(points=tmp_wave_abf, values=tmp_flux_abf, xi=wave_bf, method='linear') #now fluxes are on correct grid
        flux_bf = flux_bf*flux[mask_bf] # F_tot = F_a * F_b = e^(-tau_a-tau_b)
        # import matplotlib.pyplot as plt TESTING REMOVE
        # plt.plot(wave_bf,flux_bf)
        # #print(min(wave_bf),max(wave_bf))
        # #print(flux_bf[0],flux_bf[-20])
        # #plt.xlim(4480,4490)
        # plt.savefig("/Users/bayuwilson/Desktop/test1.pdf")
        # plt.show()
        # #sys.exit()
        # #print("!")

    else:
        wave_bf = wave[mask_bf] #beta-forest (bf) pixels
        flux_bf = flux[mask_bf]
        # import matplotlib.pyplot as plt
        # plt.plot(wave_bf,flux_bf)
        # #print(min(wave_bf),max(wave_bf)) TESTING REMOVE
        # #print(flux_bf[0],flux_bf[-20])
        # #plt.xlim(4480,4490)
        # plt.savefig("/Users/bayuwilson/Desktop/test2.pdf")
        # plt.show()
        # #sys.exit()

    ferr_af = ferr[mask_af]
    ferr_bf = ferr[mask_bf]
    dloglam_af = r_data.dloglam.values[mask_af]
    dloglam_bf = r_data.dloglam.values[mask_bf]
    res_af = r_data.res.values[mask_af]
    res_bf = r_data.res.values[mask_bf]
    z_af = wave_af/LYA_REST - 1 #spatial grid to redshift grid from
    z_bf = wave_bf/LYB_REST - 1 #lyman alpha and lyman beta transisions

    # print(min(z_af),max(z_af)) TESTING REMOVE
    # print(min(z_bf),max(z_bf))
    # sys.exit()

    return z_af,z_bf,flux_af,flux_bf,res_af,res_bf,dloglam_af,dloglam_bf

# print(catalog_table)
def create_configuration_object():
    """
    No input
    Output pandas table containing qso name, redshift, flux fluctuation array,
    pixel redshift array, normalized flux array, resolution array, and
    dloglambda array
    """
    namez_cat = opt.load_catalog(path_to_cat) # (name, zqso) #ohsdafojasdlkfpasodkaosfjhgpsdhfouvyjdnsfhuyvg hsjdnc
    fz_pix = [find_configuration_data(namez_cat[0][x],namez_cat[1][x])
                                     for x in range(nqso)]
    return pd.DataFrame({'name_qso':namez_cat[0],
                  'zqso':namez_cat[1],
                  'z_af':[fz_pix[i][0] for i in range(nqso)],
                  'z_bf':[fz_pix[i][1] for i in range(nqso)],
                  'flux_af':[fz_pix[i][2] for i in range(nqso)],
                  'flux_bf':[fz_pix[i][3] for i in range(nqso)],
                  'res_af':[fz_pix[i][4] for i in range(nqso)],
                  'res_bf':[fz_pix[i][5] for i in range(nqso)],
                  'dloglam_af':[fz_pix[i][6] for i in range(nqso)],
                  'dloglam_bf':[fz_pix[i][7] for i in range(nqso)]},
                  columns=['name_qso','zqso','z_af','z_bf','flux_af','flux_bf','res_af',
                           'res_bf','dloglam_af','dloglam_bf'])

mf_table = pd.read_csv(mf_table_path)
pixel_catalog = create_configuration_object()

# power spectrum calculation in each k,z bin
pzaa_table = [] #lya power
pzbb_table = [] #lyb power
pzab_table = [] #lya-lyb cross power [real]
qzab_table = [] #lya-lyb cross power [imaginary]

# NEW
rfa_mat = [[] for x in range(zbinlen)] #relative flucuations lya
rfb_mat = [[] for x in range(zbinlen)] #relative flucuations lyb

load_bar = 0
#cyclying through each redshift bin
for zindex in range(zbinlen):
    # reading in measurements
    pkaa_mat = [[] for x in range(kbinlen)]
    pkbb_mat = [[] for x in range(kbinlen)]
    pkab_mat = [[] for x in range(kbinlen)]
    qkab_mat = [[] for x in range(kbinlen)]

    #cyclying through each quasar
    for qindex in range(nqso):
        load_bar=load_bar+1
        opt.updt(zbinlen*nqso,load_bar)
        z_af = pixel_catalog.z_af[qindex]
        z_bf = pixel_catalog.z_bf[qindex]
        flux_af = pixel_catalog.flux_af[qindex]
        flux_bf = pixel_catalog.flux_bf[qindex]
        res_af = pixel_catalog.res_af[qindex]
        res_bf = pixel_catalog.res_bf[qindex]
        dloglam_af = pixel_catalog.dloglam_af[qindex]
        dloglam_bf = pixel_catalog.dloglam_bf[qindex]

        zmask_af = (z_af>zbin_edges[zindex])&(z_af<=zbin_edges[zindex+1])
        zmask_bf = (z_bf>zbin_edges[zindex])&(z_bf<=zbin_edges[zindex+1])

        rf_af = flux_af/mf_table.mf_lya[zindex] - 1
        rf_bf = flux_bf/mf_table.mf_tot[zindex] - 1 # IS THIS THE RIGHT THING TO DO

        # print(rf_af[zmask_af])
        #print(len(rf_af),len(zmask_af),len(flux_af), len(flux_af/mf_table.mf_lya[zindex]))
        if len(rf_af[zmask_af])>opt.min_pix:
            rfa_mat[zindex].append(rf_af[zmask_af])

            pk_table = opt.power_spectrum_fft(rf_af[zmask_af],dloglam_af[zmask_af])
            k_af,pk_af = pk_table['k'],pk_table['Pk']

            for kindex in range(kbinlen):
                kmask = (k_af>kbin_edges[kindex])&(k_af<=kbin_edges[kindex+1])&(k_af!=0)
                ksubset = k_af[kmask]
                pk_sub = pk_af[kmask]
                dv_array = opt.c_kms*np.log(10)*dloglam_af[zmask_af][kmask]
                res_wave = res_af[zmask_af][kmask]
                # print(split_tag)
                # sys.exit(':(')
                if (('corrNR' in corr_tag) or
                    ('corrN' in corr_tag) or
                    ('wN' in corr_tag) or
                    ('wNR' in corr_tag)):
                    pk_sub = pk_sub - mf_table.pn_lya[zindex]
                if (('corrNR' in corr_tag) or
                    ('corrR' in corr_tag) or
                    ('wR' in corr_tag) or
                    ('wNR' in corr_tag)):
                    pk_sub = pk_sub / opt.window(k=ksubset,p=dv_array,R=res_wave)**2
                pkaa_mat[kindex].append(pk_sub)

        if len(rf_bf[zmask_bf])>opt.min_pix:
            rfb_mat[zindex].append(rf_bf[zmask_bf])

            pk_table = opt.power_spectrum_fft(rf_bf[zmask_bf],dloglam_bf[zmask_bf])
            k_bf,pk_bf = pk_table['k'],pk_table['Pk']
            for kindex in range(kbinlen):
                kmask = (k_bf>kbin_edges[kindex])&(k_bf<=kbin_edges[kindex+1])&(k_bf!=0)
                ksubset = k_bf[kmask]
                pk_sub = pk_bf[kmask]
                dv_array = opt.c_kms*np.log(10)*dloglam_bf[zmask_bf][kmask]
                res_wave = res_bf[zmask_bf][kmask]
                if (('corrNR' in corr_tag) or
                    ('corrN' in corr_tag) or
                    ('wN' in corr_tag) or
                    ('wNR' in corr_tag)):
                    pk_sub = pk_sub - mf_table.pn_lya[zindex]
                    #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
                if (('corrNR' in corr_tag) or
                    ('corrR' in corr_tag) or
                    ('wR' in corr_tag) or
                    ('wNR' in corr_tag)):
                    pk_sub = pk_sub / opt.window(k=ksubset,p=dv_array,R=res_wave)**2
                    #print(kpower_subset[:2])
                    #pk_lyb = pk_lyb / window(k=k_lyb,p=dloglam_lyb,R=res_lyb)**2
                pkbb_mat[kindex].append(pk_sub)

        len_nan = len(np.where(np.isfinite(rf_bf[zmask_bf]))[0])
        if (len_nan>opt.min_pix)&(len(rf_af[zmask_af])>opt.min_pix):
            #print(rf_af[0],rf_bf[0])
            z_af,z_bf = z_af[zmask_af],z_bf[zmask_bf]
            rf_af,rf_bf = rf_af[zmask_af],rf_bf[zmask_bf]
            dloglam_af,dloglam_bf = dloglam_af[zmask_af],dloglam_bf[zmask_bf]
            res_af,res_bf = res_af[zmask_af],res_bf[zmask_bf]
            new_af_mask = (z_af>min(z_bf))&(z_af<max(z_bf))
            new_bf_mask = (z_bf>min(z_af))&(z_bf<max(z_af))

            #print(len(rf_af),len(rf_bf))
            pk_table = opt.cross_pk_fft(z_af[new_af_mask]      ,z_bf[new_bf_mask],
                                        rf_af[new_af_mask]     ,rf_bf[new_bf_mask],
                                        dloglam_af[new_af_mask],dloglam_bf[new_bf_mask],
                                        res_af[new_af_mask]    ,res_bf[new_bf_mask])

            k,Pab,Qab,dloglam,res = pk_table['k'],pk_table['P_real'],pk_table['P_imag'],pk_table['dloglam'],pk_table['res']
            dv_array = opt.c_kms*np.log(10)*dloglam

            # Noise Correction
            #if 'corrNR'or'corrN'or'wNR'or'wN' in split_tag:
            if (('corrNR' in corr_tag) or
                ('corrN' in corr_tag) or
                ('wN' in corr_tag) or
                ('wNR' in corr_tag)):
                #powerk = powerk - mean_flux_table.powerN[zindex]
                #pk_lya = pk_lya - mean_flux_table.powerN_lya[zindex]
                #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
                Pab = Pab - mf_table.pn_lya[zindex]
                Qab = Qab - mf_table.pn_lya[zindex]
                #print("!!!")

            # Beam Correction
            #if 'corrNR'or'corrR'or'wNR'or'wR' in split_tag:
            if (('corrNR' in corr_tag) or
                ('corrR' in corr_tag) or
                ('wR' in corr_tag) or
                ('wNR' in corr_tag)):
                Pab = Pab / opt.window(k=k,p=dv_array,R=res)**2
                Qab = Qab / opt.window(k=k,p=dv_array,R=res)**2

            for kindex in range(kbinlen):
                kmask = (k>kbin_edges[kindex])&(k<=kbin_edges[kindex+1])&(k!=0)
                ksubset = k[kmask]
                Pab_sub = Pab[kmask]
                Qab_sub = Qab[kmask]
                pkab_mat[kindex].append(Pab_sub)
                qkab_mat[kindex].append(Qab_sub)

    pk_means_lya = []
    pk_means_lyb = []
    pk_means_real = []
    pk_means_imag = []
    for columnindex in range(kbinlen):
        try:
            mean_lya = np.mean(np.concatenate(pkaa_mat[columnindex]))
        except:
            mean_lya = np.nan
        try:
            mean_lyb = np.nanmean(np.concatenate(pkbb_mat[columnindex]))
        except:
            mean_lyb = np.nan
        try:
            mean_real = np.mean(np.concatenate(pkab_mat[columnindex]))
        except:
            mean_real = np.nan
        try:
            mean_imag = np.mean(np.concatenate(qkab_mat[columnindex]))
        except:
            mean_imag = np.nan

        pk_means_lya.append(mean_lya)
        pk_means_lyb.append(mean_lyb)
        pk_means_real.append(mean_real)
        pk_means_imag.append(mean_imag)
    pzaa_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
           kbin_centers,np.array(pk_means_lya)])
    pzbb_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
           kbin_centers,np.array(pk_means_lyb)])
    pzab_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
           kbin_centers,np.array(pk_means_real)])
    qzab_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
           kbin_centers,np.array(pk_means_imag)])

lya_cfluc_squared = []
lyb_cfluc_squared = []
Nlya_list = []
Nlyb_list = []
for i in range(zbinlen):
    try:
        lya_cfluc_squared.append(np.mean(np.concatenate(rfa_mat[i])**2))
        Nlya_list.append(len(np.concatenate(rfa_mat[i])**2))
    except:
        lya_cfluc_squared.append(np.nan)
        Nlya_list.append(np.nan)
    try:
        lyb_cfluc_squared.append(np.mean(np.concatenate(rfb_mat[i])**2))
        Nlyb_list.append(len(np.concatenate(rfb_mat[i])**2))
    except:
        lyb_cfluc_squared.append(np.nan)
        Nlyb_list.append(np.nan)

varF_lya = mf_table.vf_lya.values
varF_lyb = mf_table.vf_lyb.values
mf_lya = mf_table.mf_lya.values
mf_lyb = mf_table.mf_lyb.values
# N_lya =
# N_lyb =

# print(Nlya_list)
# print(Nlyb_list)
#mf_table['err_mf_lya'] = np.sqrt(varF_lya+mf_lya*np.array(lya_cfluc_squared)/(np.array(Nlya_list)-1)) #divide by N-1
#mf_table['err_mf_lyb'] = np.sqrt(varF_lyb+mf_lyb*np.array(lyb_cfluc_squared)/(np.array(Nlyb_list)-1))
#print(mean_flux_table)


cols_lya = np.column_stack(pzaa_table)
cols_lyb = np.column_stack(pzbb_table)
cols_real = np.column_stack(pzab_table)
cols_imag = np.column_stack(qzab_table)
dict_everything={"z":cols_lya[0],"k":cols_lya[1],
                 "P_aa":cols_lya[2],
                 "P_tot":cols_lyb[2],
                 "P_ab":cols_real[2],
                 "Q_ab":cols_imag[2]}
final_t = pd.DataFrame(dict_everything)

zab_centers = opt.find_za(zbin_centers) #converting lyb zbins to the equivalent, lower, lya zbins
len_zab = len(zab_centers)

#Gives corresponding lya bin for each lyb bin. organized by increasing z.
bin_zab=np.ones(len_zab)*np.nan
for i in range(len_zab):
    for j in range(len_zab):
        if (zab_centers[i]>zbin_edges[j])&(zab_centers[i]<zbin_edges[j+1]):
            bin_zab[i] = (zbin_centers[j])

pbb_arr = np.array([])
for zidx in range(zbinlen):
    try:
        ptot = final_t.P_tot[final_t.z == zbin_centers[zidx]].values
        paa = final_t.P_aa[final_t.z == bin_zab[zidx]].values
        pbb = ptot-paa
        pbb_arr = np.append(pbb_arr,pbb)
    except:
        pbb_arr = np.append(pbb_arr,np.ones(kbinlen)*np.nan)

final_t['P_bb'] = pbb_arr

#
# # mf_lyb = np.ones(len_zab)*np.nan #nan until proven otherwise
# pk_lyb = np.ones(kbinlen*zbinlen)*np.nan #nan until proven otherwise
# # var_lyb = np.ones(len_zab)*np.nan
# # pow_lyb = np.ones(len_zab)*np.nan
# for i in range(kbinlen*zbinlen):
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
#         var_lyb[i] = mf_lyb[i]**2*(err_tot**2+err_lya**2)
#         # pow_lyb[i] = mf_lyb[i]**(-2) * var_lyb[i] * np.pi / k_nyquist_lyb


#saving_here = "results/pk_obs_{0}.csv".format(tag)
#print(everything_table)
# print(saving_here)

if verbose:
    print(final_t)
    #print(mf_table)
if save_data:
    final_t.to_csv(saving_here,index=False)
    #mean_flux_table.to_csv(mf_table,index=False)

# print("\n")




# z_lya_in_lyb = find_lya(np.unique(opt.zbin_centers))#table_tot.z.values
# z_binned_lyalyb=[]
# for i in z_lya_in_lyb:
#     for j in range(len(zbin_centers)):
#         if (i>zbin_edges[j])&(i<zbin_edges[j+1]):
#             z_binned_lyalyb.append(zbin_centers[j])
# z_binned_lyalyb=np.round(np.array(z_binned_lyalyb),2) #lya pixels in lyb forest
# z_lyalyb_edges = np.append(z_binned_lyalyb,3.6)-0.1



# # catalog = opt.load_catalog(path_to_cat)
# # print(catalog)
# x = find_configuration_data(qname_array[0],z_array[0])
# import matplotlib.pyplot as plt
# # plt.plot(x[0],np.ones_like(x[0]))
# plt.plot(x[1],np.ones_like(x[1]))
#
# # plt.plot(x[0],np.ones_like(x[0]))
# plt.plot(x[1],x[3])
# plt.show()

#
# if tag == "obs":
#     path_to_cat = "../Data/XQ-100_catalogue.txt"
#     nqso = 100
#     data_path_name = "../Data/XQ-100/released/*"
# if tag == "mocks":
#     which_catalog = int(tag.endswith("n5000"))
#     path_to_cat = ["XQ-100_catalogue_n100.mock","XQ-100_catalogue_n5000.mock"][which_catalog]
#     nqso = int(path_to_cat.split("_n")[1][:-5])
#     data_path_name = "XQ-100_{0}/lyb/ref/spectra/mock-*".format(corr_tag)
#
# mf_data_table = "results/mf_{0}.csv".format(tag)
# saving_here = "results/pk_{0}_{1}.csv".format(tag,corr_tag)
# # if tag=="mocks":
# #     data_path_name = "XQ-100_{0}/lyb/ref/spectra/mock-*".format(corr_tag)
# # if tag=="obs":
# #     data_path_name = "../Data/XQ-100/released/*"

# print("------------------ POWER SPECTRUM ------------------")
# print("Mean flux data:     ", mf_data_table)
# print("Sightline paths:    ", data_path_name)
# print("Saving P(k) here:   ", saving_here)
# print("Saving data:        ", save_data)
# print(" ")
#
# # Loading in related constants
# LYA_REST = opt.lya_rest
# LYA_MIN = opt.lya_min
# LYA_MAX = opt.lya_max
# LYB_REST = opt.lyb_rest
# LYB_MIN = opt.lyb_min
# LYB_MAX = opt.lyb_max
#
# def find_lya(lyb_z):
#     return (lyb_z + 1)*(LYB_REST/LYA_REST)-1
# def find_lyb(lya_z):
#     return (lya_z + 1)*(LYA_REST/LYB_REST)-1
#
# def load_catalog(path=path_to_cat):
#     """
#     Input catalog path
#     Output qso names and qso redshifts in catalog
#     """
#     catalog = pd.read_csv(path,delim_whitespace=True,usecols=(3,6),
#                           names = ("qso_name","redshift"))
#     return catalog["qso_name"].values,catalog["redshift"].values
#
# def find_configuration_data(name,zqso):
#     """
#     Input quasar name and quasar redshift
#     Output relative flux fluctuation, pixel redshifts, normalized flux, and dloglambda
#     in LYA FOREST ONLY & LYB FOREST
#     """
#     if tag == "obs":
#         FILE_r = "../Data/XQ-100/released/%s_uvb-vis.txt"%(name)
#         FILE_c = "../Data/XQ-100/continuum/%s_cont.txt"%(name)
#         released_filedata = pd.read_csv(FILE_r, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
#                                         names = ("wavelength", "nflux","ferr","resolution","dloglambda"))
#         cont = pd.read_csv(FILE_c, delim_whitespace=True, skiprows = 1,usecols=[0,1], names = ("wav",'flx'))
#         wavelength = released_filedata.wavelength.values
#         nflux = released_filedata.nflux.values/cont.flx.values
#         ferr = released_filedata.ferr.values/cont.flx.values
#
#
#     elif tag == "mocks":
#         FILE = "XQ-100_{0}/lyb/ref/spectra/{1}_xq{2}_{0}.txt".format(corr_tag,name,nqso)
#         released_filedata = pd.read_csv(FILE, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
#                                         names = ("wavelength", "nflux","ferr","resolution","dloglambda"))
#         wavelength = released_filedata.wavelength.values
#         nflux = released_filedata.nflux.values
#         ferr = released_filedata.ferr.values
#
#     #wavelength = released_filedata.wavelength.values
#     resolution = released_filedata.resolution.values#released_filedata["resolution"].values
#     #nflux = released_filedata.nflux.values/cont.flx.values
#     obs_LYB_MIN = LYB_MIN*(1+zqso)
#     obs_LYB_MAX = LYB_MAX*(1+zqso)
#     obs_LYA_MIN = LYA_MIN*(1+zqso)
#     obs_LYA_MAX = LYA_MAX*(1+zqso)
#     mask_lyb = (wavelength>obs_LYB_MIN)&(wavelength<obs_LYB_MAX)&(nflux>opt.min_trans)
#     mask_lya = (wavelength>obs_LYA_MIN)&(wavelength<obs_LYA_MAX)&(nflux>opt.min_trans)
#     wavelength_lya = wavelength[mask_lya]
#     wavelength_lyb = wavelength[mask_lyb]
#     nflux_lyb = nflux[mask_lyb]
#     nflux_lya = nflux[mask_lya]
#
#     # NEW
#     ferr_lya =  ferr[mask_lya]
#     ferr_lyb = ferr[mask_lyb]
#
#     res_lya = resolution[mask_lya]
#     res_lyb = resolution[mask_lyb]
#     dloglambda_lyb = released_filedata.dloglambda.values[mask_lyb]
#     dloglambda_lya = released_filedata.dloglambda.values[mask_lya]
#     z_lya = wavelength_lya/LYA_REST - 1
#     z_lyb = wavelength_lyb/LYB_REST - 1
#
#     return z_lya,z_lyb,nflux_lya,nflux_lyb,res_lya,res_lyb,dloglambda_lya,dloglambda_lyb
# def create_configuration_object():
#     """
#     No input
#     Output pandas table containing qso name, redshift, flux fluctuation array,
#     pixel redshift array, normalized flux array, resolution array, and
#     dloglambda array
#     """
#     namez_cat = load_catalog() # (name, zqso)
#     fz_pix = [find_configuration_data(namez_cat[0][x],namez_cat[1][x]) for x in range(nqso)]
#     return pd.DataFrame({'QSO_NAME':namez_cat[0],
#                   'QSO_z':namez_cat[1],
#                   'z_lya':[fz_pix[i][0] for i in range(nqso)],
#                   'z_lyb':[fz_pix[i][1] for i in range(nqso)],
#                   'nflux_lya':[fz_pix[i][2] for i in range(nqso)],
#                   'nflux_lyb':[fz_pix[i][3] for i in range(nqso)],
#                   'res_lya':[fz_pix[i][4] for i in range(nqso)],
#                   'res_lyb':[fz_pix[i][5] for i in range(nqso)],
#                   'dloglambda_lya':[fz_pix[i][6] for i in range(nqso)],
#                   'dloglambda_lyb':[fz_pix[i][7] for i in range(nqso)]},
#                   columns=['QSO_NAME','QSO_z','z_lya','z_lyb','nflux_lya','nflux_lyb','res_lya',
#                            'res_lyb','dloglambda_lya','dloglambda_lyb'])
#
# def power_spectrum_fft(c_fluc,dloglambda):
#     """
#     Input relative flux fluctuations in configuration space
#     Output dictionary containing wavenumber, k, and corresponding power spectrum for each k value
#     """
#     dv = opt.c_kms * dloglambda * np.log(10)
#     k_fluc = np.fft.fft(c_fluc)*dv
#     N = len(k_fluc)
#     V = N*dv
#     dk = 2*np.pi/(N*dv)
#     k = dk*np.arange(0,N,1)
#     Pk = np.abs(k_fluc)**2/V
#     dictionary = {'k':k,'Pk':Pk}
#     return dictionary
#
# def grid_interp(z1,z2,rf1,rf2):
#     #which_one=1 # 1 means changing z1 & rf1, so z1 and rf2 are the references
#     #              0 means changing z2 & rf2 so z1 and rf1 are the references
#     if (min(z1)<=min(z2))&(max(z1)>=max(z2)):
#         changing_this = 0
#         ref_grid, other_pts, other_data = z1,z2,rf2
#     elif (min(z2)<min(z1))&(max(z2)>max(z1)):
#         changing_this = 1
#         ref_grid, other_pts, other_data = z2,z1,rf1
#     elif len(z1)<=len(z2):
#         changing_this = 0
#         ref_grid, other_pts, other_data = z1,z2,rf2
#     elif len(z1)>len(z2):
#         changing_this = 1
#         ref_grid, other_pts, other_data = z2,z1,rf1
#     else:
#         changing_this = 0
#         ref_grid, other_pts, other_data = z1,z2,rf2
#
#     new_rf = griddata(other_pts,other_data,ref_grid,method='linear')
#     return ref_grid,new_rf,changing_this
#
# def cross_pk_fft(z_lya,z_lyb,rf_lya,rf_lyb,dloglam_lya,dloglam_lyb,res_lya,res_lyb):
#     #c_fluc_lya,c_fluc_lyb,dloglambda_lya,dloglambda_lyb):
#     """
#     Input relative flux fluctuations in configuration space
#     Output dictionary containing wavenumber, k, and corresponding power spectrum for each k value
#     """
#     x,y,which = grid_interp(z_lya,z_lyb,rf_lya,rf_lyb) #reference, changed, which one
#     finite_mask = np.isfinite(y)
#     #print(len(x),len(y),len(finite_mask))
#     if which == 1: #lyb is referece -- confirmed
#         z_lyb,rf_lyb = z_lyb[finite_mask],rf_lyb[finite_mask]
#         z_lya,rf_lya,dloglam,res = x[finite_mask],y[finite_mask],dloglam_lyb[finite_mask],res_lyb[finite_mask]
#     elif which == 0: #lya is reference -- confirmed
#         z_lya,rf_lya = z_lya[finite_mask],rf_lya[finite_mask]
#         z_lyb,rf_lyb,dloglam,res = x[finite_mask],y[finite_mask],dloglam_lya[finite_mask],res_lya[finite_mask]
# #     print(len(rf_lya),len(rf_lyb),len(dloglam))
#     dv = opt.c_kms * dloglam * np.log(10)
#
#     kf_lya = np.fft.fft(rf_lya)*dv #dloglam
#     kf_lyb = np.fft.fft(rf_lyb)*dv #dloglam
#     N=len(kf_lya)
#     V = N*dv
#     dk = 2*np.pi/(N*dv)
#     k=dk*np.arange(0,N,1)
#     P_cross = np.conj(kf_lya)*kf_lyb/V
#     P_real,P_imag = P_cross.real,P_cross.imag
# #     print(len(k),len(finite_mask),len(z_lyb))
#     dictionary = {'k':k,'P_real':P_real,'P_imag':P_imag,'dloglam':dloglam,'res':res}
#                   #'grid_mask':finite_mask}
#     return dictionary
#
# def window(k,p,R):
#     """
#     Deconvolution kernel.
#     k: velocity wavenumber
#     p: pixel width (dloglambda)
#     R: resolution
#     """
#     #print(p)
#     #p = p*opt.c_kms*np.log(10) #Check this part with somebody
#     gaussian = np.exp(-0.5*k**2*R**2)
#     tophat =  np.sinc(0.5*k*p/np.pi) #np.sin(k*p/2)/(k*p/2)
#     deconvolution_kernel = gaussian*tophat
#     return deconvolution_kernel
#
# # kbin_logcenters = np.logspace(np.log10(opt.kmin),np.log10(opt.kmax),9)
# # dlogk = kbin_logcenters[1]-kbin_logcenters[0]
# # kbin_logedges = np.logspace(np.log10(opt.kmin)-dlogk,np.log10(opt.kmax)+dlogk,10)
#
# mean_flux_table = pd.read_csv(mf_data_table)
# pixel_catalog = create_configuration_object()
#
# # zbin_centers = opt.zbin_centers
# # zbin_edges = opt.zbin_edges
# # zbinlen = len(zbin_centers)
# zbin_centers = opt.zbin_centers
# zbin_edges = opt.zbin_edges
# zbinlen = opt.zbinlen
# # kbinlen = 11
# # dk_test = (opt.kmax-opt.kmin)/kbinlen
# # kbin_centers = np.linspace(opt.kmin,opt.kmax,kbinlen) #kbin_logcenters #np.arange(opt.kmin,opt.kmax,opt.dlogk)
# # kbin_edges = np.linspace(opt.kmin-dk_test,opt.kmax+dk_test,kbinlen+1) #kbin_logedges #opt.kbin_edges
# # kbinlen = len(kbin_centers)
# kbin_centers = opt.kbin_centers
# kbin_edges = opt.kbin_edges
# kbinlen = len(kbin_centers)
#
#
# test = 0
# # power spectrum calculation in each k,z bin
# pk_table_lya = []
# pzbb_table = []
# pzab_table = []
# qzab_table = []
#
# # NEW
# rfa_mat = [[] for x in range(zbinlen)] # squared
# rfb_mat = [[] for x in range(zbinlen)] #
#
# for zindex in range(zbinlen):
#     # reading in measurements
#     pkaa_mat = [[] for x in range(kbinlen)]
#     pkbb_mat = [[] for x in range(kbinlen)]
#     pkab_mat  = [[] for x in range(kbinlen)]
#     qkab_mat  = [[] for x in range(kbinlen)]
#
#     for qindex in range(nqso):
#         test=test+1
#         opt.updt(zbinlen*nqso,test)
#         z_lya = pixel_catalog.z_lya[qindex]
#         z_lyb = pixel_catalog.z_lyb[qindex]
#         #zb_min,zb_max = np.min(z_lyb),np.max(z_lyb)
#         nflux_lya = pixel_catalog.nflux_lya[qindex]
#         nflux_lyb = pixel_catalog.nflux_lyb[qindex]
#         res_lya = pixel_catalog.res_lya[qindex]
#         res_lyb = pixel_catalog.res_lyb[qindex]
#         dloglambda_lya = pixel_catalog.dloglambda_lya[qindex]
#         dloglambda_lyb = pixel_catalog.dloglambda_lyb[qindex]
#
#         zmask_lya = (z_lya>zbin_edges[zindex])&(z_lya<=zbin_edges[zindex+1])
#         zmask_lyb = (z_lyb>zbin_edges[zindex])&(z_lyb<=zbin_edges[zindex+1])
#         #zmask_lya_subset = (z_lya>min(z_lyb))#&(z_lya<=zb_max) #[x1>min(x2)]
#         #zmask_lyb_subset = (z_lyb<max(z_lya))#&(z_lya<=zb_max) #[x1>min(x2)]
#         #print(len(zlya[zmask_lya_subset]),len(z_lyb[zmask_lyb]), qindex,zindex)
#
#         rfluc_lya = nflux_lya/mean_flux_table.meanF_lya[zindex] - 1
#         rfluc_lyb = nflux_lyb/mean_flux_table.meanF_lyb[zindex] - 1
#
#         if len(rfluc_lya[zmask_lya])>opt.min_pix:
#             rfa_mat[zindex].append(rfluc_lya[zmask_lya])
#
#             pk_ktable = power_spectrum_fft(rfluc_lya[zmask_lya],dloglambda_lya[zmask_lya])
#             k_lya,pk_lya = pk_ktable['k'],pk_ktable['Pk']
#             for kindex in range(kbinlen):
#                 kmask = (k_lya>kbin_edges[kindex])&(k_lya<=kbin_edges[kindex+1])&(k_lya!=0)
#                 ksubset = k_lya[kmask]
#                 kpower_subset = pk_lya[kmask]
#                 p_wave = dloglambda_lya[zmask_lya][kmask]
#                 res_wave = res_lya[zmask_lya][kmask]
#                 if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrN'):
#                     #print("{0:.3e} {1:.3e}".format(kpower_subset[2],mean_flux_table.powerN_lya[zindex]))
#                     kpower_subset = kpower_subset - mean_flux_table.powerN_lya[zindex]
#
#                     #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
#                 if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrR'):
#                     kpower_subset = kpower_subset / window(k=ksubset,p=p_wave,R=res_wave)**2
#                     #print(window(k=ksubset,p=p_wave,R=res_wave)**2)
#                     #pk_lyb = pk_lyb / window(k=k_lyb,p=dloglam_lyb,R=res_lyb)**2
#                 pkaa_mat[kindex].append(kpower_subset)
#                 #print(window(k=ksubset,p=p_wave,R=res_wave))
#
#         if len(rfluc_lyb[zmask_lyb])>opt.min_pix:
#             rfb_mat[zindex].append(rfluc_lyb[zmask_lyb])
#
#             pk_ktable = power_spectrum_fft(rfluc_lyb[zmask_lyb],dloglambda_lyb[zmask_lyb])
#             k_lyb,pk_lyb = pk_ktable['k'],pk_ktable['Pk']
#             for kindex in range(kbinlen):
#                 kmask = (k_lyb>kbin_edges[kindex])&(k_lyb<=kbin_edges[kindex+1])&(k_lyb!=0)
#                 ksubset = k_lyb[kmask]
#                 kpower_subset = pk_lyb[kmask]
#                 p_wave = dloglambda_lyb[zmask_lyb][kmask]
#                 res_wave = res_lyb[zmask_lyb][kmask]
#                 if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrN'):
#                     kpower_subset = kpower_subset - mean_flux_table.powerN_lya[zindex]
#                     #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
#                 if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrR'):
#                     kpower_subset = kpower_subset / window(k=ksubset,p=p_wave,R=res_wave)**2
#                     #print(kpower_subset[:2])
#                     #pk_lyb = pk_lyb / window(k=k_lyb,p=dloglam_lyb,R=res_lyb)**2
#                 pkbb_mat[kindex].append(kpower_subset)
#
#         if (len(rfluc_lyb[zmask_lyb])>opt.min_pix)&(len(rfluc_lya[zmask_lya])>opt.min_pix):
#             #print(len(rfluc_lya[zmask_lya]),len(rfluc_lyb[zmask_lyb]))
#             #print(len(rfluc_lya[zmask_lya]),len(rfluc_lyb[zmask_lyb]))
#             #k_lya,k_lyb = k_lya[zmask_lya],k_lyb[zmask_lyb]
#             z_lya,z_lyb = z_lya[zmask_lya],z_lyb[zmask_lyb]
#             rf_lya,rf_lyb = rfluc_lya[zmask_lya],rfluc_lyb[zmask_lyb]
#             dloglam_lya,dloglam_lyb = dloglambda_lya[zmask_lya],dloglambda_lyb[zmask_lyb]
#             res_lya,res_lyb = res_lya[zmask_lya],res_lyb[zmask_lyb]
#             new_lya_mask = (z_lya>min(z_lyb))&(z_lya<max(z_lyb))
#             new_lyb_mask = (z_lyb>min(z_lya))&(z_lyb<max(z_lya))
#             #z1,z2,rf1,rf2 = z1[mask1],z2[mask2],rf1[mask1],rf2[mask2]
#             pk_ktable = cross_pk_fft(z_lya[new_lya_mask],z_lyb[new_lyb_mask],
#                                      rf_lya[new_lya_mask],rf_lyb[new_lyb_mask],
#                                      dloglam_lya[new_lya_mask],dloglam_lyb[new_lyb_mask],
#                                      res_lya[new_lya_mask],res_lyb[new_lyb_mask])
#             #k,P_real,P_imag,grid_mask = pk_ktable['k'],pk_ktable['P_real'],pk_ktable['P_imag'],pk_ktable['grid_mask']
#             #k,P_real,P_imag
#             k,P_real,P_imag,dloglam,res = pk_ktable['k'],pk_ktable['P_real'],pk_ktable['P_imag'],pk_ktable['dloglam'],pk_ktable['res']
#
#             ### TEST HERE
#             # Noise Correction
#             if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrN'):
#                 #powerk = powerk - mean_flux_table.powerN[zindex]
#                 #pk_lya = pk_lya - mean_flux_table.powerN_lya[zindex]
#                 #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
#                 P_real = P_real - mean_flux_table.powerN_lya[zindex]
#                 P_imag = P_imag - mean_flux_table.powerN_lya[zindex]
#                 #print("!!!")
#
#             # Beam Correction
#             if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrR'):
#                 # print(P_real[10:13])
#                 #print(res_lya[25:30])
#                 #pk_lya = pk_lya / window(k=k_lya,p=dloglam_lya,R=res_lya)**2
#                 #pk_lyb = pk_lyb / window(k=k_lyb,p=dloglam_lyb,R=res_lyb)**2
#                 #print(len(k),len(dloglam_lya[new_lya_mask]),len(res_lya[new_lya_mask]),len(P_real))
#                 P_real = P_real / window(k=k,p=dloglam,R=res)**2
#                 P_imag = P_imag / window(k=k,p=dloglam,R=res)**2
#                 # print(P_real[10:13])
#                 #print((window(k=k,p=dloglam,R=res)**2)[:3])
#                 #print("000")
#                 #print(res_lya[25:30])
#
#             ### TEST HERE
#             for kindex in range(kbinlen):
#                 kmask = (k>kbin_edges[kindex])&(k<=kbin_edges[kindex+1])&(k!=0)
#                 ksubset = k[kmask]
#                 pk_sub_real = P_real[kmask]
#                 pk_sub_imag = P_imag[kmask]
#                 pkab_mat[kindex].append(pk_sub_real)
#                 qkab_mat[kindex].append(pk_sub_imag)
#
#     pk_means_lya = []
#     pk_means_lyb = []
#     pk_means_real = []
#     pk_means_imag = []
#     for columnindex in range(kbinlen):
#         try:
#             mean_lya = np.mean(np.concatenate(pkaa_mat[columnindex]))
#         except:
#             mean_lya=np.nan
#         try:
#             mean_lyb = np.mean(np.concatenate(pkbb_mat[columnindex]))
#         except:
#             mean_lyb = np.nan
#         try:
#             mean_real = np.mean(np.concatenate(pkab_mat[columnindex]))
#         except:
#             mean_real = np.nan
#         try:
#             mean_imag = np.mean(np.concatenate(qkab_mat[columnindex]))
#         except:
#             mean_imag = np.nan
#
#         pk_means_lya.append(mean_lya)
#         pk_means_lyb.append(mean_lyb)
#         pk_means_real.append(mean_real)
#         pk_means_imag.append(mean_imag)
#     pk_table_lya.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
#            kbin_centers,np.array(pk_means_lya)])
#     pzbb_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
#            kbin_centers,np.array(pk_means_lyb)])
#     pzab_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
#            kbin_centers,np.array(pk_means_real)])
#     qzab_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
#            kbin_centers,np.array(pk_means_imag)])
#
# lya_cfluc_squared = []
# lyb_cfluc_squared = []
# Nlya_list = []
# Nlyb_list = []
# for i in range(zbinlen):
#     try:
#         lya_cfluc_squared.append(np.mean(np.concatenate(rfa_mat[i])**2))
#         Nlya_list.append(len(np.concatenate(rfa_mat[i])**2))
#     except:
#         lya_cfluc_squared.append(np.nan)
#         Nlya_list.append(np.nan)
#     try:
#         lyb_cfluc_squared.append(np.mean(np.concatenate(rfb_mat[i])**2))
#         Nlyb_list.append(len(np.concatenate(rfb_mat[i])**2))
#     except:
#         lyb_cfluc_squared.append(np.nan)
#         Nlyb_list.append(np.nan)
#
# varF_lya = mean_flux_table.varF_lya.values
# varF_lyb = mean_flux_table.varF_lyb.values
# mf_lya = mean_flux_table.meanF_lya.values
# mf_lyb =mean_flux_table.meanF_lyb.values
# # N_lya =
# # N_lyb =
#
# # print(Nlya_list)
# # print(Nlyb_list)
# mean_flux_table['err_mf_lya'] = np.sqrt(varF_lya+mf_lya*np.array(lya_cfluc_squared)/(np.array(Nlya_list)-1)) #divide by N-1
# mean_flux_table['err_mf_lyb'] = np.sqrt(varF_lyb+mf_lyb*np.array(lyb_cfluc_squared)/(np.array(Nlyb_list)-1))
# #print(mean_flux_table)
#
#
# cols_lya = np.column_stack(pk_table_lya)
# cols_lyb = np.column_stack(pzbb_table)
# cols_real = np.column_stack(pzab_table)
# cols_imag = np.column_stack(qzab_table)
# dict_everything={"z":cols_lya[0],"k":cols_lya[1],
#                  "P_aa":cols_lya[2],
#                  "P_bb":cols_lyb[2],
#                  "P_ab":cols_real[2],
#                  "Q_ab":cols_imag[2]}
# everything_table = pd.DataFrame(dict_everything)
# #saving_here = "results/pk_obs_{0}.csv".format(tag)
# #print(everything_table)
# # print(saving_here)
# if save_data:
#     everything_table.to_csv(saving_here,index=False)
#     mean_flux_table.to_csv(mf_data_table,index=False)
#
# print("\n")
#
#
# # z_lya_in_lyb = find_lya(np.unique(opt.zbin_centers))#table_tot.z.values
# # z_binned_lyalyb=[]
# # for i in z_lya_in_lyb:
# #     for j in range(len(zbin_centers)):
# #         if (i>zbin_edges[j])&(i<zbin_edges[j+1]):
# #             z_binned_lyalyb.append(zbin_centers[j])
# # z_binned_lyalyb=np.round(np.array(z_binned_lyalyb),2) #lya pixels in lyb forest
# # z_lyalyb_edges = np.append(z_binned_lyalyb,3.6)-0.1

#!/usr/bin/env python
"""
BOOTSTRAPPING ERROR BARS with bootstrap_lyb.py
Bayu Wilson, 2018-2019, University of Washington Astronomy

ALGORITHM:
    1) Population 100 mocks
    2) Randomly sample N mocks (N=20)
    3) Do step 2 M times (M=20). Each time calculate power spectrum
    4) Calculate variance of all power spectra

Usage:
    ./bootstrap_lyb.py obs uncorr 0 1

Purpose:
    blah
"""
import numpy as np
import pandas as pd
import sys
from scipy.interpolate import griddata
import options as opt

tag=sys.argv[1]
corr_tag = str(sys.argv[2])
save_data = bool(float(sys.argv[3]))
M = int(sys.argv[4]) #M Bootsrap Samples

lya_rest_wave = opt.lya_rest # angstroms
lyb_rest_wave = opt.lyb_rest
lya_min_wave = opt.lya_min # 1045
lya_max_wave = opt.lya_max
lyb_min_wave = opt.lyb_min
lyb_max_wave = opt.lyb_max

if tag == "obs":
    path_to_cat = "../Data/XQ-100_catalogue.txt"
    N = 100
    data_path_name = "../Data/XQ-100/released/*"
if tag == "mocks":
    which_catalog = int(corr_tag.endswith("n5000"))
    path_to_cat = ["XQ-100_catalogue_n100.mock","XQ-100_catalogue_n5000.mock"][which_catalog]
    N = int(path_to_cat.split("_n")[1][:-5])
    data_path_name = "XQ-100_{0}/lyb/ref/spectra/mock-*".format(corr_tag)
saving_here = "results/pk_{0}_{1}.csv".format(tag,corr_tag)

print("Data origin:         ", tag)
print("Data path            ", data_path_name)
print("Draws                ", N)
print("Bootstrap Samples    ", M)
print("Corrections          ", corr_tag)
# print("Saving figure here:  ", saving_here)
print("Saving data:         ", save_data)
print("Saving results here: ", saving_here)

def load_catalog(path=path_to_cat):
    """
    Input catalog path
    Output qso names and qso redshifts in catalog
    """
    catalog = pd.read_csv(path,delim_whitespace=True,usecols=(3,6), names = ("qso_name","redshift"))
    return catalog["qso_name"],catalog["redshift"]

def find_configuration_data(name,zqso):
        """
        Input quasar name and quasar redshift
        Output relative flux fluctuation, pixel redshifts, and normalized flux in LYA FOREST ONLY
        """
        if tag == "obs":
            FILE_r = "../Data/XQ-100/released/%s_uvb-vis.txt"%(name)
            FILE_c = "../Data/XQ-100/continuum/%s_cont.txt"%(name)
            released_filedata = pd.read_csv(FILE_r, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
                                            names = ("wavelength", "nflux","ferr","resolution","dloglambda"))
            cont = pd.read_csv(FILE_c, delim_whitespace=True, skiprows = 1,usecols=[0,1], names = ("wav",'flx'))
            wavelength = released_filedata.wavelength.values
            nflux = released_filedata.nflux.values/cont.flx.values
            ferr = released_filedata.ferr.values/cont.flx.values


        elif tag == "mocks":
            FILE = "XQ-100_{0}/lyb/ref/spectra/{1}_xq{2}_{0}.txt".format(corr_tag,name,N)
            released_filedata = pd.read_csv(FILE, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
                                            names = ("wavelength", "nflux","ferr","resolution","dloglambda"))
            wavelength = released_filedata.wavelength.values
            nflux = released_filedata.nflux.values
            ferr = released_filedata.ferr.values

        #wavelength = released_filedata.wavelength.values
        resolution = released_filedata.resolution.values#released_filedata["resolution"].values
        #nflux = released_filedata.nflux.values/cont.flx.values
        obs_LYB_MIN = lyb_min_wave*(1+zqso)
        obs_LYB_MAX = lyb_max_wave*(1+zqso)
        obs_LYA_MIN = lya_min_wave*(1+zqso)
        obs_LYA_MAX = lya_max_wave*(1+zqso)
        mask_lyb = (wavelength>obs_LYB_MIN)&(wavelength<obs_LYB_MAX)&(nflux>opt.min_trans)
        mask_lya = (wavelength>obs_LYA_MIN)&(wavelength<obs_LYA_MAX)&(nflux>opt.min_trans)
        wavelength_lya = wavelength[mask_lya]
        wavelength_lyb = wavelength[mask_lyb]
        nflux_lyb = nflux[mask_lyb]
        nflux_lya = nflux[mask_lya]

        # NEW
        ferr_lya =  ferr[mask_lya]
        ferr_lyb = ferr[mask_lyb]

        res_lya = resolution[mask_lya]
        res_lyb = resolution[mask_lyb]
        dloglambda_lyb = released_filedata.dloglambda.values[mask_lyb]
        dloglambda_lya = released_filedata.dloglambda.values[mask_lya]
        z_lya = wavelength_lya/lya_rest_wave - 1
        z_lyb = wavelength_lyb/lyb_rest_wave - 1

        return z_lya,z_lyb,nflux_lya,nflux_lyb,res_lya,res_lyb,dloglambda_lya,dloglambda_lyb

def create_configuration_object_boot(nz_catalog,N):
    """
    No input
    Output pandas table containing qso name, redshift, fluctuation array, redshift array and normalized flux array
    """
    #namez_cat = load_catalog() # (name, zqso)
    # print(np.array(nz_catalog[0][0]))
    # single_name_arr = np.array(nz_catalog[0][x])
    # single_z_arr = np.array(nz_catalog[1][x])
    # if len(single_name_arr)>1:
    #     single_name_arr = single_name_arr[0]
    #     single_z_arr = single_z_arr[0]
    # np.unique(nz_catalog)
    #print()
    #N = len(np.unique(nz_catalog[0]))
    #print(np.unique(nz_catalog[0]))
    fz_pix = [find_configuration_data(nz_catalog[0].values[x],nz_catalog[1].values[x]) for x in range(N)]

    return pd.DataFrame({'QSO_NAME':nz_catalog[0],
                  'QSO_z':nz_catalog[1],
                  'z_lya':[fz_pix[i][0] for i in range(N)],
                  'z_lyb':[fz_pix[i][1] for i in range(N)],
                  'nflux_lya':[fz_pix[i][2] for i in range(N)],
                  'nflux_lyb':[fz_pix[i][3] for i in range(N)],
                  'res_lya':[fz_pix[i][4] for i in range(N)],
                  'res_lyb':[fz_pix[i][5] for i in range(N)],
                  'dloglambda_lya':[fz_pix[i][6] for i in range(N)],
                  'dloglambda_lyb':[fz_pix[i][7] for i in range(N)]},
                  columns=['QSO_NAME','QSO_z','z_lya','z_lyb','nflux_lya','nflux_lyb','res_lya',
                           'res_lyb','dloglambda_lya','dloglambda_lyb'])

def power_spectrum_fft(c_fluc,dloglambda):
    """
    Input relative flux fluctuations in configuration space
    Output dictionary containing wavenumber, k, and corresponding power spectrum for each k value
    """
    dv = opt.c_kms * dloglambda * np.log(10)
    k_fluc = np.fft.fft(c_fluc)*dv
    N = len(k_fluc)
    V = N*dv
    dk = 2*np.pi/(N*dv)
    k = dk*np.arange(0,N,1)
    Pk = np.abs(k_fluc)**2/V
    dictionary = {'k':k,'Pk':Pk}
    return dictionary

def grid_interp(z1,z2,rf1,rf2):
    #which_one=1 # 1 means changing z1 & rf1, so z1 and rf2 are the references
    #              0 means changing z2 & rf2 so z1 and rf1 are the references
    if (min(z1)<=min(z2))&(max(z1)>=max(z2)):
        changing_this = 0
        ref_grid, other_pts, other_data = z1,z2,rf2
    elif (min(z2)<min(z1))&(max(z2)>max(z1)):
        changing_this = 1
        ref_grid, other_pts, other_data = z2,z1,rf1
    elif len(z1)<=len(z2):
        changing_this = 0
        ref_grid, other_pts, other_data = z1,z2,rf2
    elif len(z1)>len(z2):
        changing_this = 1
        ref_grid, other_pts, other_data = z2,z1,rf1
    else:
        changing_this = 0
        ref_grid, other_pts, other_data = z1,z2,rf2

    new_rf = griddata(other_pts,other_data,ref_grid,method='linear')
    return ref_grid,new_rf,changing_this

def cross_pk_fft(z_lya,z_lyb,rf_lya,rf_lyb,dloglam_lya,dloglam_lyb,res_lya,res_lyb):
    #c_fluc_lya,c_fluc_lyb,dloglambda_lya,dloglambda_lyb):
    """
    Input relative flux fluctuations in configuration space
    Output dictionary containing wavenumber, k, and corresponding power spectrum for each k value
    """
    x,y,which = grid_interp(z_lya,z_lyb,rf_lya,rf_lyb) #reference, changed, which one
    finite_mask = np.isfinite(y)
    #print(len(x),len(y),len(finite_mask))
    if which == 1: #lyb is referece -- confirmed
        z_lyb,rf_lyb = z_lyb[finite_mask],rf_lyb[finite_mask]
        z_lya,rf_lya,dloglam,res = x[finite_mask],y[finite_mask],dloglam_lyb[finite_mask],res_lyb[finite_mask]
    elif which == 0: #lya is reference -- confirmed
        z_lya,rf_lya = z_lya[finite_mask],rf_lya[finite_mask]
        z_lyb,rf_lyb,dloglam,res = x[finite_mask],y[finite_mask],dloglam_lya[finite_mask],res_lya[finite_mask]
        #        print(len(rf_lya),len(rf_lyb),len(dloglam))
    dv = opt.c_kms * dloglam * np.log(10)

    kf_lya = np.fft.fft(rf_lya)*dv #dloglam
    kf_lyb = np.fft.fft(rf_lyb)*dv #dloglam
    N=len(kf_lya)
    V = N*dv
    dk = 2*np.pi/(N*dv)
    k=dk*np.arange(0,N,1)
    P_cross = np.conj(kf_lya)*kf_lyb/V
    P_real,P_imag = P_cross.real,P_cross.imag
    #     print(len(k),len(finite_mask),len(z_lyb))
    dictionary = {'k':k,'P_real':P_real,'P_imag':P_imag,'dloglam':dloglam,'res':res}
                  #'grid_mask':finite_mask}
    return dictionary

def window(k,p,R):
    """
    Deconvolution kernel.
    k: velocity wavenumber
    p: pixel width (dloglambda)
    R: resolution
    """
    #p = p*opt.c_kms*np.log(10) #Check this part with somebody
    gaussian = np.exp(-0.5*k**2*R**2)
    tophat =  np.sinc(0.5*k*p) #np.sin(k*p/2)/(k*p/2)
    deconvolution_kernel = gaussian*tophat
    return deconvolution_kernel

# mf_data_table = "results/mf_{0}.csv".format(tag)
mean_flux_table = pd.read_csv("results/mf_{0}.csv".format(tag))
nz_catalog = load_catalog(path=path_to_cat)
kbin_centers = opt.kbin_centers
zbin_centers = opt.zbin_centers
zbin_edges = opt.zbin_edges
kbin_edges = opt.kbin_edges
zbinlen = opt.zbinlen
kbinlen = opt.kbinlen

np.random.seed(1)
Paa_table = []
Pbb_table = []
Pab_table = []
Qab_table = []

#test = 0
for i in range(M):
    opt.updt(M, i)
    mask = np.floor(np.random.rand(N)*N).astype(int)
    #print(mask)
    #np.random.choice(100,N,False)#np.random.randint(0,100,N)
    iter_catalog = nz_catalog[0][mask],nz_catalog[1][mask]
    #print(iter_catalog)
    pixel_catalog = create_configuration_object_boot(iter_catalog,N).reset_index()
    #print(pixel_catalog)

    pk_table_lya = []
    pk_table_lyb = []
    pk_real_cross_table = []
    pk_imag_cross_table = []

    # NEW
    rlya_fluc_matrix = [[] for x in range(zbinlen)] # squared
    rlyb_fluc_matrix = [[] for x in range(zbinlen)] #

    for zindex in range(zbinlen):
        # reading in measurements
        kpower_matrix_lya = [[] for x in range(kbinlen)]
        kpower_matrix_lyb = [[] for x in range(kbinlen)]
        real_power_matrix  = [[] for x in range(kbinlen)]
        imag_power_matrix  = [[] for x in range(kbinlen)]

        for qindex in range(N):
            #test=test+1
            #opt.updt(zbinlen*N,test)
            z_lya = pixel_catalog.z_lya[qindex]
            z_lyb = pixel_catalog.z_lyb[qindex]
            #zb_min,zb_max = np.min(z_lyb),np.max(z_lyb)
            nflux_lya = pixel_catalog.nflux_lya[qindex]
            nflux_lyb = pixel_catalog.nflux_lyb[qindex]
            res_lya = pixel_catalog.res_lya[qindex]
            res_lyb = pixel_catalog.res_lyb[qindex]
            dloglambda_lya = pixel_catalog.dloglambda_lya[qindex]
            dloglambda_lyb = pixel_catalog.dloglambda_lyb[qindex]

            zmask_lya = (z_lya>zbin_edges[zindex])&(z_lya<=zbin_edges[zindex+1])
            zmask_lyb = (z_lyb>zbin_edges[zindex])&(z_lyb<=zbin_edges[zindex+1])
            #zmask_lya_subset = (z_lya>min(z_lyb))#&(z_lya<=zb_max) #[x1>min(x2)]
            #zmask_lyb_subset = (z_lyb<max(z_lya))#&(z_lya<=zb_max) #[x1>min(x2)]
            #print(len(zlya[zmask_lya_subset]),len(z_lyb[zmask_lyb]), qindex,zindex)

            rfluc_lya = nflux_lya/mean_flux_table.meanF_lya[zindex] - 1
            rfluc_lyb = nflux_lyb/mean_flux_table.meanF_lyb[zindex] - 1

            if len(rfluc_lya[zmask_lya])>opt.min_pix:
                rlya_fluc_matrix[zindex].append(rfluc_lya[zmask_lya])

                pk_ktable = power_spectrum_fft(rfluc_lya[zmask_lya],dloglambda_lya[zmask_lya])
                k_lya,pk_lya = pk_ktable['k'],pk_ktable['Pk']
                for kindex in range(kbinlen):
                    kmask = (k_lya>kbin_edges[kindex])&(k_lya<=kbin_edges[kindex+1])&(k_lya!=0)
                    ksubset = k_lya[kmask]
                    kpower_subset = pk_lya[kmask]
                    p_wave = dloglambda_lya[zmask_lya][kmask]
                    res_wave = res_lya[zmask_lya][kmask]
                    if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrN'):
                        #print("{0:.3e} {1:.3e}".format(kpower_subset[2],mean_flux_table.powerN_lya[zindex]))
                        kpower_subset = kpower_subset - mean_flux_table.powerN_lya[zindex]

                        #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
                    if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrR'):
                        kpower_subset = kpower_subset / window(k=ksubset,p=p_wave,R=res_wave)**2
                        #pk_lyb = pk_lyb / window(k=k_lyb,p=dloglam_lyb,R=res_lyb)**2
                    kpower_matrix_lya[kindex].append(kpower_subset)
                    #print(window(k=ksubset,p=p_wave,R=res_wave))

            if len(rfluc_lyb[zmask_lyb])>opt.min_pix:
                rlyb_fluc_matrix[zindex].append(rfluc_lyb[zmask_lyb])

                pk_ktable = power_spectrum_fft(rfluc_lyb[zmask_lyb],dloglambda_lyb[zmask_lyb])
                k_lyb,pk_lyb = pk_ktable['k'],pk_ktable['Pk']
                for kindex in range(kbinlen):
                    kmask = (k_lyb>kbin_edges[kindex])&(k_lyb<=kbin_edges[kindex+1])&(k_lyb!=0)
                    ksubset = k_lyb[kmask]
                    kpower_subset = pk_lyb[kmask]
                    p_wave = dloglambda_lyb[zmask_lyb][kmask]
                    res_wave = res_lyb[zmask_lyb][kmask]
                    if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrN'):
                        kpower_subset = kpower_subset - mean_flux_table.powerN_lya[zindex]
                        #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
                    if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrR'):
                        kpower_subset = kpower_subset / window(k=ksubset,p=p_wave,R=res_wave)**2
                        #print(kpower_subset[:2])
                        #pk_lyb = pk_lyb / window(k=k_lyb,p=dloglam_lyb,R=res_lyb)**2
                    kpower_matrix_lyb[kindex].append(kpower_subset)

            if (len(rfluc_lyb[zmask_lyb])>opt.min_pix)&(len(rfluc_lya[zmask_lya])>opt.min_pix):
                #print(len(rfluc_lya[zmask_lya]),len(rfluc_lyb[zmask_lyb]))
                #print(len(rfluc_lya[zmask_lya]),len(rfluc_lyb[zmask_lyb]))
                #k_lya,k_lyb = k_lya[zmask_lya],k_lyb[zmask_lyb]
                z_lya,z_lyb = z_lya[zmask_lya],z_lyb[zmask_lyb]
                rf_lya,rf_lyb = rfluc_lya[zmask_lya],rfluc_lyb[zmask_lyb]
                dloglam_lya,dloglam_lyb = dloglambda_lya[zmask_lya],dloglambda_lyb[zmask_lyb]
                res_lya,res_lyb = res_lya[zmask_lya],res_lyb[zmask_lyb]
                new_lya_mask = (z_lya>min(z_lyb))&(z_lya<max(z_lyb))
                new_lyb_mask = (z_lyb>min(z_lya))&(z_lyb<max(z_lya))
                #z1,z2,rf1,rf2 = z1[mask1],z2[mask2],rf1[mask1],rf2[mask2]
                pk_ktable = cross_pk_fft(z_lya[new_lya_mask],z_lyb[new_lyb_mask],
                                         rf_lya[new_lya_mask],rf_lyb[new_lyb_mask],
                                         dloglam_lya[new_lya_mask],dloglam_lyb[new_lyb_mask],
                                         res_lya[new_lya_mask],res_lyb[new_lyb_mask])
                #k,P_real,P_imag,grid_mask = pk_ktable['k'],pk_ktable['P_real'],pk_ktable['P_imag'],pk_ktable['grid_mask']
                #k,P_real,P_imag
                k,P_real,P_imag,dloglam,res = pk_ktable['k'],pk_ktable['P_real'],pk_ktable['P_imag'],pk_ktable['dloglam'],pk_ktable['res']

                ### TEST HERE
                # Noise Correction
                if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrN'):
                    #powerk = powerk - mean_flux_table.powerN[zindex]
                    #pk_lya = pk_lya - mean_flux_table.powerN_lya[zindex]
                    #pk_lyb = pk_lyb - mean_flux_table.powerN_lyb[zindex]
                    P_real = P_real - mean_flux_table.powerN_lya[zindex]
                    P_imag = P_imag - mean_flux_table.powerN_lya[zindex]
                    #print("!!!")

                # Beam Correction
                if np.any(corr_tag=='corrNR') or np.any(corr_tag=='corrR'):
                    # print(P_real[10:13])
                    #print(res_lya[25:30])
                    #pk_lya = pk_lya / window(k=k_lya,p=dloglam_lya,R=res_lya)**2
                    #pk_lyb = pk_lyb / window(k=k_lyb,p=dloglam_lyb,R=res_lyb)**2
                    #print(len(k),len(dloglam_lya[new_lya_mask]),len(res_lya[new_lya_mask]),len(P_real))
                    P_real = P_real / window(k=k,p=dloglam,R=res)**2
                    P_imag = P_imag / window(k=k,p=dloglam,R=res)**2
                    # print(P_real[10:13])
                    #print((window(k=k,p=dloglam,R=res)**2)[:3])
                    #print("000")
                    #print(res_lya[25:30])

                ### TEST HERE
                for kindex in range(kbinlen):
                    kmask = (k>kbin_edges[kindex])&(k<=kbin_edges[kindex+1])&(k!=0)
                    ksubset = k[kmask]
                    pk_sub_real = P_real[kmask]
                    pk_sub_imag = P_imag[kmask]
                    real_power_matrix[kindex].append(pk_sub_real)
                    imag_power_matrix[kindex].append(pk_sub_imag)

        pk_means_lya = []
        pk_means_lyb = []
        pk_means_real = []
        pk_means_imag = []
        for columnindex in range(kbinlen):
            try:
                mean_lya = np.mean(np.concatenate(kpower_matrix_lya[columnindex]))
            except:
                mean_lya=np.nan
            try:
                mean_lyb = np.mean(np.concatenate(kpower_matrix_lyb[columnindex]))
            except:
                mean_lyb = np.nan
            try:
                mean_real = np.mean(np.concatenate(real_power_matrix[columnindex]))
            except:
                mean_real = np.nan
            try:
                mean_imag = np.mean(np.concatenate(imag_power_matrix[columnindex]))
            except:
                mean_imag = np.nan

            pk_means_lya.append(mean_lya)
            pk_means_lyb.append(mean_lyb)
            pk_means_real.append(mean_real)
            pk_means_imag.append(mean_imag)
        Paa_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
               kbin_centers,np.array(pk_means_lya)])
        Pbb_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
               kbin_centers,np.array(pk_means_lyb)])
        Pab_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
               kbin_centers,np.array(pk_means_real)])
        Qab_table.append([zbin_centers[zindex]*np.ones_like(kbin_centers),
               kbin_centers,np.array(pk_means_imag)])

Paa_cols = np.column_stack(Paa_table)
Pbb_cols = np.column_stack(Pbb_table)
Pab_cols = np.column_stack(Pab_table)
Qab_cols = np.column_stack(Qab_table)

k_exact = list(set(Paa_cols[1]))
z_exact = list(set(Paa_cols[0]))
z_exact.sort()
k_exact.sort()
mydict = {"var_Paa":[],"var_Pbb":[],"var_Pab":[],"var_Qab":[],"k":[],"z":[]}#{"Paa_var":[],"k":[],"z":[],"pk":[]}
for i in k_exact:
    for j in z_exact:
        masking_Paa = (Paa_cols[0] == j)&(Paa_cols[1]==i)
        Paa_masked = Paa_cols[2][masking_Paa]
        masking_Pbb = (Pbb_cols[0] == j)&(Pbb_cols[1]==i)
        Pbb_masked = Pbb_cols[2][masking_Pbb]
        masking_Pab = (Pab_cols[0] == j)&(Pab_cols[1]==i)
        Pab_masked = Pab_cols[2][masking_Pab]
        masking_Qab = (Qab_cols[0] == j)&(Qab_cols[1]==i)
        Qab_masked = Qab_cols[2][masking_Qab]

        #mydict["pk"].append(np.mean(pk_masked))
        mydict["var_Paa"].append(np.var(Paa_masked))
        mydict["var_Pbb"].append(np.var(Pbb_masked))
        mydict["var_Pab"].append(np.var(Pab_masked))
        mydict["var_Qab"].append(np.var(Qab_masked))
        mydict["k"].append(i)
        mydict["z"].append(j)

save_table = pd.DataFrame(mydict)
save_table.sort_values(by=['z','k'], inplace=True)
save_table.reset_index(drop=True,inplace=True)
# save_table.to_csv(
# "../../results/boot_lya_mocks_{0}_M{1}.csv".format(tag,M))

# pk_vals = Paa_cols[2]
new_table = pd.read_csv(saving_here)
# print(new_table)
# print(save_table)
new_table['err_Paa'] = np.sqrt(save_table.var_Paa)
new_table['err_Pbb'] = np.sqrt(save_table.var_Pbb)
new_table['err_Pab'] = np.sqrt(save_table.var_Pab)
new_table['err_Qab'] = np.sqrt(save_table.var_Qab)
if save_data:
    new_table.to_csv(saving_here,index=False)
#print(new_table)
#save_pk = pd.DataFrame(pk_vals)
#save_pk.to_csv(saving_here)
opt.updt(M, M)
# print("File is located here: ")
#print("../../results/boot_lya_mocks_{0}_M{1}.csv".format(tag,M))

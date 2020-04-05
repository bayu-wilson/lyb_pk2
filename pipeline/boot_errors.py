#!/usr/bin/env python

import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
import pandas as pd

#OPTIONS

add_additional_boostrap_error = 0 #NOT CURRENTLY DONE
add_mean_flux_error = 0  #currently no sample variance error and we find this error is small
add_res_error = 1
add_metal_error = 1  #also requires inis.subtract_metal_power 

add_additional_boostrap_error_factor = 1.2
extra_error_factor_meanflux = 1.2
dlogsigma = 0.05 #fractional error on resolution sigma
metal_error_Factor = 1.0  #factor by which amplitude of metal power is uncertain (motivated by BOSS measurements)

###################################################################
#Read in power spectra and mean flux measurements
######################################################################
mfest = pd.read_csv(inis.save_mf_path)
pkest = pd.read_csv(inis.save_pk_path)

########################################################################
#Read in bootstraps
#########################################################################
mf = np.loadtxt(inis.save_boot_mf_path)
p = np.loadtxt(inis.save_boot_pk_path)

#############################

if inis.subtract_metal_power and add_metal_error:
    metal_power = np.concatenate([np.loadtxt('../data/obs/pk_xs_avg.txt')]*opt.zbinlen) #metal pwoer spectrum
    meta_power_error = metal_error_Factor
    
    #x.Paa = x.Paa.values-metal_power
    #x.Ptot = x.Ptot.values-metal_power


###################################################################
#Make covariances
####################################################################
mfboot = np.average(mf, axis=1)
error_mf = np.sqrt(np.cov(mf).diagonal())
covpk = np.cov(p)
pkboot = np.average(p, axis=1)

covpk0 = np.cov(p) #uncorrected

zmfest_array =  np.concatenate([mfest.z.values.T,  mfest.z.values.T])
mfest_array =  np.concatenate([mfest.mf_a.values.T,  mfest.mf_tot.values.T])


pkest_array =  np.concatenate([pkest.paa.values.T,  pkest.ptt.values.T, pkest.pab.values.T])
kest_arr =  np.concatenate([pkest.k.values.T,  pkest.k.values.T, pkest.k.values.T])
zpkest_arr = np.concatenate([pkest.z.values.T,  pkest.z.values.T, pkest.z.values.T])
#print("pkest =", (pkest_array- pkboot)/pkest_array, len(zpkest_arr)) #check to make sure understand things



############################################################################################
#add in mean flux error
##########################################################################################

def getIndex(i, zbin, num):
    if i < num:
        iX = zbin
        iY = zbin
        lambdaX = opt.lya_rest
        lambdaY = opt.lya_rest
    elif i < 2*num:
        iX = zbin+ opt.zbinlen
        iY = zbin+ opt.zbinlen
        lambdaX = opt.lyb_rest
        lambdaY = opt.lyb_rest        
    else:
        iX = zbin
        iY = zbin + opt.zbinlen
        lambdaX = opt.lya_rest
        lambdaY = opt.lyb_rest   

    armX =  'VIS' if opt.overlap_maxwav < lambdaX*(1+opt.zbin_centers[zbin]) else 'UV'
    armY =  'VIS' if opt.overlap_maxwav < lambdaY*(1+opt.zbin_centers[zbin]) else 'UV'
    sigmaX = opt.R_VIS_carswell  if opt.overlap_maxwav < lambdaX*(1+opt.zbin_centers[zbin]) else opt.R_UV_carswell  #Bob's resolution seems like good approximation
    sigmaY = opt.R_VIS_carswell  if opt.overlap_maxwav < lambdaY*(1+opt.zbin_centers[zbin]) else opt.R_UV_carswell
    return iX, iY, sigmaX, sigmaY, armX, armY


N = len(pkest_array)
    
print(N, opt.kbinlen*opt.zbinlen)
for i in range(N):
    for j in range(N):
        #print(i, j, N)
        iX,iY, sigiX, sigiY, armiX, armiY = getIndex(i, (i//opt.kbinlen)%opt.zbinlen, N//3)
        jX,jY, sigjX, sigjY, armjX, armjY = getIndex(j, (j//opt.kbinlen)%opt.zbinlen, N//3)
        ki = opt.kbin_centers[i%opt.kbinlen];kj = opt.kbin_centers[j%opt.kbinlen];
            
        if add_mean_flux_error:
            #print(i, j, iX, iY, jX, jY, (((iX==jX)+(iX==jY))),  (error_mf[iX]/mfest_array[iX]) )
            x = (((iX==jX)+(iX==jY))*(error_mf[iX]/mfest_array[iX])**2 + ((iY==jX)+(iY==jY))*(error_mf[iY]/mfest_array[iY])**2)*pkest_array[i]*pkest_array[j]
            x *= extra_error_factor_meanflux*extra_error_factor_meanflux
            covpk[i, j] += x#only diagonal terms in mean flux correlate
            #if i == j:
            #    print(opt.zbin_centers[(i//opt.kbinlen)%opt.zbinlen], i, j, np.sqrt(x)/pkest_array[i], np.sqrt(covpk[i, j]/pkest_array[i]), error_mf[iX], mfest_array[iX], error_mf[iY], mfest_array[iY])
        if add_res_error: #only correlates between arms
            x = dlogsigma**2*ki**2*kj**2*pkest_array[i]*pkest_array[j]*(sigiX**4*((armiX == armjX)+(armiX == armjY))+sigiY**4*((armiY == armjX)+(armiY == armjY))) 
            covpk[i, j] += x

            #if i == j:
            #    print(opt.zbin_centers[(i//opt.kbinlen)%opt.zbinlen], i, j, sigiX, sigiY, (sigiX**4*((armiX == armjX)+(armiX == armjY))+sigiY**4*((armiY == armjX)+(armiY == armjY))))
########################################################################
#add uncertainty to mean flux in bootstap errors
########################################################################
#if add_meanfluxerror_bootstrap:
####################################################################
#add 10% uncertainty for resolution
####################################################################

        

N_kz = opt.zbinlen*opt.kbinlen
pk_err_diag = covpk.diagonal()**0.5
pk0_err_diag = covpk0.diagonal()**0.5
err_paa = pk_err_diag[0:N_kz]; err0_paa = pk0_err_diag[0:N_kz]
print("paa = ", (pk_err_diag-pk0_err_diag)/pkest_array) #, (pk_err_diag-pk0_err_diag)/pk0_err_diag)
err_ptt = pk_err_diag[N_kz:N_kz*2]
err_pab = pk_err_diag[N_kz*2:N_kz*3]

err_paa_sub = err_paa[0*opt.kbinlen:3*opt.kbinlen]
err_ptt_sub = err_ptt[4*opt.kbinlen:7*opt.kbinlen]
err_pbb_sub = np.sqrt(err_paa_sub**2+err_ptt_sub**2)
err_pbb = np.ones(N_kz)*np.nan
err_pbb[4*opt.kbinlen:7*opt.kbinlen] = err_pbb_sub

if inis.save_pk_with_err:
    columns = ['k','z','paa','err_paa','N_aa','ptt', 'err_ptt','N_tt','pab','err_pab','N_ab','pbb','err_pbb']
    pk_everything = np.column_stack((pkest.k,pkest.z,
                    pkest.paa, err_paa, pkest.npix_aa,
                    pkest.ptt, err_ptt, pkest.npix_tt,
                    pkest.pab,err_pab, pkest.npix_ab,
                    pkest.pbb, err_pbb,))
    df_pk = pd.DataFrame(pk_everything,columns=columns)

    optionsstr = '';
    if add_mean_flux_error:
        optionsstr = '_mf' + str(extra_error_factor_meanflux);
    if add_res_error:
        optionsstr = optionsstr + '_dlogsigerr' + str(dlogsigma)
    if inis.subtract_metal_power and add_metal_error:
        optionsstr = optionsstr + '_wmetalerror'
        
    outputfile = inis.save_pk_with_err_path  + optionsstr + '.csv'
    df_pk.to_csv(outputfile, index=False)


    print("Saving new datatables here:\n{0}".format(outputfile))

# ### INCLUDING RESOLUTION UNCERTAINTY
# sigma_res_aa = opt.find_res_uncertainty(pkest.k,pkest.z,pkest.Paa)
# err_paa = np.sqrt(err_paa**2+sigma_res_aa**2)
# sigma_res_tt = opt.find_res_uncertainty(pkest.k,pkest.z,pkest.Ptot)
# err_ptt = np.sqrt(err_ptt**2+sigma_res_tt**2)
# sigma_res_ab = opt.find_res_uncertainty(pkest.k,pkest.z,pkest.Pab)
# err_pab = np.sqrt(err_pab**2+sigma_res_ab**2)
# sigma_res_bb = opt.find_res_uncertainty(pkest.k,pkest.z,pkest.Pbb)
# err_pbb = np.sqrt(err_pbb**2+sigma_res_bb**2)

# N = 91
# boot_mat = np.reshape(boot_arr,(M,N))
# # [np.mean(boot_mat[:,i]) for i in range(91)]
# # np.reshape(boot_arr,(10,91))[:,90]
# # np.array([paa[:,i] - np.mean(paa,axis=1) for i in range(3000)]).T
# # np.array([boot_mat.T[:,i] - np.mean(boot_mat.T,axis=1) for i in range(10)]).T
#
# paa_minus_paamean = np.array([boot_mat[i] - np.mean(boot_mat,axis=0) for i in range(M)]).T
# Cij = np.zeros((N,N))
# for i in range(N):
#     for j in range(N):
#         Cij[i][j] = np.mean(paa_minus_paamean[i]*paa_minus_paamean[j])
# rij = np.zeros((N,N))
# for i in range(N):
#     for j in range(N):
#         rij[i][j] = Cij[i][j]/np.sqrt(np.abs(Cij[i][i]*Cij[j][j]))
# plt.imshow(rij)
# plt.show()

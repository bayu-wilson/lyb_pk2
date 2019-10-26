# #!/usr/bin/env python
#
#
# ################################################
# ####### BOOT REALIZATIONS. NOT BOOT INDO  ######
# ################################################
#
#
import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
import pandas as pd


realization_table_pk = np.zeros((50,21,13)) # 50 realizations
realization_table_var = np.zeros((50,21,13)) # 50 realizations


table_mf =  pd.read_csv("../output/qsos_mf_lyb_nocorr_n5000.txt",delim_whitespace=True,
                             names=['qso_index','z','sum_F_lya','pix_F_lya','sum_F_lyb','pix_F_lyb'])
#quasar index, z-index, sum of lya flux in z-bin, sum/num of lya pixels in z-bin,
#sum of lyb (total) flux in z-bin, sum/num of lyb (total) pixels in z-bin
table_pk = pd.read_csv("../output/qsos_pk_lyb_nocorr_n5000.txt",delim_whitespace=True,
                             names=['qso_index','msrmnt','z','k','pk_bin','pix_bin'])
#quasar index, type of power measurement (numbered from 0-2 for paa,ptt, and pab),
#redshift, wavenumber (k), sum of power in k,z bin, sum of pixels in k,z bin.

# opt.zbinlen = 5 #sept12
M = 1000 #inis.M #number of bootstrap samples
nqso = 100 #inis.nqso

# #realizations = np.zeros((50,100),dtype=int)
for r in range(50):
    fpath_realization="../data/realizations/QSO_catalogue/xq100_realizations/XQ-100_catalogue_r{:02}.mock".format(r)
    rtable = pd.read_csv(fpath_realization,delim_whitespace=True,usecols=(3,6),names = ['name','redshift'])
    realization_idcs = np.array([rtable.name.values[i][5:] for i in range(len(rtable))],dtype=int)
    mask_mf = np.concatenate([np.where(table_mf.qso_index.values == realization_idcs[q])[0] for q in range(100)])
    mask_pk = np.concatenate([np.where(table_pk.qso_index.values == realization_idcs[q])[0] for q in range(100)])

    subset_table_mf = table_mf.iloc[mask_mf]
    subset_table_pk = table_pk.iloc[mask_pk]

    np.random.seed(1)
    rand_matrix = [np.floor(np.random.rand(nqso)*nqso).astype(int) for i in range(M)]
    boot_mf_arr = np.zeros((M,(opt.zbinlen)*2)) #empty matrix to be filled with mf values for each bootstrap sample
    boot_pk_arr = np.zeros((M,opt.zbinlen*3,opt.kbinlen))

    for m in range(M):
        opt.updt(50*M,r*M + m)
        rand_realization = realization_idcs[rand_matrix[m]]
        #rand_sorted_idcs_mf = subset_table_mf.qso_index.values[::7][rand_matrix[m]]

        mask_rand_mf = np.concatenate([np.where(subset_table_mf.qso_index == rsi)[0] for rsi in rand_realization])
        temp_mf_table = subset_table_mf.iloc()[mask_rand_mf].groupby(['z']).sum()
        temp_mf_table.sum_F_lya = temp_mf_table.sum_F_lya/temp_mf_table.pix_F_lya
        temp_mf_table.sum_F_lyb = temp_mf_table.sum_F_lyb/temp_mf_table.pix_F_lyb
        boot_mf_arr[m] = np.concatenate((temp_mf_table.sum_F_lya.values,temp_mf_table.sum_F_lyb))


        mask_rand_pk = np.concatenate([np.where(subset_table_pk.qso_index == rsi)[0] for rsi in rand_realization])
        temp_pk_table = subset_table_pk.iloc()[mask_rand_pk].groupby(['msrmnt','z','k']).sum()
        temp_pk_table.pk_bin = temp_pk_table.pk_bin/temp_pk_table.pix_bin
        #boot_pk_arr[m] = temp_pk_table.pk_bin.values
#         temp_mf_table.sum_F_lyb = temp_mf_table.sum_F_lyb/temp_mf_table.pix_F_lyb
#         boot_mf_arr[m] = np.concatenate((temp_mf_table.sum_F_lya.values,temp_mf_table.sum_F_lyb))

        pk_arr = np.zeros((opt.zbinlen*3,opt.kbinlen))
        #pk_arr[pidx*opt.zbinlen+zidx]
        for pidx in range(3):
            for zidx in range(opt.zbinlen):
                for kidx in range(opt.kbinlen):
                    new = temp_pk_table.reset_index()
                    mask = [(new.msrmnt.values==pidx)&
                            (new.z.values==opt.zbin_centers[zidx])&
                            (np.round(new.k.values,5)==np.round(opt.kbin_centers[kidx],5))]
                    try:
                        pk_arr[pidx*opt.zbinlen+zidx][kidx] = new.pk_bin.values[mask]
                    except:
                        pk_arr[pidx*opt.zbinlen+zidx][kidx] = np.nan
        boot_pk_arr[m] = pk_arr


        realization_table_pk[r] = np.mean(boot_pk_arr,axis=0)
        realization_table_var[r] = np.var(boot_pk_arr,axis=0)

opt.updt(50*M,50*M)

np.savetxt('../output/realizations/mean_pk.txt',np.column_stack(np.reshape(realization_table_pk,(50,21*13))))
np.savetxt('../output/realizations/var_pk.txt',np.column_stack(np.reshape(realization_table_var,(50,21*13))))



#         temp_table = (subset_table_mf.reset_index().iloc[rand_matrix[m]]).groupby(['z']).sum()
#         #temp_table = temp_table.groupby(['z']).sum()
#         temp_table.sum_F_lya = temp_table.sum_F_lya/boot_table_mf.pix_F_lya
#         temp_table.sum_F_lyb = temp_table.sum_F_lyb/boot_table_mf.pix_F_lyb
#         boot_mf_arr[m] = np.concatenate((temp_table.sum_F_lya.values,
#                                         temp_table.sum_F_lyb.values))
        #subset_table_mf['new_idx'] = np.arange(len(subset_table_mf))
        #subset_table_mf[subset_table_mf.new_idx[rand_matrix[m]]]
#         x = subset_table_mf.groupby(['z']).sum()
#         x.sum_F_lya = x.sum_F_lya/x.pix_F_lya
#         x.sum_F_lyb = x.sum_F_lyb/x.pix_F_lyb


        #subset_table_pk.grou
        #table_pk.groupby(['msrmnt','z','k']).sum()
# del temp_pk_table['qso_index']

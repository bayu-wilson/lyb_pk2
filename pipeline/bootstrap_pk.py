#!/usr/bin/env python

import numpy as np
import inis
import pandas as pd
import options as opt
from QuasarSpectrum import QuasarSpectrum
from scipy.interpolate import griddata


tag = inis.tag
cat_name = inis.cat_name
rescale_flux = inis.rescale_flux

QuasarSpectrum.load_cat(cat_name)
nqso = QuasarSpectrum.nqso

print("Loading Data")
qso_list = []
for i in range(nqso):
    q = QuasarSpectrum.load_qso_data(i,tag=tag,rescale_flux=rescale_flux)
    if "noB" in tag:
        q.get_new_forest(rescale_flux=rescale_flux,wrange =
            (opt.lyb_min,opt.lyb_max,opt.lyb_rest, opt.xs_beta))
        q.get_new_forest(rescale_flux=rescale_flux,wrange =
            (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d1, opt.xs_ovi))
        q.get_new_forest(rescale_flux=rescale_flux,wrange =
            (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d2, opt.xs_ovi))
    qso_list.append(q)
print("Done!\n")

# BOOTSTRAP SPECIFIC
np.random.seed(1)
nz_arr = QuasarSpectrum.all_redshifts,QuasarSpectrum.all_names
M = inis.M # M Samples



mf_msrmnts = ['mf_a', 'nvar_a','dloglambda_a', 'npow_a',
           'mf_tot', 'nvar_tot','dloglambda_tot', 'npow_tot','z','mf_b',
           'var_atot','npow_atot']
n_mf_msrmnts = len(mf_msrmnts)

pk_msrmnts = ['k','Paa','Ptot','Pab','Qab','Pbb','num','z']
n_pk_msrmnts = len(pk_msrmnts)

mf_bootstrap = np.zeros((M,opt.zbinlen,n_mf_msrmnts))
pk_bootstrap = np.zeros((M,opt.kbinlen*opt.zbinlen, n_pk_msrmnts))

# Paa_table = []
# Pbb_table = []
# Pab_table = []
# Qab_table = []


for m in range(M):
    opt.updt(M, m)
    qmask = np.floor(np.random.rand(nqso)*nqso).astype(int)
    qso_arr = np.array(qso_list)[qmask]

    zbin_msr_matrix = np.zeros((opt.zbinlen,n_mf_msrmnts))
    # flux_pdf = [[] for tmp in range(opt.zbinlen)]
    # pdf_bins = np.arange(-0.025,1.05,.05)
    for zidx in range(opt.zbinlen):
        msrmnt_in_zbin =  np.zeros(n_mf_msrmnts)
        count_in_zbin = np.zeros(n_mf_msrmnts)
        #opt.updt(opt.zbinlen, zidx)

        zbin_msrmnt = [[] for idx in range(n_mf_msrmnts)]
        for i in range(nqso):
            #i=85
            #zidx = 6
            zpix_a = qso_arr[i].get_zpix(opt.lya_rest)

            name = qso_arr[i].name
            mask = qso_arr[i].get_zmask(forest=(opt.lya_min,opt.lya_max,opt.lya_rest),
                                         zpix=zpix_a,zidx=zidx,zedges=opt.zbin_edges,name=name)
            #FLUX PDF
            #flux_pdf[zidx].append(np.histogram(qso_arr[i].flux[mask],bins=pdf_bins)[0])#/np.nansum(mask))

            zpix_b = qso_arr[i].get_zpix(opt.lyb_rest) # Here is where I want to change optical depth of lyb pixels
            mask_b = qso_arr[i].get_zmask(forest=(opt.lyb_min,opt.lyb_max,opt.lyb_rest),
                                         zpix=zpix_b,zidx=zidx,zedges=opt.zbin_edges,name=name)

            za = zpix_a[mask]#qso_arr[i].wavelength[mask]/opt.lya_rest-1
            ztot = zpix_b[mask_b]#qso_arr[i].wavelength[mask_b]/opt.lyb_rest-1
            try:
                new_af_mask = (za>np.min(ztot))&(za<np.max(ztot))
                new_bf_mask = (ztot>np.min(za))&(ztot<np.max(za))
                ferra = qso_arr[i].err_flux[mask][new_af_mask]
                ferrtot = qso_arr[i].err_flux[mask_b][new_bf_mask]

                # Interpolating to the smaller one
                if len(ferrtot)<=len(ferra):
                    ferra = griddata(za[new_af_mask],ferra,ztot[new_bf_mask],method='linear')
                    ferra = ferra[np.isfinite(ferra)]
                else:
                    ferrtot = griddata(ztot[new_bf_mask],ferrtot,za[new_af_mask],method='linear')
                    ferrtot = ferrtot[np.isfinite(ferrtot)]
                #print(np.nansum(ferra*ferrtot))
                msrmnt_in_zbin[10]+= np.sum(ferra*ferrtot)*0 #CHANGED 5/3/19 after meeting with Matt
                # var lya-tot 10
                count_in_zbin[10] += len(ferra)                         # len lya_tot 10
            except:
                pass

            msrmnt_in_zbin[0]+= np.sum(qso_arr[i].flux[mask])          # mf lya 0
            count_in_zbin[0] += np.sum(mask)                           # len lya 0
            msrmnt_in_zbin[1]+= np.sum(qso_arr[i].err_flux[mask]**2)   # var lya 1
            count_in_zbin[1] += np.sum(mask)                           # len lya 1
            msrmnt_in_zbin[2]+= np.sum(qso_arr[i].dloglambda[mask])    # dloglam lya 2
            count_in_zbin[2] += np.sum(mask)                           # len lya 2
            msrmnt_in_zbin[4]+= np.sum(qso_arr[i].flux[mask_b])        # mf tot 4
            count_in_zbin[4] += np.sum(mask_b)                         # len tot 4
            msrmnt_in_zbin[5]+= np.sum(qso_arr[i].err_flux[mask_b]**2) # var tot 5
            count_in_zbin[5] += np.sum(mask_b)                         # len tot 5
            msrmnt_in_zbin[6]+= np.sum(qso_arr[i].dloglambda[mask_b])  # dloglam tot 6
            count_in_zbin[6] += np.sum(mask_b)                         # len tot 6
        zbin_msr_matrix[zidx] = msrmnt_in_zbin/count_in_zbin
        #print(count_in_zbin[0])
        #print(msrmnt_in_zbin[0])
    #opt.updt(opt.zbinlen, opt.zbinlen)

    #print("Done!\n")
    # npow alpha 3
    zbin_msr_matrix.T[3] = list(QuasarSpectrum.get_npow(mf=zbin_msr_matrix.T[0],
                                            nvar=zbin_msr_matrix.T[1],
                                            dloglambda=zbin_msr_matrix.T[2]))
    # npow total 7
    zbin_msr_matrix.T[7] = list(QuasarSpectrum.get_npow(mf=zbin_msr_matrix.T[4],
                                            nvar=zbin_msr_matrix.T[5],
                                            dloglambda=zbin_msr_matrix.T[6]))

    # zbins 8
    zbin_msr_matrix.T[8] = opt.zbin_centers

    # npow lya-tot 11
    zbin_msr_matrix.T[11] = ((zbin_msr_matrix.T[4]*zbin_msr_matrix.T[0])**(-1)*
                              zbin_msr_matrix.T[10] * np.pi / (opt.kmax-opt.kmin))
    #print('1',zbin_msr_matrix.T[10])

    mf_output_df = pd.DataFrame(zbin_msr_matrix)
    mf_output_df.columns = mf_msrmnts

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
        if bin_zab[i] in mf_output_df.z.values:
            za_idx = mf_output_df.z == bin_zab[i]
            ztot_idx = i
            mf_lyb[i] = mf_output_df.mf_tot[ztot_idx]/mf_output_df.mf_a[za_idx]
    mf_output_df['mf_b'] = mf_lyb

    mf_bootstrap[m] = mf_output_df

    # POWER SPECTRUM



    znk_matrix = np.zeros((opt.zbinlen,n_pk_msrmnts,opt.kbinlen)) #  7 zbins,6 measurements, 20 kbins
    #print("Pk")
    for zidx in range(opt.zbinlen):
        #opt.updt(opt.zbinlen, zidx)
        msrmnt_in_kbin = np.zeros((n_pk_msrmnts,opt.kbinlen))
        count_in_kbin = np.zeros((n_pk_msrmnts,opt.kbinlen))
        msrmnt_in_kbin[0] = opt.kbin_centers
        count_in_kbin[0] = np.ones_like(opt.kbin_centers)
        count_in_kbin[6] = np.ones_like(opt.kbin_centers)
        msrmnt_in_kbin[7] = np.ones_like(opt.kbin_centers) * opt.zbin_centers[zidx]
        count_in_kbin[7] = np.ones_like(opt.kbin_centers)

        for qidx in range(nqso):
            # LYA FOREST: P ALPHA ALPHA
            zpix_a = qso_arr[qidx].get_zpix(opt.lya_rest)
            zmask_a = qso_arr[qidx].get_zmask(forest=(opt.lya_min,opt.lya_max,opt.lya_rest),
                                            zpix=zpix_a,zidx=zidx,zedges=opt.zbin_edges,name=name)
            if np.sum(zmask_a)>opt.min_pix:
                kpix,pk = qso_arr[qidx].get_autopower(mf_output_df.mf_a[zidx],zmask_a)
                for kidx in range(opt.kbinlen):
                    npow = mf_output_df.npow_a.values[zidx]
                    kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
                    pk_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=pk,zmask=zmask_a,kmask=kmask,
                                                       corr_tag=tag,npow=npow)
                    #znk_matrix[zidx][1,kidx] += np.sum(pk_sub)
                    #znk_matrix[zidx][6,kidx] += len(pk_sub)
                    msrmnt_in_kbin[1,kidx] += np.sum(pk_sub) #Paa
                    #msrmnt_in_kbin[6,kidx] += len(pk_sub)
                    count_in_kbin[1,kidx] += len(pk_sub) #num is Paa

            # LYB FOREST: P TOTAL TOTAL
            zpix_tot = qso_arr[qidx].get_zpix(opt.lyb_rest)
            zmask_tot = qso_arr[qidx].get_zmask(forest=(opt.lyb_min,opt.lyb_max,opt.lyb_rest),
                                            zpix=zpix_tot,zidx=zidx,zedges=opt.zbin_edges,name=name)

            if (np.sum(zmask_tot)>opt.min_pix):

                kpix,pk = qso_arr[qidx].get_autopower(mf_output_df.mf_tot[zidx],zmask_tot)
                for kidx in range(opt.kbinlen):
                    npow = mf_output_df.npow_tot.values[zidx]
                    kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
                    pk_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=pk,zmask=zmask_tot,kmask=kmask,
                                                       corr_tag=tag,npow=npow)
                    msrmnt_in_kbin[2,kidx] += np.nansum(pk_sub)
                    count_in_kbin[2,kidx] += len(pk_sub)

            #Cross power
            if (np.sum(zmask_a)>opt.min_pix)&(np.sum(zmask_tot)>opt.min_pix):
                kpix,pab,qab,dlam,res = qso_arr[qidx].cross_pk_fft(mask_lya=zmask_a,mask_lyb=zmask_tot,
                                      mf_lya=mf_output_df.mf_a[zidx],
                                      mf_lyb=mf_output_df.mf_tot[zidx])

                npow = mf_output_df.npow_atot.values[zidx]
                for kidx in range(opt.kbinlen):
                    kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
                    pab_sub,qab_sub = qso_arr[qidx].get_xpk_subsets(kpix,pab,qab,dlam,res,tag,npow,kmask)
                    #msrmnt_in_kbin[2,kidx] += np.nansum(pk_sub)
                    #count_in_kbin[2,kidx] += len(pk_sub)
                    # msrmnt_in_kbin[3,kidx] += np.nansum(pab_sub) #remove!
                    #print(np.nansum(pab_sub))
                    msrmnt_in_kbin[3,kidx] += np.nansum(pab_sub)
                    count_in_kbin[3,kidx] += len(pab_sub)
                    msrmnt_in_kbin[4,kidx] += np.sum(qab_sub) #remove!
                    count_in_kbin[4,kidx] += len(qab_sub)
                    msrmnt_in_kbin[6,kidx] += len(pab_sub)
        znk_matrix[zidx] = msrmnt_in_kbin/count_in_kbin

    #opt.updt(opt.zbinlen, opt.zbinlen)
    #print("Done!\n")

    # Finding Lyman beta power
    for i in range(len_zab):
        if bin_zab[i] in opt.zbin_centers:

            za_idx = np.where(opt.zbin_centers == bin_zab[i])[0][0]
            znk_matrix[i][5] = znk_matrix[i][2]-znk_matrix[za_idx][1]

    # Making 3d pk matrix into 2d pk data frame
    x = pd.DataFrame(znk_matrix[0].T,columns=pk_msrmnts)
    for i in range(1,opt.zbinlen):
        x = x.append(pd.DataFrame(znk_matrix[i].T,columns=pk_msrmnts))
    #print(np.shape(x))
    pk_bootstrap[m] = x


mf_bootstrap_2d = np.concatenate(mf_bootstrap.T)
pk_bootstrap_2d = np.concatenate(pk_bootstrap.T)
#np.reshape(mf_bootstrap.T,(M,opt.zbinlen*n_mf_msrmnts))
if inis.save_boot_mf:
    np.savetxt(inis.save_boot_mf_path,mf_bootstrap_2d)#).to_csv(inis.save_boot_mf_path, index=False)

if inis.save_boot_pk:
    np.savetxt(inis.save_boot_pk_path,pk_bootstrap_2d)

# print(pk_bootstrap[m])
opt.updt(M, M)




# PLOTTING ROUTINE
# import matplotlib.pyplot as plt
# qwer = mf_bootstrap.T[0]
# # Add a colorbar
# fig,ax = plt.subplots(1, 2,gridspec_kw={'width_ratios': [1,1]})
# fig.set_size_inches(12,8)
#
# im = ax[0].imshow(np.corrcoef(qwer), cmap = plt.cm.jet)
# fig.colorbar(im, ax=ax[0])
# # set the color limits - not necessary here, but good to know how.
# im.set_clim(0.0, 1.0)
# #plt.show()
#
# x = [np.var(qwer.T[:i].T) for i in range(1,M)]
# ax[1].plot(x)
# med = np.median(x)
# ax[1].set_ylim(med - med*0.01,med + med*0.01)
# fig.tight_layout()
# fig.savefig('../plot/figures/test.pdf')
# plt.clf()











# plt.show()
# print()
# qwer = np.reshape(np.concatenate(mf_bootstrap.T,axis=1)[0],(M,n_mf_msrmnts))
# print("HERE")
# print(qwer)
# print(qwer)
# print(np.concatenate(mf_bootstrap.T,axis=1)[0])
# qwer = np.reshape(np.concatenate(mf_bootstrap.T,axis=1)[0],(M,n_mf_msrmnts))
# # print(qwer)
# print(np.corrcoef(qwer))
# import matplotlib.pyplot as plt
# plt.imshow(np.corrcoef(qwer))
# plt.show()




# print()
# print(np.nanmean(mf_bootstrap.T,axis=2))
# print(np.nanvar(mf_bootstrap.T,axis=2))
# print(np.corrcoef(np.reshape(np.concatenate(mf_bootstrap,axis=1)[0],(M,n_mf_msrmnts))))
# print(np.mean(mf_bootstrap.T[1],axis=1))/
# for i in range(n_mf_msrmnts):
#     #print(mf_bootstrap[0].T[i])
#     print(np.nanmean(mf_bootstrap[0].T[i]))
#     print(np.nanvar(mf_bootstrap[0].T[i]))
#
# print(np.nanmean(mf_bootstrap[0].T,axis=1))

# print('0')
# print(mf_bootstrap[0].T)
# print("1")
# print(mf_bootstrap[0].T[0])
# print("2")
# print(np.nanmean(mf_bootstrap[0].T[0]))
# print("3")
# print(np.nanmean(mf_bootstrap[0],axis=1))
# print("4")

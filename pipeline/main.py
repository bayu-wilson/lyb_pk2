#!/usr/bin/env python

from QuasarSpectrum import QuasarSpectrum
import inis
import numpy as np
import pandas as pd
import options as opt
from scipy.interpolate import griddata
import copy

cat_name = inis.cat_name
tag = inis.tag
rescale_flux = inis.rescale_flux #0.1 #0.1#0.1 #0.1 #0.0 #0.1
print(cat_name)
print(tag)

QuasarSpectrum.load_cat(cat_name)
nqso = QuasarSpectrum.nqso
print("Loading Data")
qso_arr = []

for i in range(nqso):
    q = QuasarSpectrum.load_qso_data(i,tag=tag,rescale_flux=rescale_flux)
    if ("noB" in tag)&(inis.add_beta): # adding Lyb
        q.get_new_forest(rescale_flux=rescale_flux,wrange =
            (opt.lyb_min,opt.lyb_max,opt.lyb_rest, opt.xs_beta))
    if inis.cat_name.startswith('mocks')&(inis.add_ovi): # adding OVI
        q.get_new_forest(rescale_flux=rescale_flux,wrange =
            (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d1, opt.xs_ovi))
        q.get_new_forest(rescale_flux=rescale_flux,wrange =
            (opt.ovi_min,opt.ovi_max,opt.ovi_rest_d2, opt.xs_ovi))
    if inis.cat_name.startswith('mocks')&(inis.add_sithree): #adding SiIII
        q.get_new_forest(rescale_flux=rescale_flux,wrange =
            (opt.sithree_min,opt.sithree_max,opt.sithree_rest_d1, opt.xs_sithree))
        # q.get_new_forest(rescale_flux=rescale_flux,wrange =
        #     (opt.sithree_min,opt.sithree_max,opt.sithree_rest_d2, opt.xs_sithree))
    qso_arr.append(q)

print("Done!\n")
print("nchunks: ", len(qso_arr))

mf_msrmnts = opt.mf_msrmnts
n_mf_msrmnts = len(mf_msrmnts)
zbin_msr_matrix = np.zeros((opt.zbinlen,n_mf_msrmnts))
# flux_pdf = [[] for tmp in range(opt.zbinlen)]
# pdf_bins = np.arange(-0.025,1.05,.05)
print("mf")
f = open(inis.save_kzq_mf_path,'w') #july10
for zidx in range(opt.zbinlen):
    msrmnt_in_zbin =  np.zeros(n_mf_msrmnts)
    count_in_zbin = np.zeros(n_mf_msrmnts)
    opt.updt(opt.zbinlen, zidx)

    zbin_msrmnt = [[] for idx in range(n_mf_msrmnts)]
    for i in range(nqso):
        name = qso_arr[i].name
        zpix_a = qso_arr[i].get_zpix(opt.lya_rest)
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

        s = "{0:g} {1:g} {2:g} {3:g} {4:g} {5:g}".format(i,opt.zbin_centers[zidx],np.sum(qso_arr[i].flux[mask]), np.sum(mask),np.sum(qso_arr[i].flux[mask_b]),np.sum(mask_b)) #july10
        #"{4:g} 0 {0:g} {1:g} {2:g} {3:g}".format(opt.zbin_centers[zidx],opt.kbin_centers[kidx],
        #                                np.sum(pk_sub),len(pk_sub),qidx) #july1 #average of pab_sub??
        f.write(s+'\n') #july10
    zbin_msr_matrix[zidx] = msrmnt_in_zbin/count_in_zbin
opt.updt(opt.zbinlen, opt.zbinlen)
f.close() #july10
print("Done!\n")

zbin_msr_matrix.T[3] = list(QuasarSpectrum.get_npow(mf=zbin_msr_matrix.T[0], # npow alpha 3
                                        nvar=zbin_msr_matrix.T[1],
                                        dloglambda=zbin_msr_matrix.T[2]))
zbin_msr_matrix.T[7] = list(QuasarSpectrum.get_npow(mf=zbin_msr_matrix.T[4], # npow total 7
                                        nvar=zbin_msr_matrix.T[5],
                                        dloglambda=zbin_msr_matrix.T[6]))
zbin_msr_matrix.T[8] = opt.zbin_centers # zbins 8
zbin_msr_matrix.T[11] = ((zbin_msr_matrix.T[4]*zbin_msr_matrix.T[0])**(-1)* # npow lya-tot 11
                          zbin_msr_matrix.T[10] * np.pi / (opt.kmax-opt.kmin))

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

### CONTINUUM CORRECTION
# F_true = F_est * (1-deltaC/C_true)
# print(mf_output_df['mf_a'])
if inis.continuum_correction:
    C_a_ratio = opt.continuum_correction(opt.zbin_centers)
    C_t_ratio = opt.continuum_correction(opt.zbin_centers-0.2)
    C_b_ratio = np.abs(np.abs(C_a_ratio)-np.abs(C_t_ratio))
    mf_output_df['mf_a'] = mf_output_df['mf_a']*(1-C_a_ratio)
    mf_output_df['mf_tot'] = mf_output_df['mf_tot']*(1-C_t_ratio)
    mf_output_df['mf_b'] = mf_output_df['mf_b']*(1-C_b_ratio)



# print(mf_output_df['mf_a'])

# mf_output_df.mf_a  = mf_output_df.mf_a*0 + 1.0 #aug15 set mf to zero for test for matt
# mf_output_df.mf_tot  = mf_output_df.mf_tot*0 + 1.0 #aug15 set mf to zero for test for matt
# mf_output_df.mf_b  = mf_output_df.mf_b*0 + 1.0 #aug15 set mf to zero for test for matt

# mf_output_df = pd.read_csv("../output/test_dla_jun20.csv")

################################################################################
################        Power Spectrum          ################################
################################################################################

pk_msrmnts = opt.pk_msrmnts
n_pk_msrmnts = len(pk_msrmnts)

znk_matrix = np.zeros((opt.zbinlen,n_pk_msrmnts,opt.kbinlen)) #  7 zbins,6 measurements, 20 kbins
print("Pk")
f = open(inis.save_kzq_pk_path,'w') #"../output/qsos_pk.txt",'w') #july1
for zidx in range(opt.zbinlen): #aug15 #aug20
    opt.updt(opt.zbinlen, zidx)
    msrmnt_in_kbin = np.zeros((n_pk_msrmnts,opt.kbinlen))
    count_in_kbin = np.zeros((n_pk_msrmnts,opt.kbinlen))
    msrmnt_in_kbin[0] = opt.kbin_centers
    count_in_kbin[0] = np.ones_like(opt.kbin_centers) # these have to be ones.I just divide by one later.
    count_in_kbin[6] = np.ones_like(opt.kbin_centers)
    count_in_kbin[7] = np.ones_like(opt.kbin_centers)
    count_in_kbin[8] = np.ones_like(opt.kbin_centers)
    msrmnt_in_kbin[-1] = np.ones_like(opt.kbin_centers) * opt.zbin_centers[zidx]
    count_in_kbin[-1] = np.ones_like(opt.kbin_centers)
    #if zidx == 6: # Checking why z4.2 sucks. 6/13/19
    #    f_big = open("N_cross.txt",'w')
    for qidx in range(nqso): #aug20
        #################### LYA FOREST: P ALPHA ALPHA ####################
        name = qso_arr[qidx].name
        zpix_a = qso_arr[qidx].get_zpix(opt.lya_rest)
        zmask_a = qso_arr[qidx].get_zmask(forest=(opt.lya_min,opt.lya_max,opt.lya_rest),
                                        zpix=zpix_a,zidx=zidx,zedges=opt.zbin_edges,name=name)
        zpix_tot = qso_arr[qidx].get_zpix(opt.lyb_rest)
        zmask_tot = qso_arr[qidx].get_zmask(forest=(opt.lyb_min,opt.lyb_max,opt.lyb_rest),
                                        zpix=zpix_tot,zidx=zidx,zedges=opt.zbin_edges,name=name)
        nchunks_a = opt.how_many_chunks(zmask_a)
        nchunks_tot = opt.how_many_chunks(zmask_tot)

        # if zidx == 2: #aug20
        #     debug = np.column_stack((qso_arr[qidx].wavelength,qso_arr[qidx].flux,#aug20
        #                             zmask_a,zmask_tot))#aug20
        #
        #     np.savetxt("../output/qsos/q{:02d}_z3.8_masks_full_forest.txt".format(qidx),debug) #ff #aug20
            #np.savetxt("../output/qsos/q{:02d}_z3.8_masks_match_lyb.txt".format(qidx),debug)


        # SWITCH #aug14
        #nchunks_a=2 # TEST CHANGE THIS
        for chunk_idx_a in range(nchunks_a):
            #if qidx == 2999: #july 17
            #    print(np.sum(zmask_a),np.where(zmask_a)[0],len(zmask_a),chunk_idx_a,opt.get_chunks(zmask_a),nchunks_a)
            zmask_a_sub = opt.get_chunks(zmask_a)[chunk_idx_a]
            # if (qidx == 64)&(zidx==2):
            #     #zmask_a_sub = zmask_a_sub[23:]
            #     #zmask_a_sub = zmask_a_sub[22:]
            #     print("\n remove dla?: ", inis.remove_dla)
            #     print("nchunks: ", nchunks_a)
            #     print("index chunk: ",chunk_idx_a)
            #     print("npix: ",len(zmask_a_sub))
            #     print(zmask_a_sub)


            if np.sum(zmask_a_sub)>opt.min_pix:
                kpix,pk = qso_arr[qidx].get_autopower(mf_output_df.mf_a[zidx],zmask_a_sub)
                for kidx in range(opt.kbinlen): #aug15 #aug20
                    npow = mf_output_df.npow_a.values[zidx]
                    kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
                    pk_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=pk,zmask=zmask_a_sub,kmask=kmask,
                                                       corr_tag=tag,npow=npow)
                    msrmnt_in_kbin[1,kidx] += np.sum(pk_sub) #Paa
                    count_in_kbin[1,kidx] += len(pk_sub) #num is Paa
                    msrmnt_in_kbin[6,kidx] += len(pk_sub) # npix_aa

                    s = "{4:g} 0 {0:g} {1:g} {2:g} {3:g}".format(opt.zbin_centers[zidx],opt.kbin_centers[kidx],
                                                    np.sum(pk_sub),len(pk_sub),qidx) #july1 #average of pab_sub??
                    f.write(s+'\n') #july1

        # SWITCH #aug14
        # Loop through chunks again (using zmask_tot)
        # LYB FOREST: P TOTAL TOTAL
        for chunk_idx_tot in range(nchunks_tot):
            zmask_tot_sub = opt.get_chunks(zmask_tot)[chunk_idx_tot]
            #if (qidx == 60)&(zidx==4):
                #zmask_a_sub = zmask_a_sub[23:]
                #zmask_a_sub = zmask_a_sub[22:]
                # print("\n remove dla?: ", inis.remove_dla)
                # print("nchunks: ", nchunks_tot)
                # print("index chunk: ",chunk_idx_tot)
                # print("npix: ",len(zmask_tot_sub))
                # print(zmask_tot_sub)


            if (np.sum(zmask_tot_sub)>opt.min_pix):
                kpix,pk = qso_arr[qidx].get_autopower(mf_output_df.mf_tot[zidx],zmask_tot_sub)
                for kidx in range(opt.kbinlen): #aug15 #aug20
                    npow = mf_output_df.npow_tot.values[zidx]
                    kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
                    pk_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=pk,zmask=zmask_tot_sub,kmask=kmask,
                                                       corr_tag=tag,npow=npow)
                    msrmnt_in_kbin[2,kidx] += np.sum(pk_sub) #aug15
                    count_in_kbin[2,kidx] += len(pk_sub)
                    msrmnt_in_kbin[7,kidx] += len(pk_sub) #npix_tt #sept25

                    s = "{4:g} 1 {0:g} {1:g} {2:g} {3:g}".format(opt.zbin_centers[zidx],opt.kbin_centers[kidx],
                                                    np.sum(pk_sub),len(pk_sub), qidx) #july1 #average of pab_sub??
                    f.write(s+'\n') #july1


        # Loop through more chunks here. Check zmask_tot and zmask_a CROSS POWER
        idx_a = 0
        idx_tot = 0
        while (idx_tot < nchunks_tot)&(idx_a < nchunks_a):
            if (nchunks_a>0)&(nchunks_tot>0):
                ztot = zpix_tot[opt.get_chunks(zmask_tot)[idx_tot]]
                za = zpix_a[opt.get_chunks(zmask_a)[idx_a]]
                ztot_min = np.min(ztot)
                ztot_max = np.max(ztot)
                za_min = np.min(za)
                za_max = np.max(za)
                mask_chk_a = (zpix_a>np.max([ztot_min,za_min]))&(zpix_a<np.min([ztot_max,za_max]))
                mask_chk_tot = (zpix_tot>np.max([ztot_min,za_min]))&(zpix_tot<np.min([ztot_max,za_max]))
                if (np.sum(mask_chk_a)>opt.min_pix)&(np.sum(mask_chk_tot)>opt.min_pix):
                    kpix,pab,qab,dlam,resa,resb = qso_arr[qidx].cross_pk_fft(mask_lya=mask_chk_a,mask_lyb=mask_chk_tot,
                                          mf_lya=mf_output_df.mf_a[zidx],
                                          mf_lyb=mf_output_df.mf_tot[zidx])
                    npow = mf_output_df.npow_atot.values[zidx]
                    for kidx in range(opt.kbinlen): #aug20
                        kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
                        pab_sub,qab_sub = qso_arr[qidx].get_xpk_subsets(kpix,pab,qab,dlam,resa,resb,tag,npow,kmask)
                        ##########################################
                        # HERE IS WHERE I NEED TO WRITE OUT POWER VALUES TO ANOTHER FILE!!

                        #for i in range(len(qso_arr[qidx].wavelength[zmask_a])): #july1
                        s = "{4:g} 2 {0:g} {1:g} {2:g} {3:g}".format(opt.zbin_centers[zidx],opt.kbin_centers[kidx],
                                                        np.sum(pab_sub),len(pab_sub), qidx) #july1 #average of pab_sub??
                        f.write(s+'\n') #july1
                        ##########################################

                        msrmnt_in_kbin[3,kidx] += np.sum(pab_sub)
                        count_in_kbin[3,kidx] += len(pab_sub)
                        msrmnt_in_kbin[4,kidx] += np.sum(qab_sub)
                        count_in_kbin[4,kidx] += len(qab_sub)
                        msrmnt_in_kbin[8,kidx] += len(pab_sub) #npix_ab #sept25
                if za_max<ztot_min: # no overlap
                    idx_a +=1
                    idx_tot +=0
                elif ztot_max<za_min:
                    idx_a +=0
                    idx_tot +=1
                else:
                    idx_tot +=1

            else:
                break
    znk_matrix[zidx] = msrmnt_in_kbin/count_in_kbin
f.close() #july1
# f_new.close() #july31

# f_big.close() # Checking why z4.2 sucks. 6/13/19
opt.updt(opt.zbinlen, opt.zbinlen)
print("Done!\n")

# Finding Lyman beta power
for i in range(len_zab):
    if bin_zab[i] in opt.zbin_centers:
        za_idx = np.where(opt.zbin_centers == bin_zab[i])[0][0]
        znk_matrix[i][5] = znk_matrix[i][2]-znk_matrix[za_idx][1]

# Making 3d pk matrix into 2d pk data frame
x = pd.DataFrame(znk_matrix[0].T,columns=pk_msrmnts)
for i in range(1,opt.zbinlen):
    x = x.append(pd.DataFrame(znk_matrix[i].T,columns=pk_msrmnts))


########## SUBTRACTING METAL POWER ##########
if inis.subtract_metal_power:
    metal_power = np.concatenate([np.loadtxt('../data/obs/pk_xs_avg.txt')]*opt.zbinlen) #subtracting metals again sept25

    x.Paa = x.Paa.values-metal_power
    x.Ptot = x.Ptot.values-metal_power
    #x.Pbb = x.Pbb.values-metal_power
#############################################

# Saving routine
if inis.save_mf:
    mf_output_df.to_csv(inis.save_mf_path,index=False)

if inis.save_pk:
    x.to_csv(inis.save_pk_path,index=False)
if inis.save_dla:
    x.to_csv(inis.save_dla_path,index=False)

print("OVI Added?: ", inis.add_ovi)
print("SiIII Added?: ", inis.add_sithree)
if (inis.remove_dla)&inis.cat_name.startswith('obs'):
    print("DLA's removed?: ", inis.remove_dla)
else:
    print("DLA's removed?: False")
print("Metals Subtracted?: ", inis.subtract_metal_power) #see line 310
print("Continuum Corrected?: ", inis.continuum_correction)
print("saved pre-boot mf file here: ", inis.save_kzq_mf_path)
print("saved pre-boot pk file here: ", inis.save_kzq_pk_path)
print("saved mf here: ", inis.save_mf_path)
print("saved pf here: ", inis.save_pk_path)
print("---- Script complete ----")


# print(qso_arr[99].resolution)
    #print(i)
#     if inis.remove_dla:
#         #DLAs
#         m = q.mask_dla
#         if opt.how_many_chunks(m) == 0:
#             qso_arr.append(q)
#         else:
#             #print(i)
#             for j in range(opt.how_many_chunks(m)):
#                 q_test = copy.deepcopy(q)
#                 q_test.get_chunks(j)
#                 qso_arr.append(q_test)
#     else:
#         qso_arr.append(q)
# if inis.remove_dla:
#     nqso = len(qso_arr)

                    #print(idx_a)
                #if (ztot_max<za_max):
                # idx_tot +=1
                    #print(idx_tot)
                    #print(idx_tot)
                #if np.min([ztot_max,za_max]) == za_max:
                #    idx_tot+=1
                #    print(idx_tot)
                #if np.min([ztot_max,za_max]) == za_max:
                #    idx_a +=1

        # if nchunks_a == nchunks_tot: # DONT DELETE
        #     for chunk_idx_x in range(nchunks_a):
        #         zmask_a_sub = opt.get_chunks(zmask_a)[chunk_idx_x]
        #         zmask_tot_sub = opt.get_chunks(zmask_tot)[chunk_idx_x]
        #         if (np.sum(zmask_a_sub)>opt.min_pix)&(np.sum(zmask_tot_sub)>opt.min_pix):
        #             kpix,pab,qab,dlam,res = qso_arr[qidx].cross_pk_fft(mask_lya=zmask_a_sub,mask_lyb=zmask_tot_sub,
        #                                   mf_lya=mf_output_df.mf_a[zidx],
        #                                   mf_lyb=mf_output_df.mf_tot[zidx])
        #             npow = mf_output_df.npow_atot.values[zidx]
        #             for kidx in range(opt.kbinlen):
        #                 kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
        #                 pab_sub,qab_sub = qso_arr[qidx].get_xpk_subsets(kpix,pab,qab,dlam,res,tag,npow,kmask)
        #                 msrmnt_in_kbin[3,kidx] += np.nansum(pab_sub)
        #                 count_in_kbin[3,kidx] += len(pab_sub)
        #                 msrmnt_in_kbin[4,kidx] += np.sum(qab_sub)
        #                 count_in_kbin[4,kidx] += len(qab_sub)
        #                 msrmnt_in_kbin[7,kidx] += len(pab_sub) #npix_ab

        # if nchunks_a != nchunks_tot: # DONT DELETE
        #     min_nchunk = np.min(nchunks_a,nchunks_tot)
        #     for chunk_idx_x in range(min_nchunk):
        #         zmask_a_sub = opt.get_chunks(zmask_a)[chunk_idx_x]
        #         zmask_tot_sub = opt.get_chunks(zmask_tot)[chunk_idx_x]
        #         if (np.sum(zmask_a_sub)>opt.min_pix)&(np.sum(zmask_tot_sub)>opt.min_pix):
        #             kpix,pab,qab,dlam,res = qso_arr[qidx].cross_pk_fft(mask_lya=zmask_a_sub,mask_lyb=zmask_tot_sub,
        #                                   mf_lya=mf_output_df.mf_a[zidx],
        #                                   mf_lyb=mf_output_df.mf_tot[zidx])
        #             npow = mf_output_df.npow_atot.values[zidx]
        #             for kidx in range(opt.kbinlen):
        #                 kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
        #                 pab_sub,qab_sub = qso_arr[qidx].get_xpk_subsets(kpix,pab,qab,dlam,res,tag,npow,kmask)
        #                 msrmnt_in_kbin[3,kidx] += np.nansum(pab_sub)
        #                 count_in_kbin[3,kidx] += len(pab_sub)
        #                 msrmnt_in_kbin[4,kidx] += np.sum(qab_sub)
        #                 count_in_kbin[4,kidx] += len(qab_sub)
        #                 msrmnt_in_kbin[7,kidx] += len(pab_sub) #npix_ab



# (mf_output_df.mf_tot*mf_output_df.mf_a)**(-1)*np.sum()
# print(mf_output_df)
# print(x)
# print(x[x.z==4.2])


# test_dlambda = np.log10(qso_arr[0].wavelength[2])-np.log10(qso_arr[0].wavelength[1])
# print(test_dlambda)
# print(np.log10(test_dlambda))
# print(qso_arr[0].ferr)
# sys.exit()

        #kmask = qso_arr[qidx].get_kmask(kpix=kpix,kidx=kidx,kedges=opt.kbin_edges)
        #    if (('corrNR' in tag) or
        #         ('corrR' in tag) or
        #         ('wR' in tag) or
        #         ('wNR' in tag)):
        #         #dloglam = np.median(qso_arr[qidx].dloglambda)
        #         dv = opt.c_kms*np.log(10)*dlam[kmask]
        #         res_new = res[kmask]
        #         #R = np.median(qso_arr[qidx].resolution)
        #         pab_sub,qab_sub =pab[kmask],qab[kmask]
        #         pab_sub = pab_sub/opt.window(k=kpix[kmask],p=dv,R=res_new)**2
        #         qab_sub = qab_sub/opt.window(k=kpix[kmask],p=dv,R=res_new)**2
        #         #print(opt.window(k=kpix[kmask],p=dv,R=res_new)**2)
        #     if (('corrNR' in tag) or
        #         ('corrN' in tag) or
        #         ('wN' in tag) or
        #         ('wNR' in tag)):
        #             a
        #
        #     else:
        #         pab_sub,qab_sub =pab[kmask],qab[kmask]


            #pab_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=pab,zmask=zmask_tot,kmask=kmask,
            #                                   corr_tag=tag,npow=None)#npow)
            #print(len(kmask),len(zmask_tot),len(pab),len(kpix))
            #sys.exit()
            #qab_sub = qso_arr[qidx].get_pk_subsets(kpix=kpix,pk=qab,zmask=zmask_tot,kmask=kmask,
            #                                   corr_tag=tag,npow=npow)
            #qab_sub =qab[kmask]


            # msrmnt_in_kbin[3,kidx] += np.nansum(pab_sub) #remove!
            # #print(np.nansum(pab_sub))
            # count_in_kbin[3,kidx] += len(pab_sub)
            # msrmnt_in_kbin[4,kidx] += np.sum(qab_sub) #remove!
            # count_in_kbin[4,kidx] += len(qab_sub)
            #
            # msrmnt_in_kbin[6,kidx] += len(pab_sub)




            # Checking if `how_many_chunks` works and if chunks are okay 6/19/19
            # if opt.how_many_chunks(zmask_a)>1: #mask is basically mask_a
            #     # print(zmask_a[zmask_a][0])
            #     # print("{0}".format(zmask_a[zmask_a][0]))
            #     f = open("test{0}{1}.txt".format(qidx,zidx),'w')
            #     for i in range(len(qso_arr[qidx].wavelength)):
            #         s = "{0:g} {1:g} {2}".format(qso_arr[qidx].wavelength[i],
            #                                     qso_arr[qidx].flux[i],
            #                                     zmask_a[i])
            #         #print(zmask_a[zmask_a][i])

            #     f.write(s+'\n')
            # f.close()
            # np.save
            # import sys
            # sys.exit()
            # print(i,zidx)


            # TEST_NAME = "J0234-1806" # Checking why z4.2 sucks. 6/13/19
            #chunks_a[chnk_idx]
            #     if (qso_arr[qidx].name == TEST_NAME) & (zidx==6):
            #         f = open("test.txt",'w')
            #         for i in range(len(qso_arr[qidx].wavelength[zmask_a])):
            #             s = "{0:g} {1:g}".format(qso_arr[qidx].wavelength[zmask_a][i],
            #                                         qso_arr[qidx].flux[zmask_a][i])
            #             f.write(s+'\n')
            #         f.close()




        # if (qso_arr[qidx].name == TEST_NAME) & (zidx==6): # Checking why z4.2 sucks. 6/13/19
        #     print(opt.how_many_chunks(qso_arr[qidx].mask_dla))
        #     print(inis.remove_dla, np.sum(zmask_a),np.sum(zmask_tot))
        #Cross power
            #print(np.sum(zmask_a),np.sum(zmask_tot))
            # if (qso_arr[qidx].name == TEST_NAME) & (zidx==6): # Checking why z4.2 sucks. 6/13/19
            #     print("asdfasdf")
            #     if inis.remove_dla:
            #         f = open("test_cross_remove_dla.txt",'w')
            #     else:
            #         f = open("test_cross_keep_dla.txt",'w')
            #     for i in range(len(qso_arr[qidx].wavelength[zmask_tot])):
            #         s = "{0:g} {1:g}".format(qso_arr[qidx].wavelength[zmask_tot][i],
            #                                     qso_arr[qidx].flux[zmask_tot][i])
            #         f.write(s+'\n')
            #     f.close()

                #print(np.array(za)[0])
                # print(ztot)
                # print(np.sum(mask_chk_a))
                # print(np.sum(mask_chk_tot))
                # f = open("test_a.txt",'w')
                # for i in range(len(za)):
                #     #print(za[i],mask_chk_a[i])
                #     s = "{0} {1}".format(za[i],qso_arr[qidx].flux[opt.get_chunks(zmask_a)[idx_a]][i])
                #     f.write(s+'\n')
                # f.close()
                # f = open("test_tot.txt",'w')
                # for i in range(len(ztot)):
                #     #print(ztot[i],mask_chk_tot[i])
                #     s = "{0} {1}".format(ztot[i],qso_arr[qidx].flux[opt.get_chunks(zmask_tot)[idx_tot]][i])
                #     f.write(s+'\n')
                # f.close()
                # import sys
                # sys.exit()

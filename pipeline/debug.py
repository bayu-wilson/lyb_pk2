#!/usr/bin/env python

import numpy as np
import inis
import options as opt
import matplotlib.pyplot as plt
import pandas as pd
import glob

qso_step = 1

qsos_ffi = glob.glob("../output/qsos/q*_z3.8_masks_full_forest.txt")
qsos_mbi = glob.glob("../output/qsos/q*_z3.8_masks_match_lyb.txt")
qsos_ffii = glob.glob("../output/qsos/q*_z4.0_masks_full_forest.txt")
qsos_mbii = glob.glob("../output/qsos/q*_z4.0_masks_match_lyb.txt")
qsos_ffi.sort(),qsos_mbi.sort(),qsos_ffii.sort(),qsos_mbii.sort()

def get_autopower(flux,mf):
    dloglambda = 1e-5
    rf_x = flux/mf - 1
    dv = opt.c_kms*dloglambda*np.log(10)
    rf_k = np.fft.fft(rf_x)*dv
    N = len(rf_k)
    V = N*dv
    dk = 2*np.pi/V
    k = dk*np.arange(0,N,1)
    pk = np.abs(rf_k)**2/V
    return k,pk
# print(np.loadtxt(qsos_ffi[0])[:,0])
#
# counter = 0
# fig,ax=plt.subplots(nrows=1,ncols=2)
# fig.set_size_inches(10,20)
mf_a_ff, npix_a_ff = 0,0
mf_t_ff, npix_t_ff = 0,0
mf_a_mb, npix_a_mb = 0,0
mf_t_mb, npix_t_mb = 0,0
for qidx in range(0,100,qso_step):
    t_ff = np.loadtxt(qsos_ffii[qidx])
    t_mb = np.loadtxt(qsos_mbii[qidx])
    # ax[0].plot(t_ff[:,0],t_ff[:,1]+counter, color='blue')
    # ax[1].plot(t_mb[:,0],t_mb[:,1]+counter, color='blue')
    mask_a_ff = np.array(t_ff[:,2],dtype='bool')
    mask_t_ff = np.array(t_ff[:,3],dtype='bool')
    mask_a_mb = np.array(t_mb[:,2],dtype='bool')
    mask_t_mb = np.array(t_mb[:,3],dtype='bool')

    mf_a_ff   += np.sum(t_ff[:,1][mask_a_ff])
    npix_a_ff += np.sum(mask_a_ff)

    mf_t_ff   += np.sum(t_ff[:,1][mask_t_ff])
    npix_t_ff += np.sum(mask_t_ff)

    mf_a_mb   += np.sum(t_mb[:,1][mask_a_mb])
    npix_a_mb += np.sum(mask_a_mb)

    mf_t_mb   += np.sum(t_mb[:,1][mask_t_mb])
    npix_t_mb += np.sum(mask_t_mb)
mf_a_ff = mf_a_ff/npix_a_ff
mf_t_ff = mf_t_ff/npix_t_ff
mf_a_mb = mf_a_mb/npix_a_mb
mf_t_mb = mf_t_mb/npix_t_mb
print(mf_a_ff,mf_t_ff,mf_a_mb,mf_t_mb)

pk_a_ff, npix_pk_a_ff = (np.zeros(opt.kbinlen),np.zeros(opt.kbinlen))
pk_t_ff, npix_pk_t_ff = (np.zeros(opt.kbinlen),np.zeros(opt.kbinlen))
pk_a_mb, npix_pk_a_mb = (np.zeros(opt.kbinlen),np.zeros(opt.kbinlen))
pk_t_mb, npix_pk_t_mb = (np.zeros(opt.kbinlen),np.zeros(opt.kbinlen))
for qidx in range(0,100,qso_step):
    t_ff = np.loadtxt(qsos_ffi[qidx])
    t_mb = np.loadtxt(qsos_mbi[qidx])
    mask_a_ff = np.array(t_ff[:,2],dtype='bool')
    mask_t_ff = np.array(t_ff[:,3],dtype='bool')
    mask_a_mb = np.array(t_mb[:,2],dtype='bool')
    mask_t_mb = np.array(t_mb[:,3],dtype='bool')
    if np.sum(mask_a_ff)>0:
        kpix,pk = get_autopower(t_ff[:,1][mask_a_ff],mf_a_ff)
        for kidx in range(opt.kbinlen):
            mask = ((kpix>=opt.kbin_edges[kidx])&(kpix<opt.kbin_edges[kidx+1]))
            pk_a_ff[kidx] += np.sum(pk[mask])
            npix_pk_a_ff[kidx] += np.sum(mask)
    if np.sum(mask_a_mb)>0:
        kpix,pk = get_autopower(t_mb[:,1][mask_a_mb],mf_a_mb)
        for kidx in range(opt.kbinlen):
            mask = ((kpix>=opt.kbin_edges[kidx])&(kpix<opt.kbin_edges[kidx+1]))
            pk_a_mb[kidx] += np.sum(pk[mask])
            npix_pk_a_mb[kidx] += np.sum(mask)
    if np.sum(mask_t_ff)>0:
        kpix,pk = get_autopower(t_ff[:,1][mask_t_ff],mf_t_ff)
        for kidx in range(opt.kbinlen):
            mask = ((kpix>=opt.kbin_edges[kidx])&(kpix<opt.kbin_edges[kidx+1]))
            pk_t_ff[kidx] += np.sum(pk[mask])
            npix_pk_t_ff[kidx] += np.sum(mask)
    if np.sum(mask_t_mb)>0:
        kpix,pk = get_autopower(t_ff[:,1][mask_t_mb],mf_t_ff)
        for kidx in range(opt.kbinlen):
            mask = ((kpix>=opt.kbin_edges[kidx])&(kpix<opt.kbin_edges[kidx+1]))
            pk_t_mb[kidx] += np.sum(pk[mask])
            npix_pk_t_mb[kidx] += np.sum(mask)

print(pk_a_ff,npix_pk_a_ff)
pk_a_ff = pk_a_ff/npix_pk_a_ff
pk_t_ff = pk_t_ff/npix_pk_t_ff
pk_a_mb = pk_a_mb/npix_pk_a_mb
pk_t_mb = pk_t_mb/npix_pk_t_mb
kpix = opt.kbin_centers
k_x = np.log10(kpix)

plt.plot(k_x,kpix*pk_a_ff/np.pi,color='b')
plt.plot(k_x,kpix*pk_t_ff/np.pi,color='r')
plt.plot(k_x,kpix*pk_a_mb/np.pi,color='g')
plt.plot(k_x,kpix*pk_t_mb/np.pi,color='purple',ls='dotted')
plt.show()

    # ax[0].plot(t_ff[:,0][mask_a_ff],t_ff[:,1][mask_a_ff]+counter, color='red')
    # ax[0].plot(t_ff[:,0][mask_b_ff],t_ff[:,1][mask_b_ff]+counter, color='green')
    # ax[1].plot(t_mb[:,0][mask_a_mb],t_mb[:,1][mask_a_mb]+counter, color='gold')
    # ax[1].plot(t_mb[:,0][mask_b_mb],t_mb[:,1][mask_b_mb]+counter, color='purple')
    # counter+=1


    #print(np.sum(t_ff[:,2]),np.sum(t_mb[:,2]))

#     t = np.loadtxt(qsos_ffi[qidx])
#     ax.plot(t[:,0],t[:,1]+counter, color='blue')
#
# # print("1")
# plt.savefig("/Users/bayuwilson/Desktop/debug.pdf")
# # plt.show()
# plt.clf()

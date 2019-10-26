#!/usr/bin/env python

import inis
import options as opt
import numpy as np
import pandas as pd


mfdata = pd.read_csv(inis.save_mf_path)
pkdata = pd.read_csv(inis.save_pk_path)
if inis.mock_or_obs == "obs":
    mf = np.loadtxt(inis.save_boot_mf_path)
    err_mfa = np.cov(mf).diagonal()[0:opt.zbinlen]**0.5
    err_mft = np.cov(mf).diagonal()[opt.zbinlen:opt.zbinlen*2]**0.5
    err_mfb = opt.find_err_mf_beta(mfdata.mf_b.values,mfdata.mf_a.values,err_mfa,
                                                      mfdata.mf_tot.values,err_mft)
    # mf_bootstrap = np.reshape(mf,(12,opt.zbinlen,inis.M))
    # err_mfa = np.cov(mf_bootstrap[0]).diagonal()**0.5
    # err_mft = np.cov(mf_bootstrap[4]).diagonal()**0.5
    # err_mfb = np.cov(mf_bootstrap[9]).diagonal()**0.5
    if inis.save_mf_with_err:
        columns = ["z","mfa","err_mfa","mft","err_mft","mfb","err_mfb"]
        mf_everything = np.column_stack((mfdata.z,
                        mfdata.mf_a, err_mfa,
                        mfdata.mf_tot, err_mft,
                        mfdata.mf_b, err_mfb))
        df_meanflux = pd.DataFrame(mf_everything,columns=columns)
        #df_meanflux.mft[1] = np.nan
        df_meanflux.to_csv(inis.save_mf_with_err_path, index=False)
        #pd.savetxt(inis.save_mf_with_err_path,mf_everything)


    p = np.loadtxt(inis.save_boot_pk_path)
    N_kz = opt.zbinlen*opt.kbinlen
    pk_err_diag = np.cov(p).diagonal()**0.5
    err_paa = pk_err_diag[0:N_kz]
    err_ptt = pk_err_diag[N_kz:N_kz*2]
    err_pab = pk_err_diag[N_kz*2:N_kz*3]

    err_paa_sub = err_paa[0*opt.kbinlen:3*opt.kbinlen]
    err_ptt_sub = err_ptt[4*opt.kbinlen:7*opt.kbinlen]
    err_pbb_sub = np.sqrt(err_paa_sub**2+err_ptt_sub**2)
    err_pbb = np.ones(N_kz)*np.nan
    err_pbb[4*opt.kbinlen:7*opt.kbinlen] = err_pbb_sub
    # pk_bootstrap = np.reshape(p,(8,opt.kbinlen*opt.zbinlen,inis.M))
    # err_paa = np.reshape(np.cov(pk_bootstrap[1]).diagonal()**0.5,
    #          (opt.zbinlen,opt.kbinlen))
    # err_ptt = np.reshape(np.cov(pk_bootstrap[2]).diagonal()**0.5,
    #          (opt.zbinlen,opt.kbinlen))
    # err_pab = np.reshape(np.cov(pk_bootstrap[3]).diagonal()**0.5,
    #          (opt.zbinlen,opt.kbinlen))
    if inis.save_pk_with_err:
        columns = ['k','z','paa','err_paa','ptt', 'err_ptt', 'pab','err_pab','pbb','err_pbb']
        pk_everything = np.column_stack((pkdata.k,pkdata.z,
                        pkdata.Paa, err_paa,
                        pkdata.Ptot, err_ptt,
                        pkdata.Pab,err_pab,
                        pkdata.Pbb, err_pbb,))
        df_pk = pd.DataFrame(pk_everything,columns=columns)
        df_pk.to_csv(inis.save_pk_with_err_path, index=False)
        # np.savetxt(inis.save_pk_with_err_path,pk_everything)

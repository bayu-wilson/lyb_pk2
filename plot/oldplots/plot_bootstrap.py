#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import inis
import options as opt

##### Control Area ######
# PLOTTING ratio of P?/P_nocorr
show_plot = False
# ['k','Paa','Ptot','Pab','Qab','Pbb','num','z']
# ['mf_a', 'nvar_a','dloglambda_a', 'npow_a',
#           'mf_tot', 'nvar_tot','dloglambda_tot', 'npow_tot','z','mf_b',
#           'var_atot','npow_atot']
#########################

y = np.loadtxt(inis.save_boot_mf_path)
p = np.loadtxt(inis.save_boot_pk_path)
mf_bootstrap = np.reshape(y,(12,opt.zbinlen,inis.M))
pk_bootstrap = np.reshape(p,(8,opt.kbinlen*opt.zbinlen,inis.M))
mf_a = mf_bootstrap[0]
mf_t = mf_bootstrap[4]
mf_b = mf_bootstrap[9]
pk_a = pk_bootstrap[1]
pk_t = pk_bootstrap[2]
pk_ab = pk_bootstrap[3]
mf_list = [mf_a,mf_t,mf_b]
pk_list = [pk_a,pk_t,pk_ab]
# print(np.nansum(pk_list[-1].T[5]))
mf_labels  = [r"$\overline{F}_{\alpha}$",
              r"$\overline{F}_{tot}$",
              r"$\overline{F}_{\beta}$"]
pk_labels  = [r"$P_{\alpha\alpha}$",r"$P_{TT}$",r"$P_{T\alpha}$"]
# Add a colorbar
fig,ax = plt.subplots(6, 2,gridspec_kw={'width_ratios': [1,1]})
fig.set_size_inches(12,8*4)

### PLOTTING BOOTSTRAP CONVERGENCE MF AA
# x = [np.var(mf_a.T[:i].T) for i in range(1,inis.M)]
# ax[0,1].plot(x)
# med = np.median(x)
# ax[0,1].set_ylim(med - med*0.01,med + med*0.01)
# ax[0,1].set_ylabel(r"Variance of $\overline{F}$ measurement")

# ### PLOTTING BOOTSTRAP CONVERGENCE FOR PK
# x = [np.var(pk_a.T[:i].T) for i in range(1,inis.M)]
# ax[1,1].plot(x)
# med = np.median(x)
# ax[1,1].set_ylim(med - med*0.01,med + med*0.01)
# ax[1,1].set_xlabel("Bootstrap Iterations")
# ax[1,1].set_ylabel("Variance of P(k) measurement")

# ### PLOTTING BOOTSTRAP CONVERGENCE MF TT
# x = [np.nanvar(mf_t.T[:i].T) for i in range(1,inis.M)]
# ax[2,1].plot(x)
# med = np.median(x)
# ax[2,1].set_ylim(med - med*0.01,med + med*0.01)
# ax[2,1].set_xlabel("Bootstrap Iterations")
# ax[2,1].set_ylabel("Variance of P(k) measurement")

#### PLOTTING MEAN FLUX COV MATRICES
for j,axi in enumerate(ax[:len(mf_list),0]):
    im = axi.imshow(np.corrcoef(mf_list[j]), cmap = plt.cm.jet)
    fig.colorbar(im, ax=axi)
    im.set_clim(0.0, 1.0)

    axi.set_xticks(np.arange(0, opt.zbinlen, step=1))
    axi.set_xticklabels(opt.zbin_centers)
    axi.set_xticks(np.arange(1/2, opt.zbinlen, step=1),'minor')

    axi.set_yticks(np.arange(0, opt.zbinlen, step=1))
    axi.set_yticklabels(opt.zbin_centers)
    axi.set_yticks(np.arange(1/2, opt.zbinlen, step=1),'minor')
    axi.set_ylabel("r({0})".format(mf_labels[j]))

#### PLOTTING PK COV MATRICES
for j,axi in enumerate(ax[len(mf_list):,0]):
    idx = j+1
    im = axi.imshow(np.corrcoef(pk_bootstrap[idx]), cmap = plt.cm.jet)
    fig.colorbar(im, ax=axi)
    im.set_clim(0.0, 1.0)

    axi.set_xticks(np.arange(opt.kbinlen/2, opt.kbinlen*opt.zbinlen, step=opt.kbinlen))
    axi.set_xticklabels(opt.zbin_centers)
    axi.set_xticks(np.arange(0, opt.kbinlen*opt.zbinlen, step=opt.kbinlen),'minor')

    axi.set_yticks(np.arange(opt.kbinlen/2, opt.kbinlen*opt.zbinlen, step=opt.kbinlen))
    axi.set_yticklabels(opt.zbin_centers)
    axi.set_yticks(np.arange(0, opt.kbinlen*opt.zbinlen, step=opt.kbinlen),'minor')
    axi.set_ylabel("r({0})".format(pk_labels[j]))

### PLOTTING BOOTSTRAP CONVERGENCE MF
for j,axi in enumerate(ax[:len(mf_list),1]):
    x = [np.nanvar(mf_list[j].T[:i].T) for i in range(1,inis.M)]
    axi.plot(x)
    med = np.median(x)
    axi.set_ylim(med - med*0.01,med + med*0.01)
    axi.set_ylabel(r"Variance of {0}".format(mf_labels[j]))
#
### PLOTTING BOOTSTRAP CONVERGENCE FOR PK
for j,axi in enumerate(ax[len(mf_list):,1]):
    x = [np.nanvar(pk_list[j].T[:i].T) for i in range(1,inis.M)]
    axi.plot(x)
    med = np.median(x)
    axi.set_ylim(med - med*0.01,med + med*0.01)
    axi.set_xlabel("Bootstrap Iterations")
    axi.set_ylabel("Variance of {0}".format(pk_labels[j]))



fig.tight_layout()
if inis.save_boot_fig:
    fig.savefig(inis.save_boot_fig_path)
if show_plot:
    plt.show()
plt.clf()

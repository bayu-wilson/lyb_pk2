#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import options as opt
import inis

##### Control Area ######
plot_inis = True
show_plot = True
# fontsize = 15
#########################

"""
F(z)
['mf_a', 'nvar_a','dloglambda_a', 'npow_a',
         'mf_tot', 'nvar_tot','dloglambda_tot', 'npow_tot','z','mf_b',
         'var_atot','npow_atot']

P(k,z)
['k','Paa','Ptot','Pab','Qab','Pbb','num','z']
"""
#mask
cover_mult = np.zeros((3*opt.zbinlen*opt.kbinlen,3*opt.zbinlen*opt.kbinlen),dtype=bool)
mask = np.ones((3*opt.zbinlen*opt.kbinlen,3*opt.zbinlen*opt.kbinlen),dtype=bool)
for i in range(int(273/opt.kbinlen)):
    for j in range(int(273/opt.kbinlen)):
        #getting rid of asymmetry, when z1!=z2
        if (i%int(91/opt.kbinlen)) == (j%int(91/opt.kbinlen)):
            #print(i,j)
            cover_mult[i*opt.kbinlen:(i+1)*opt.kbinlen,
                  j*opt.kbinlen:(j+1)*opt.kbinlen]= True

        #
        if (i%int(91/opt.kbinlen)<2)|((j%int(91/opt.kbinlen)<2)): #91/13= 7
            mask[i*opt.kbinlen:(i+1)*opt.kbinlen,
                  j*opt.kbinlen:(j+1)*opt.kbinlen]= False




labels = [r'$P_{\alpha \alpha}$',r'$P_{TT}$',r'$P_{\alpha \beta}$']

plt.rc('xtick', labelsize=15)    # fontsize of the tick labels
plt.rc('ytick', labelsize=15)

fig,ax = plt.subplots(1)
fig.set_size_inches(7,7)


if plot_inis:
    m = np.loadtxt(inis.save_boot_mf_path)
    p = np.loadtxt(inis.save_boot_pk_path)
    N = opt.kbinlen*(opt.zbinlen-2)
    #Nx = opt.kbinlen*(opt.zbinlen-2)
    tick_labels = np.concatenate([opt.zbin_centers[2:]]*3)
    ax.tick_params(axis='both',which='major',bottom='off',left='off',pad=15)
    ax.tick_params(axis='both',which='minor',bottom='off',left='off')
    # plt.tick_params(axis='both', which='both', bottom='off', top='off',
    #labelbottom='off', right='off', left='off', labelleft='off')

    # mask = np.ones(3*opt.kbinlen*opt.zbinlen,dtype=bool)# ones
    # mask[[i+(0+opt.zbinlen*0)*opt.kbinlen for i in range(0,opt.kbinlen)]]=False
    # mask[[i+(1+opt.zbinlen*0)*opt.kbinlen for i in range(0,opt.kbinlen)]]=False
    # mask[[i+(0+opt.zbinlen*1)*opt.kbinlen for i in range(0,opt.kbinlen)]]=False
    # mask[[i+(1+opt.zbinlen*1)*opt.kbinlen for i in range(0,opt.kbinlen)]]=False
    # mask[[i+(0+opt.zbinlen*2)*opt.kbinlen for i in range(0,opt.kbinlen)]]=False
    # mask[[i+(1+opt.zbinlen*2)*opt.kbinlen for i in range(0,opt.kbinlen)]]=False
    # r = np.corrcoef(p[N:-N*4])
    # r = r[mask,:][:,mask]
    #r_mf = np.corrcoef(m)
    r_pk = np.corrcoef(p)*cover_mult #np.cov(p)#np.corrcoef(p)
    new_N = int(np.sqrt(np.sum(mask)))
    r_pk = np.reshape(r_pk[mask],(new_N,new_N))
    im = ax.imshow(r_pk, cmap = plt.cm.jet)
    fig.colorbar(im)
    im.set_clim(-1.0, 1.0) #july31
    # LINES
    for i in range(1,15):
        ax.axvline(i*opt.kbinlen,ls='dotted',color='gray',lw=0.75)
        ax.axhline(i*opt.kbinlen,ls='dotted',color='gray',lw=0.75)
    ax.axvline(N, color='white', ls='-', lw=1)
    ax.axvline(2*N, color='white', ls='-', lw=1)
    ax.axhline(N, color='white', ls='-', lw=1)
    ax.axhline(2*N, color='white', ls='-', lw=1)

    #ax.set_xticks(np.arange(0, Nx*4, step=Nx),'minor')
    ax.set_xticks(np.arange(N/2, N*3, step=N)) #major
    ax.set_xticklabels(labels)
    ax.set_xticks(np.arange(opt.kbinlen/2, N*3, step=opt.kbinlen),'minor')
    ax.set_xticklabels(tick_labels,minor=True,fontsize=5)

    ax.set_yticks(np.arange(N/2, N*3, step=N)) #major
    ax.set_yticklabels(labels)
    ax.set_yticks(np.arange(opt.kbinlen/2, N*3, step=opt.kbinlen),'minor')
    ax.set_yticklabels(tick_labels,minor=True,fontsize=5)


    #ax.grid(which="minor", ls='dotted',color='gray',lw=0.75)

# np.savetxt("/Users/bayuwilson/Desktop/cov_matrix.txt",r_pk)

if inis.save_paper_covmatrix:
    plt.savefig(inis.save_paper_covmatrix_path)
print("Saving Here:\n",inis.save_paper_covmatrix_path)
if show_plot:
    plt.show()
plt.clf()

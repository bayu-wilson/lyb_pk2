#!/usr/bin/env python
import numpy as np
import options as opt
import inis
import matplotlib.pyplot as plt

show_plot = False
fontsize=25


def total_path(z):
    #dz = 3e-5*np.log(10)
    dv = opt.c_kms*(3.e-5)*np.log(10)
    return dv/opt.H_0/np.sqrt(opt.omega_m0*(1+z))
    #(2*opt.c_kms/(70*np.sqrt(0.3)) *((1+z-dz)**(-0.5)-(1+z+dz)**(-0.5)))

print(inis.save_kzq_pk_path)
qsos_table = np.loadtxt(inis.save_kzq_pk_path).T
#qso index, power type, z,k, sum of power in bin, num of pix in bin
#power type: paa=0, ptt=1, pab=2

lya_pix_dist = np.zeros(opt.zbinlen)
lyb_pix_dist = np.zeros(opt.zbinlen)
lyx_pix_dist = np.zeros(opt.zbinlen)
for zidx in range(opt.zbinlen):
    zmask_a = (qsos_table[2] == opt.zbin_centers[zidx])&(qsos_table[1] == 0)
    zmask_b = (qsos_table[2] == opt.zbin_centers[zidx])&(qsos_table[1] == 1)
    zmask_x = (qsos_table[2] == opt.zbin_centers[zidx])&(qsos_table[1] == 2)

    lya_pix_dist[zidx] = np.sum(qsos_table[5][zmask_a])
    lyb_pix_dist[zidx] = np.sum(qsos_table[5][zmask_b])
    lyx_pix_dist[zidx] = np.sum(qsos_table[5][zmask_x])

fig,ax = plt.subplots(1)
fig.set_size_inches(9,6)

path_factor = total_path(opt.zbin_centers)
ax.bar(opt.zbin_centers, path_factor*lya_pix_dist,width=0.15,color='red', label=r"Ly$\alpha$ forest")
ax.bar(opt.zbin_centers, path_factor*lyb_pix_dist,width=0.12,color='blue', label=r"Ly$\beta$ forest")
# ax.bar(opt.zbin_centers, lyx_pix_dist,width=0.1,color='green', label=r"Ly$\alpha$-Ly$\beta$ overlap")
ax.set_xlabel("z", fontsize=fontsize)
ax.set_ylabel("Mpc", fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=20)

ax.legend(loc='upper right', fontsize=fontsize)
plt.tight_layout()

if inis.save_paper_pixel_dist:
    plt.savefig(inis.save_paper_pixel_dist_path)
print("Saving Here:\n",inis.save_paper_pixel_dist_path)
if show_plot:
    plt.show()
plt.clf()

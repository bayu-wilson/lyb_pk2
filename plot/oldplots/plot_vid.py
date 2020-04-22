#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

## get T_COLD and initial power spectrum of the simulation

path = "../data/sims/symp_model/"
testme = "pkT_LCDM_COLDg1.0_f1.0_z3.800.txt"
pi = np.pi
with open(path+testme) as header:
    lines = header.readlines()
    T0 = float(lines[1].rstrip().split(' ')[-2])
    gamma = float(lines[1].rstrip().split(' ')[-1])
pk_cold_g10_z38 = np.loadtxt(path+testme,skiprows=4)
# print(pk_cold_g10_z38)
kp = pk_cold_g10_z38[:,0]#pk_ref_z38[:,0]
## print what T_COLD and gamma in the initial simulation are
print(T0,gamma)

## defines the window function
def wink(k,T1,g,line='lya'):
    if (line=='lya'):
        dT = T1*3**(g-1.)
    elif (line=='lyb'):
        dT = T1*12**(g-1.)
    sigma = 14.0 * np.sqrt(dT/1e4) # in km/s
    print(np.sqrt(T0*T0 + dT*dT), np.sqrt(T0*T0 + T1*T1))
    return np.exp(-k*k*sigma*sigma)

## alpha-alpha
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1],color='red',linestyle='-.',alpha=1)
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1]*wink(kp,5000,1.0,'lya'),color='red',linestyle=':',alpha=1)
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1]*wink(kp,5000,1.6,'lya'),color='red',linestyle='-',alpha=1)
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,1]*wink(kp,10000,1.6,'lya'),color='red',linestyle='--',alpha=1)

## T-alpha
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3],color='blue',linestyle='-.',alpha=1)
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3]*wink(kp,5000,1.0,'lyb'),color='blue',linestyle=':',alpha=1)
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3]*wink(kp,5000,1.6,'lyb'),color='blue',linestyle='-',alpha=1)
plt.plot(np.log10(kp),kp/pi*pk_cold_g10_z38[:,3]*wink(kp,10000,1.6,'lyb'),color='blue',linestyle='--',alpha=1)

plt.yscale('log')
plt.xlim(-2.4,-1.2)
plt.ylim(1e-2,3e-1)
# plt.xscale('log')

plt.show()

#!/usr/bin/env python

# VID'S RELATIVE ERROR PLOTTING ROUTINE. (X,Y) AXES: (ERR/P,K)
# OUTPUT: [MANUAL]

from scipy import *
from scipy import interpolate
from scipy import signal
import numpy
import pylab as plt

def get_figsize(fig_width_cm):
    inches_per_cm = 0.393700787             # Convert cm to inch
    golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_width = fig_width_cm*inches_per_cm  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

# params2 = {'backend': 'pdf',
#           'axes.labelsize': 24,
#           'text.fontsize': 24,
#           'xtick.labelsize': 24,
#           'ytick.labelsize': 24,
# 	      #'legend.draw_frame': False,
#           'legend.fontsize': 18,
#           'lines.markersize': 6,
#           'font.size': 28,
# #          'font.family':'cursive',
#            'font.family':'serif',
#            'font.family':'Utopia',
#           'mathtext.fontset':'cm',
#           'text.usetex': True,
#           'figure.figsize': get_figsize(32)}
# plt.rcParams.update(params2)
########################################################
## Bayu's data
import pandas as pd
path='../output/'
data = pd.read_csv(path+'pk_errboot_obs_corrNR.txt').values
#loadtxt(path+'pk_errboot_obs_corrNR.txt',skiprows=1,delimiter=',')
zlist_by = unique(data[:,1])
klist_by = unique(data[:,0])
nz_by = len(zlist_by)
nk_by = len(klist_by)
bayu = []
for z in zlist_by:
    bayu.append(data[where(data[:,1]==z)[0]])

path='../data/obs/vid_2020/'
# data = pd.read_csv(path+'pk_xs_final.txt').values
data = loadtxt(path+'pk_xs_final.txt')
zlist_xq = unique(data[:,0])
klist_xq = unique(data[:,1])
nz_xq = len(zlist_xq)
nk_xq = len(klist_xq)
xq = []
for z in zlist_xq:
    xq.append(data[where(data==z)[0]])


########################################################
plt.clf()

bs_xq = sqrt(1.3)

c=['red','blue','green']
c2=['green','magenta','orangered']
c3=['red','blue','green','purple','orangered']

## xq100
kstart = 0
kend = len(xq[0][:,1])
j=0
for i in range(nz_xq):
    if (zlist_xq[i] == 3.4 or zlist_xq[i] == 3.6 or zlist_xq[i] == 3.8 or zlist_xq[i] == 4.0 or zlist_xq[i] == 4.2):
        myx = log10(xq[i][kstart:kend,1])
        yy = xq[i][kstart:kend,2]
        ee = xq[i][kstart:kend,3]
        plt.plot(myx,ee/yy,color=c3[j],linestyle='--')
        j=j+1

## bayu
kstart = 0
kend = len(bayu[0][:,0])
j=0
for i in range(nz_by):
    if (zlist_by[i] == 3.4 or zlist_by[i] == 3.6 or zlist_by[i] == 3.8 or zlist_by[i] == 4.0 or zlist_by[i] == 4.2):
        myx = log10(bayu[i][kstart:kend,0])
        yy = bayu[i][kstart:kend,2]
        ee = bayu[i][kstart:kend,3]
        plt.plot(myx,ee/yy,color=c3[j],linestyle='-')
        j=j+1

plt.plot(0,0,color='black',linestyle='--',label='Irsic+17')
plt.plot(0,0,color='black',linestyle='-',label='Wilson+19')

#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')

#plt.gca().set_xlim(2.9e-3,0.1)
plt.gca().set_xlim(-2.4,-1.1)
#plt.gca().set_ylim(8e-2,3e-1)
plt.gca().set_ylim(0.03,0.20)

plt.xlabel('$k\\;[\\mathrm{km^{-1}\,s}]$')
plt.ylabel('$\\sigma_{P}/P_F$')
plt.legend(loc=0,numpoints=1,ncol=2).draw_frame(False)
plt.show()
# plt.savefig('slike/pk_bayu_error.pdf')

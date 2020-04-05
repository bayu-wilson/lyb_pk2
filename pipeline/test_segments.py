import numpy as np
###########################
#effect of nearest grid point interpolation to add exposures on power spectrum and its cross

import numpy as np
from scipy import interpolate, integrate, signal





##################################################################################
#velocity range
###################################################################################
dv = 2 #km/s way less than pixel to generate large array and include aliasing effects from scales down to thermal
vrange = 1e7*dv
N = int(vrange/dv)

####################################################################################
#K bin range
###################################################################################
kmin = 10**-2.5 #3e-3
kmax = 10**-.5 #6e-2
kbinlen = 20
kbin_centers = np.logspace(np.log10(kmin),np.log10(kmax),kbinlen)
dk = (np.log10(kmax) - np.log10(kmin))/(kbinlen-1) #feb4
kbin_edges = np.logspace(np.log10(kmin)-dk/2,np.log10(kmax)+dk/2,kbinlen+1)

#parameters for generating power spectra
vDoppler = 14  #km/s
tau0 = .5 #optical depth at mean density of universe
random_seed = 12 #for initializing random number generator

def onedimDensityPower(kvec):
    return 20*kvec**-0.2;  #adjusted amplitude to get reasonable mean flux 

def ThermalWindowSq(kvec):
    return np.exp(-vDoppler**2*kvec**2)

def get_autopower(transmission):
    """
    Computes auto-power spectrum.
    """
    rf_x = transmission - 1 #relative flux fluctuations
    rf_k = np.fft.fft(rf_x)*dv # relative fft of flux fluctuations
    N = len(rf_k)
    V = N*dv #one dimensional volume
    dk = 2*np.pi/V
    k = dk*np.arange(0,N,1)
    pk = np.abs(rf_k)**2/V
    return k,pk

def get_kmask(kpix,kidx,kedges): # line 190 main.py
    """
    Mask to get k-bins.
    """
    mask = ((kpix>=kedges[kidx])&(kpix<kedges[kidx+1]))#&(kpix!=0))
    return mask

  
########################################################
#Make Gaussian random field
########################################################
print("Making one large Gaussian Random Field...")
numk =N//2+1
vvec = np.array([dv*i for i in range(0, N)])
kvec = np.array([2*np.pi/vrange*i for i in range(0, numk)])
ivec = np.arange(N) 
np.random.seed(random_seed)
densityfieldk = np.random.normal(size=numk) + np.random.normal(size=numk)*1j
densityfieldk[-1] = np.real(densityfieldk[-1]) #last component should not be imaginary
densityfieldk[1:]= np.sqrt(vrange*onedimDensityPower(kvec[1:])*ThermalWindowSq(kvec[1:])/2.)*densityfieldk[1:]

print("FFTing to reals space..."),
densityfieldx= 1+ np.fft.irfft(densityfieldk)*N/vrange
print("done.")
print(densityfieldx)

var = np.var(densityfieldx)
transmission = np.exp(-tau0*np.exp(densityfieldx-var/2)**2)  #lognormal mapping where var/2 is to make mean of density field equal to one 

meanFlux = np.mean(transmission)
print("meanFlux = ", meanFlux)


#true power spectrum
pkfullarr = np.zeros(kbinlen);
countfullarr = np.zeros(kbinlen);
kpix, Pk = get_autopower(transmission) #entire array
for kidx in range(kbinlen):
    mask = get_kmask(kpix,kidx,kbin_edges)
    pk_sub = Pk[mask]
    pkfullarr[kidx] += np.sum(pk_sub) #adding ptot in k,z bin for the qidx'th quas
    countfullarr[kidx] += len(pk_sub) #adding number of ptot pixels

pkfullarr/=countfullarr #devide out number of samples

###################################################################
# Fourier transform segments and entire
###############################################################
pksegarr = np.zeros(kbinlen);
countsegarr = np.zeros(kbinlen);

numsegments = 1000
velocitybetweensegments = 2000
numpixelstep = int(velocitybetweensegments/dv)
ipix=0
for i in range(numsegments):
    pixels_in_segment = 1000 #equal to 100 20km/s  

    kpix, Pk = get_autopower(transmission[ipix:ipix+pixels_in_segment])
    for kidx in range(kbinlen):
        mask = get_kmask(kpix,kidx,kbin_edges)
        pk_sub = Pk[mask]
        pksegarr[kidx] += np.sum(pk_sub) #adding ptot in k,z bin for the qidx'th quas
        countsegarr[kidx] += len(pk_sub) #adding number of ptot pixels
        
    ipix += pixels_in_segment+numpixelstep     

pksegarr/=countsegarr #devide out number of samples

    

   # [np.mean(transmission[numsamplingpixel*i:numsamplingpixel*(i+1)]) for i in xrange(Npix)]
#######################################################################
#make real space plot of smoothed fields
#######################################################################
import matplotlib.pylab as plt
params = {'xtick.labelsize': 13, 'ytick.labelsize': 13, 'font.family': 'serif','axes.formatter.limits': (-4,4), 'xtick.major.size':8, 'ytick.major.size':8, 'xtick.minor.size':4, 'ytick.minor.size': 4}
plt.rcParams.update(params)
fig = plt.figure(7,figsize=(7.0, 9.))
ax1 = fig.add_subplot(111)

#title = identification_string #"total numk = ", N/2+1, "reconst. numk = ", numk
#ax1.set_title(title, fontsize=12)
ax1.set_xlabel("k [s/km]")
ax1.set_ylabel(r"kP(k)/\pi")

#ax1.loglog(kbin_centers, kbin_centers*pkfullarr/np.pi, label="full")
#ax1.loglog(kbin_centers, kbin_centers*pksegarr/np.pi, label="segments")
ax1.semilogx(kbin_centers, pksegarr/pkfullarr, label="full")
plt.xlim(2*np.pi/2000, 1)
#plt.ylim(.001,.1)
#ax1.legend(loc="lower left", borderaxespad=0, fontsize=16, ncol=3)
figname = 'test.eps'
fig.savefig(figname)
print("output ", figname)


###############################################################
#Let's calclate histograms of bin lengths
###############################################################
fig.clf()


#npzfile = np.load('Lyaskewerlengths.npz')
npzfile = np.load('Lytskewerlengths.npz')
print(npzfile.files)

plt.rcParams.update(params)
fig = plt.figure(7,figsize=(7.0, 9.))
ax1 = fig.add_subplot(111)

#title = identification_string #"total numk = ", N/2+1, "reconst. numk = ", numk
#ax1.set_title(title, fontsize=12)
ax1.set_xlabel("num pixels segment")
ax1.set_ylabel(r"N")

ax1.hist(npzfile['arr_0'][1:], edgecolor='black', ls='dashed',lw=3, facecolor="None", label='z=3.0')
ax1.hist(npzfile['arr_1'][1:],edgecolor='blue', ls='dashed', lw=3, facecolor="None", label='z=3.2')
ax1.hist(npzfile['arr_2'][1:], edgecolor='green', ls='dashed',lw=3, facecolor="None", label='z=3.4')
ax1.hist(npzfile['arr_3'][1:], edgecolor='cyan', ls='dashed',lw=3, facecolor="None", label='z=3.6')
ax1.hist(npzfile['arr_4'][1:], edgecolor='red', ls='dashed',lw=3, facecolor="None", label='z=3.8')
ax1.hist(npzfile['arr_5'][1:], edgecolor='magenta', ls='dashed',lw=3, facecolor="None", label='z=4.0')
ax1.hist(npzfile['arr_6'][1:], edgecolor='orange', ls='dashed', lw=3, facecolor="None", label='z=4.2')
print(npzfile['arr_4'])
#plt.xlim(2*np.pi/2000, 1)
#plt.ylim(.001,.1)
ax1.legend(loc="upper right", borderaxespad=0, fontsize=12)
figname = 'testhist.eps'
fig.savefig(figname)
print("output ", figname)

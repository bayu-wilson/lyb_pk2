###########################
#effect of nearest grid point interpolation to add exposures on power spectrum and its cross

import numpy as np
from scipy import interpolate, integrate, signal

dvpixel = 7 #Pixel size
dvXS = 20.
numsamplingpixel = 7
dv = dvpixel/numsamplingpixel #way less than pixel to generate large array and include aliasing effects from scales down to thermal
vrange = 1e6*dv
N = int(vrange/dv)
Npix = int(vrange/dvpixel)
NXS = int(vrange/dvXS)
Numexposures = 12
vDoppler = 14  #km/s
tau0 = .5 #optical depth at mean density of universe
random_seed = 12 #for initializing random number generator

def onedimDensityPower(kvec):
    return 5*kvec**-0.5;  #adjusted amplitude to get reasonable mean flux 

def ThermalWindowSq(kvec):
    return np.exp(-vDoppler**2*kvec**2)
  
########################################################
#Make Gaussian random field
########################################################
print("Making one large Gaussian Random Field...")
numk =N//2+1
vvec = np.array([dv*i for i in range(0, N)])
vpixvec = np.array([dvpixel*i for i in range(0, Npix)])
vXSvec = np.array([dvXS*i for i in range(0, NXS)])
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

#############################################################
#gridding
#############################################################
transdvpixel_master = [np.mean(transmission[numsamplingpixel*i:numsamplingpixel*(i+1)]) for i in range(Npix)]
onearr = np.ones(Npix)

NumSkewers = 12
vbins = [np.array]*NumSkewers
binnedSpectrum = np.zeros(len(vXSvec )); numberinBin = np.zeros(len(vXSvec ))


for i in range(NumSkewers):
    nshift = int(dvXS*(np.random.rand()-0.5)/dv)
    print("nshift = ", nshift)
    vbins[i] = vpixvec + nshift*dv
    #putting shifted spectra in bins
    transdvpixel =  [np.mean(transmission[(numsamplingpixel*j+nshift)%N:(numsamplingpixel*(j+1)+nshift)%N]) for j in range(Npix)]
    
    inds = np.digitize(vbins[i], vXSvec)

    #print(len(binnedSpectra[inds-1]), len(transdvpixel))

    binnedSpectrum[inds-1] += transdvpixel  #Not doing what I want it to
    numberinBin[inds-1] += 1

binnedSpectrum/=numberinBin
    

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
ax1.set_xlabel("position [km/s]")
ax1.set_ylabel(r"\delta")

ax1.plot(vpixvec + 0.5*dvpixel, transdvpixel_master, label="input")
ax1.plot(vXSvec + 0.5*dvXS, binnedSpectrum, label="input")

plt.xlim(100, 300);
#plt.ylim(-3*np.sqrt(sigmasq),3*np.sqrt(sigmasq))
#ax1.legend(loc="lower left", borderaxespad=0, fontsize=16, ncol=3)
figname = 'test_spectra.eps'
fig.savefig(figname)
print("output ", figname)

#delta = []*Numexposures

#generates Nexposures
#for(i = 0; i< Nexposure; i++)
#    phase = dv*np.random.rand()
#    delta[i] = *npexp(ij*phase*kvec)


#claclate inverse FFT


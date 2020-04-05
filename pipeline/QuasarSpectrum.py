import pandas as pd
import matplotlib.pyplot as plt
import options as opt
import numpy as np
import sys
from scipy.interpolate import griddata
import glob
import inis

class QuasarSpectrum(object):
    all_names = None
    all_redshifts = None
    dla_table = None
    seeing_table = None
    nqso = None
    """Container for quasar spectra"""
    def __init__(self,wavelength=None, flux=None, err_flux=None,
                 resolution=None, mask_dla=None, dloglambda=None, kmax=None, name=None, redshift=None,
                 unnormalized_flux=None):
        """
        Parameters
        ----------
        wavelength : array-like
            Wavelength of each flux measurement
        flux : array-like
            Transmitted fluxes at each wavelength (normalized)
        err_flux : array-like
            Uncertainties of each flux measurement
        resolution : array-like
            Spectroscopic resolution (in Angstroms)
        dloglambda : array-like
            Difference of wavelength pixels in log space.
        name : string
            Name of the target
        redshift : float
            Redshift of target
        unnormalized_flux : array-like
            Radiative flux of observed quasar.
        """
        self.wavelength = wavelength
        self.flux = flux
        self.err_flux = err_flux
        self.name = name
        self.resolution = resolution
        self.dloglambda = dloglambda
        self.redshift = redshift
        self.unnormalized_flux = unnormalized_flux
        self.mask_dla = mask_dla  #will be set if inis.remove_dla == True
        self.kmax = kmax #set in resolution if there is a maximum wavenumber 
    def plot_spectrum(self,norm=True):
        """
        Plots normalized or unnormalized qso spectrum.

        Parameters
        ----------
        norm : bool
            If `True`, then normalized (transmitted) flux is used. If `False`, then
            unnormalized (radiative) flux is used. Currently onle applicable to observed data.
        """
        fig,ax = plt.subplots(1)
        fig.set_size_inches(15,4)
        fontsize = 15
        ax.grid()
        ax.set_xlabel("Wavelength [Å]",fontsize=fontsize)

        if (not norm) & (not self.name.startswith("mock")): # if normalized and observed data
            ax.plot(self.wavelength, self.unnormalized_flux/1e-15)
            ax.set_ylabel(r"Flux [$10^{-15}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$]",fontsize=fontsize)
        else:
            ax.plot(self.wavelength, self.flux)
            ax.set_ylabel(r"Flux [erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$]",fontsize=fontsize)
        plt.show()
    @classmethod
    def load_cat(cls,tag):
        """
        Load a catalog of sightline names and source redshifts from a raw text file.
        Possible inputs are obs, mocks, mocks_n5000

        Parameters
        ----------
        name : str
            Name of the quasar in text file.
        """
        if tag.startswith('obs'): #What matters is what the tag starts/ends with. This is loading a catalog for OBS
            path_to_cat = "../data/obs/XQ-100_catalogue.txt"
            catalog_table = pd.read_csv(path_to_cat,delim_whitespace=True,usecols=(3,6),
                                names = ("qso_name","redshifts"), comment='#')
            path_to_dlas = "../data/obs/XQ-100_DLA_catalogue.txt"
            dla_colnames = ['idx','name','z','NHI','err_NHI']
            ####
            dla_table = pd.read_csv(path_to_dlas,delim_whitespace=True,names=dla_colnames)
            path_to_seeing = "../data/obs/XQ-100_catalogue_seeing.txt"
            seeing_table  = pd.read_csv(path_to_seeing,delim_whitespace=True,names=['qsonum', 'catelog', 'longname', 'name', 'RA', 'DEC', 'z_qso', 'aperture', 'seeing_min', 'seeing_max'])
            #print("seeing min= ", seeing_table['seeing_min'][seeing_table['name']=='J1024+1819'].values)  #Matt M: testing that this works
            #exit();
            
            qname_array,z_array = catalog_table['qso_name'].values,catalog_table['redshifts'].values
            nqso = 100
            cls.dla_table = dla_table
            cls.seeing_table = seeing_table
        elif tag.startswith('mocks'): #This is loading a catalog for MOCKS
            cat_mask = int(tag.split("_n")[1]) == np.array([100,600,700,5000])
            #which_catalog = int(tag.endswith("n5000")) # index for `path_to_cat`
            path_to_cat = np.array(["../data/mocks/XQ-100_catalogue_n100.mock",
                           "../data/mocks/XQ-100_catalogue_n600.mock",
                           "../data/mocks/XQ-100_catalogue_n700.mock",
                           "../data/mocks/XQ-100_catalogue_n5000.mock"])[cat_mask][0]

            catalog_table = pd.read_csv(path_to_cat,delim_whitespace=True,usecols=(3,6),
                                names = ("qso_name","redshifts"),float_precision='high')
            nqso = int(path_to_cat.split("_n")[1][:-5])
        else:
            print("\n Error! Please input a valid tag. \n")
        cls.all_names = catalog_table.qso_name.values
        cls.all_redshifts = catalog_table.redshifts.values
        cls.nqso = nqso
    @classmethod
    def load_qso_data(cls,cat_index, tag, rescale_flux):
        """
        Load a quasar spectrum from a raw text file.

        Parameters
        ----------
        name : str
            Name of the quasar in text file. i.e. "J0003-2603" or "mock-0000"
        """
        name = cls.all_names[cat_index]
        redshift = cls.all_redshifts[cat_index]

        kmax_dict = {'UV': 100., 'VIS': 100.} #initialize to large number
        if name.startswith("mock"):
            ### MOCKS ###
            all_paths = glob.glob("../data/mocks/XQ-100_{0}/ly*/ref/spectra/mock-*.txt.gz".format(
                                    tag))
            all_paths.sort()
            path = all_paths[cat_index]
            qso_table = pd.read_csv(path, delim_whitespace=True,compression='gzip',skiprows = 1,usecols=(0,1,2,3,4),
                                    names = ("wav", "flx", "ferr", "res","dloglam"),float_precision='high')
            mask_wave = ((qso_table.wav>opt.lyb_min*(1+redshift))&
                        (qso_table.wav<opt.lya_rest/opt.lyb_rest*opt.lyb_max*(1+redshift))&
                        (qso_table.flx>opt.min_trans))
            qso_table=qso_table[mask_wave]
            wavelength = qso_table.wav.values
            flux = qso_table.flx.values
            err_flux = qso_table.ferr.values
            unnormalized_flux=np.ones_like(qso_table.flx.values)*np.nan
            resolution = qso_table.res.values # june 24
            mask_dla =  np.ones_like(wavelength,'bool')  #no DLAs curently in mocks so this doens't do anything
        else:

            ### LOADING IN DATA NORMALLY ###
            released_path = "../data/obs/XQ-100/released/{0}_uvb-vis.txt".format(name)
            continuum_path = "../data/obs/XQ-100/continuum/{0}_cont.txt".format(name)
            qso_table = pd.read_csv(released_path, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
                                        names = ("wav", "flx", "ferr", "res","dloglam"))
            cont_table = pd.read_csv(continuum_path, delim_whitespace=True, skiprows = 1,usecols=(0,1),
            names = ("wav",'flx'))
            mask_wave = (qso_table.wav.values>opt.lyb_min*(1+redshift))&(
                qso_table.wav.values<opt.lya_rest/opt.lyb_rest*opt.lyb_max*(1+redshift))&(
                (qso_table.flx.values/cont_table.flx.values)>opt.min_trans)  #Matt M: Second criteria is because opt.lyb_max is the highest redshift pixel that is used (and I suspect we keep Lyman-alpha for it in case we want to regenerate forest)
            # ### TESTING FEBRUARY 25 #feb25 #THERE ARE 2 BAD PIXELS AT Z= 3.25 AND 1.98
            # bad_pix_mask = qso_table.flx.values/cont_table.flx.values<opt.min_trans
            # asdf = np.sum(bad_pix_mask)
            # if asdf != 0:
            #     z_bad = qso_table.wav.values[bad_pix_mask]/opt.lya_rest-1
            #     print(asdf,z_bad)
            # ###
            qso_table = qso_table[mask_wave]
            cont_table = cont_table[mask_wave]
            wavelength = qso_table.wav.values
            unnormalized_flux = qso_table.flx.values
            flux = unnormalized_flux/cont_table.flx.values
            err_flux = qso_table.ferr.values/cont_table.flx.values
            resolution = qso_table.res.values #this will be overwritten if we use Carswell+18 resolution
            mask_dla =  np.ones_like(wavelength,'bool')

            if 'ADC' in tag and name in opt.ADC_off_qsos:
                print("ADC off: quasar ", name, " is not being used")
                mask_dla *= False 
            
            m1 = wavelength>=opt.overlap_maxwav #Carswell+18 #this is exactly where the overlap stops and it becomes the VIS arm
            if inis.carswell_res: #feb27
                resolution[m1] = opt.R_VIS_carswell
                resolution[~m1] = opt.R_UV_carswell
            else:
                slitVIS = .9; slitUV = 1.;

                #print("name = ", name, cls.seeing_table[cls.seeing_table.name==name].name.values)
                seeing_min = cls.seeing_table[cls.seeing_table.name == name].seeing_min.values
                seeing_max = cls.seeing_table[cls.seeing_table.name == name].seeing_max.values
                #print("seeing = ", seeing_min, seeing_max)
                FWHMseeing = 0.5*(seeing_min+seeing_max)

                if ("goodseeing" in tag and FWHMseeing <0.75) or ("badseeing" in tag and FWHMseeing >0.75) :
                    mask_dla *= False  #make it so quasar is not used
                    if "goodseeing" in tag:
                        print("goodseeing", FWHMseeing)
                    else:
                         print("badseeing", FWHMseeing)

                    
                #don't use quasar
                if np.isnan(FWHMseeing):
                    mask_dla *= False  #make it so quasar is not used
                    FWHMseeing = 1 #since not using it anyway this is easy way to continue

        
                         
                #print("FWHM = ", FWHMseeing, cls.seeing_table[cls.seeing_table.name == name])
                for arm in ['VIS', 'UV']:
                    x =(FWHMseeing/slitVIS if arm == 'VIS' else FWHMseeing/slitUV)
                    x2 =(seeing_min/slitVIS if arm == 'VIS' else seeing_min/slitUV)
                    x3 =(seeing_max/slitVIS if arm == 'VIS' else seeing_max/slitUV)
                    
                    if x < 0.8:
                        res = 14.1375 + 8.88945*(-0.65 + x) - 13.9827*(-0.65 + x)**2 + 16.4146*(-0.65 + x)**3
                        res2 = 14.1375 + 8.88945*(-0.65 + x2) - 13.9827*(-0.65 + x2)**2 + 16.4146*(-0.65 + x2)**3
                        res3 = 14.1375 + 8.88945*(-0.65 + x3) - 13.9827*(-0.65 + x3)**2 + 16.4146*(-0.65 + x3)**3
                    else:
                        res = 16.0069 + 3.02506*(-1 + x) - 4.15522*(-1 + x)**2 +   4.3044*(-1 + x)**3
                        res2 = 16.0069 + 3.02506*(-1 + x2) - 4.15522*(-1 + x2)**2 +   4.3044*(-1 + x2)**3
                        res3 = 16.0069 + 3.02506*(-1 + x3) - 4.15522*(-1 + x3)**2 +   4.3044*(-1 + x3)**3
                    if arm == 'UV':
                        res*=5400/4700
                        res2*= 5400/4700
                        res3*= 5400/4700
                        resolution[~m1] = res
                    else:
                        res*=5400/8900
                        res2*= 5400/8900
                        res3*= 5400/8900
                        resolution[m1] = res  #using ESO resolution; this is how it has to scale assuming intrinsic dispersion scales with FWHM of tophat
                                                            #(where intrinsic is Gaussian that is convolved with tophat)
                    if not np.isnan(seeing_min+seeing_max): 
                        kmax_dict[arm] = np.sqrt(.2/(res3**2 - res2**2 + 1e-5))  #maximum wavenumber for 10% error
                    print(arm, " FWHM = ", seeing_min, seeing_max, res, res2, kmax_dict[arm], opt.R_UV_carswell,  kmax_dict[arm]);
                    
        if "wR2" in tag and not inis.wR2:
            # Using new column
            resolution = np.ones_like(qso_table.res.values)*11 * 0.2
        if ("noB" in tag)&(inis.add_beta): #july21
            flux = flux**rescale_flux #np.exp(rescale_flux*np.log(flux)) #flux**rescale_flux #july20
        dloglambda = qso_table.dloglam.values

        
        if inis.remove_dla:
            cls.getDLAmask(name, mask_dla, wavelength, flux)
            
        return cls(name=name,redshift=redshift,
                   wavelength=wavelength, flux=flux, err_flux=err_flux,
                   resolution=resolution, mask_dla=mask_dla, dloglambda=dloglambda, kmax=kmax_dict,
                   unnormalized_flux=unnormalized_flux)

    @classmethod
    def getDLAmask(cls, name, mask_dla, wavelength, flux):
        """Masking DLAs
        """
        if name in cls.dla_table.name.values:
            sub_table = cls.dla_table[cls.dla_table.name == name]
            num_dlas = len(sub_table)
            #mask_dla = np.ones_like(wavelength,'bool')
            for n in range(num_dlas):
                row = sub_table.iloc[n]
                z_dla,NHI = row.z, row.NHI

                lambda_dla = opt.lya_rest*(1+ z_dla)
                lambda_dlb = opt.lyb_rest*(1+ z_dla)
                deltalambda_dla = opt.find_EW(NHI,"alpha",z_dla)*opt.DLA_cut_factor/2  #/opt.lya_rest
                deltalambda_dlb = opt.find_EW(NHI,"beta",z_dla)*opt.DLA_cut_factor/2
                dla_profile = opt.get_dla_profile(wavelength/opt.lya_rest/(1+z_dla) - 1.,10**NHI)
                dlb_profile = opt.get_dlb_profile(wavelength/opt.lyb_rest/(1+z_dla) - 1.,10**NHI)

                #mask central most observed part
                mask_dla *= (wavelength<(lambda_dla-deltalambda_dla))|(wavelength>(lambda_dla+deltalambda_dla)) #lyman alpha
                mask_dla *= (wavelength<(lambda_dlb-deltalambda_dlb))|(wavelength>(lambda_dlb+deltalambda_dlb)) #lyman beta
                
                #devide out effect in flux
                flux[mask_dla] = flux[mask_dla]/dla_profile[mask_dla]
                flux[mask_dla] = flux[mask_dla]/dlb_profile[mask_dla] 
    
  

                #import sys
                #import matplotlib.pyplot as plt
                #print("DLA specs", name, z_dla, NHI)
                #plt.plot(wavelength,flux)
                #plt.plot(wavelength[mask_dla], flux[mask_dla]/dla_profile[mask_dla])
                #plt.xlim(1216*(1+z_dla-.1), 1216*(1+z_dla+.1))
                #plt.show()
                    # sys.exit()
                    # plt.plot(zpix[mask*mask_dla],self.flux[mask*mask_dla])
                    # plt.plot(zpix[mask*mask_dla],self.flux[mask*mask_dla]/dla_profile[mask*mask_dla],color='red')


                    # wave_dla = (1+z_dla)*opt.lya_rest
                    # wave_dlb = (1+z_dla)*opt.lyb_rest
                    # lo_wave_dla = wave_dla-opt.find_EW(NHI,'alpha',z_dla)*opt.DLA_cut_factor/2
                    # hi_wave_dla = wave_dla+opt.find_EW(NHI,'alpha',z_dla)*opt.DLA_cut_factor/2
                    # lo_wave_dlb = wave_dlb-opt.find_EW(NHI,'beta',z_dla)*opt.DLA_cut_factor/2
                    # hi_wave_dlb = wave_dlb+opt.find_EW(NHI,'beta',z_dla)*opt.DLA_cut_factor/2

                    ### Dividing out the DLA profile
                    # mask_dla*=(self.wavelength<lo_wave_dla)|(self.wavelength>hi_wave_dla)
                    # mask_dla*=(self.wavelength<lo_wave_dlb)|(self.wavelength>hi_wave_dlb)

                    # dla_profile = opt.get_dla_profile((1.+z_dla)/(1.+zpix)-1,10**NHI)
                    # plt.plot(zpix[mask*mask_dla],self.flux[mask*mask_dla])
                    # plt.plot(zpix[mask*mask_dla],self.flux[mask*mask_dla]/dla_profile[mask*mask_dla],color='red')
                    # #plt.axvline(z_dla)
                    # plt.show()
                    # #sys.exit()
                    # if zidx==2:
                    #     sys.exit()
                    #print(str(np.mean(self.flux[mask*mask_dla]/dla_profile[mask*mask_dla]))+' '+str(np.mean(self.flux[mask*mask_dla])))

    def get_new_forest(self,rescale_flux,wrange): # line 26 in main.py
        """
        Basically what is going on is that I extrapolating the values of flux for the beta (or OVI,SiII,etc.) forest
        via its relationship with the lya forest wavelength and flux.

        We know how the wavelength relate:
        lambda_beta = lambda_alpha*lyb_rest/lya_rest (there are similar relationships for OVI and SiIII)

        And the fluxes relate:
        flux_beta = flux_alpha**(cross_section_beta/cross_section_alpha)
        """
        wmin,wmax,wrest,wxs = wrange
        obmax = wmax*(1+self.redshift)
        obmin = wmin*(1+self.redshift)
        oabmax = obmax*opt.lya_rest/wrest #lyb max converted to lya wavelength
        mask_abf = (self.wavelength>obmin)&(self.wavelength<oabmax)&(self.flux>opt.min_trans) # alpha wavelength (large mask)
        mask_bf = (self.wavelength>obmin)&(self.wavelength<obmax)&(self.flux>opt.min_trans)
        wave_abf = self.wavelength[mask_abf]
        tmp_wave_abf = wave_abf*wrest/opt.lya_rest #beta wavelength but incorrect grid
        flux_abf = self.flux[mask_abf]  #alpha fluxes (larger mask)
        tmp_flux_abf = flux_abf**(wxs/opt.xs_alpha)
        tmp_err_flux_abf = self.err_flux[mask_abf]
        wave_bf = self.wavelength[mask_bf]

        flux_bf = griddata(points=tmp_wave_abf, values=tmp_flux_abf, xi=wave_bf, method='linear')
        err_flux_bf = griddata(points=tmp_wave_abf, values=tmp_err_flux_abf, xi=wave_bf, method='linear')
        flux_bf = flux_bf**(1/rescale_flux) #july20
        err_flux_bf = err_flux_bf**(1/rescale_flux) #july20
        self.flux[mask_bf] = flux_bf*self.flux[mask_bf] #now fluxes are on correct grid and we can get our new forest
    def get_zpix(self,wave_rest):
        """
        Given the rest wavelength of the desired transition, this returns the corresponding redshift of pixels for the wavelength array of the quasar.
        """
        return self.wavelength/wave_rest - 1
    def get_zmask(self,forest,zpix,zidx,zedges,name):
        """
        Masks a quasar based on the absorption line (lya,lyb etc), redshift bin.
        """
        (rwave_min,rwave_max,rwave_rest) = forest
        owave_min = rwave_min*(1+self.redshift)
        owave_max = rwave_max*(1+self.redshift)
        owave_rest = rwave_rest*(1+self.redshift)
        mask = ((self.wavelength>owave_min)&(self.wavelength<owave_max)&
                (self.flux>opt.min_trans)&(zpix>=zedges[zidx])&(zpix<zedges[zidx+1]))

        #Use only resolution of coming from one arm. Either UV or VIS. The carswell resolutions are 17 or 10 while the original resolutions are 20 or 11...
        #for UV and VIS respectively
        if inis.onlyVIS:
            if zidx==3: #feb6 #feb11 commented again
                ### FOR REFERENCE
                # Pixels lost by just using VIS arm in z=3.6: 56%
                # Pixels lost by just using UV arm in z=3.6: 85%
                ###
                #res_mask = self.wavelength >= opt.overlap_maxwav #only use the VIS region. If you want to do the UV you must manually change code.
                mask = ((self.wavelength>owave_min)&(self.wavelength<owave_max)&
                    (self.flux>opt.min_trans)&(zpix>=zedges[zidx])&(zpix<zedges[zidx+1])&
                    (self.wavelength >= opt.overlap_maxwav))
        elif inis.no_overlap:
            if zidx==3:
                res_mask = (self.wavelength>opt.overlap_maxwav)|(self.wavelength<opt.overlap_minwav)
                mask = ((self.wavelength>owave_min)&(self.wavelength<owave_max)&
                    (self.flux>opt.min_trans)&(zpix>=zedges[zidx])&(zpix<zedges[zidx+1])&
                    (res_mask))

 
        #print("self.mask_dla ", self.mask_dla, np.sum(mask), np.sum(self.mask_dla*mask))
        return mask*self.mask_dla #july 16
        

    @staticmethod
    def get_npow(mf,nvar,dloglambda): #line 121 main.py
        """
        Computes noise power.
        mf**(-2) * nvar * dv
        """
        mf,nvar,dloglambda = np.array(mf),np.array(nvar),np.array(dloglambda)
        dv = opt.c_kms*np.log(10)*dloglambda
        return mf**(-2) * nvar * dv #np.pi / k_nyquist_lya

    def get_autopower(self,mf,mask,tag = "None"):
        """
        Computes auto-power spectrum.
        """
        rf_x = self.flux[mask]/mf - 1 #relative flux fluctuations
        dv = opt.c_kms*self.dloglambda[mask]*np.log(10) #delta v (km/s)
        rf_k = np.fft.fft(rf_x)*dv # relative fft of flux fluctuations
        N = len(rf_k)
        V = N*dv #one dimensional volume
        dk = 2*np.pi/V
        k = dk*np.arange(0,N,1)
        pk = np.abs(rf_k)**2/V
        return k,pk
    @staticmethod
    def get_kmask(kpix,kidx,kedges, kmax): # line 190 main.py
        """
        Mask to get k-bins.
        """
        mask = ((kpix>=kedges[kidx])&(kpix<kedges[kidx+1])) & (kpix<kmax) #&(kpix!=0))
        return mask
        
    def get_pk_subsets(self,kpix,pk,zmask,kmask,corr_tag,npow):
        """
        Slices previously computed power spectrum by k-bin (specified by kmask) and (potentially) applies noise and resolution correction
        """
        kpix_sub,pk_sub = kpix[kmask], pk[kmask]
        R = self.resolution[zmask][kmask] 
        dv_array = opt.c_kms*np.log(10)*self.dloglambda[zmask][kmask]   #Matt M: COME back and understand this

        #print("resolution = ", np.mean(R))
        if (('corrNR' in corr_tag) or
            ('corrN' in corr_tag) or
            ('wN' in corr_tag) or
            ('wNR' in corr_tag)):
            pk_sub = pk_sub - npow*np.sinc(0.5*kpix_sub*dv_array/np.pi)**2 #noise correction   #Matt M:  subtracting noise power prior to correcting for resolution (next step), which is correct order, I added this
        if (('corrNR' in corr_tag) or
            ('corrR' in corr_tag) or
            ('wR' in corr_tag) or
            ('wNR' in corr_tag)):
            pk_sub = pk_sub / opt.window(k=kpix_sub,p=dv_array,R=R)**2 #resolution correction
        return pk_sub
    @staticmethod
    def get_xpk_subsets(kpix,pab,qab,dlam,resa,resb,corr_tag,npow,kmask): # 247 main.py
        """
        Slices previously computed CROSS power spectrum by k-bin (specified by kmask) and (potentially) applies noise and resolution correction
        """
        kpix_sub,pab_sub,qab_sub = kpix[kmask], pab[kmask], qab[kmask]
        dv_array = opt.c_kms*np.log(10)*dlam[kmask]
        #if (('corrNR' in corr_tag) or
        #    ('corrN' in corr_tag) or
        #    ('wN' in corr_tag) or
        #    ('wNR' in corr_tag)):
        #    pab_sub = pab_sub - npow
        #    qab_sub = qab_sub - npow

        if (('corrNR' in corr_tag) or
            ('corrR' in corr_tag) or
            ('wR' in corr_tag) or
            ('wNR' in corr_tag)):
            pab_sub = pab_sub / opt.window(k=kpix_sub,p=dv_array,R=resa[kmask]
                               )/opt.window(k=kpix_sub,p=dv_array,R=resb[kmask])
            qab_sub = qab_sub / opt.window(k=kpix_sub,p=dv_array,R=resa[kmask]
                              )/opt.window(k=kpix_sub,p=dv_array,R=resb[kmask])
        return pab_sub,qab_sub
    def cross_pk_fft(self,mask_lya,mask_lyb,mf_lya,mf_lyb):
        """
        Computes cross-power spectrum.
        """
        za = self.wavelength[mask_lya]/opt.lya_rest-1
        zb = self.wavelength[mask_lyb]/opt.lyb_rest-1
        rfa = self.flux[mask_lya]/mf_lya -1
        rfb = self.flux[mask_lyb]/mf_lyb -1

        ferra = self.err_flux[mask_lya]
        ferrb = self.err_flux[mask_lyb]

        dlama = self.dloglambda[mask_lya]
        dlamb = self.dloglambda[mask_lyb]
        resa = self.resolution[mask_lya]
        resb = self.resolution[mask_lyb]
        new_af_mask = (za>np.min(zb))&(za<np.max(zb))
        new_bf_mask = (zb>np.min(za))&(zb<np.max(za))
        
        #new_bf_mask = (zb>np.min(zb))&(zb<np.max(zb))
        #print("dv = ", zb[1], 3e5*(np.min(za[new_af_mask])-za[0])/(1+za[1]), 3e5*(zb[1]-np.min(zb[new_bf_mask]))/(1+zb[0]))
        print("v offset", (np.min(za[new_af_mask])- np.min(zb[new_bf_mask]))/(1+np.min(za[new_af_mask]))*3e5+7.25)
        
        N_a = np.sum(new_af_mask)#aug5
        N_b = np.sum(new_bf_mask)#aug5
        if N_a>N_b: # beta is reference
            dlam = np.copy(dlamb[new_bf_mask])
            dv = opt.c_kms * dlam * np.log(10)
            rfa_k = np.fft.fft(rfa[new_af_mask])[:N_b]*dv
            rfb_k = np.fft.fft(rfb[new_bf_mask])*dv
            resa = resa[new_af_mask][:N_b]
            resb = resb[new_bf_mask]
        elif N_a<=N_b: # alpha is reference
            dlam = np.copy(dlama[new_af_mask])
            dv = opt.c_kms * dlam * np.log(10)
            rfa_k = np.fft.fft(rfa[new_af_mask])*dv
            rfb_k = np.fft.fft(rfb[new_bf_mask])[:N_a]*dv
            resa = resa[new_af_mask]
            resb = resb[new_bf_mask][:N_a]

        N=len(rfa_k)
        V = N*dv
        dk = 2*np.pi/(N*dv)
        k=dk*np.arange(0,N,1)
        P_cross = np.conj(rfa_k)*rfb_k/V
        P_real,P_imag = P_cross.real,P_cross.imag
        return k,P_real,P_imag,dlam,resa,resb

    # @staticmethod
    # def grid_interp(z1,z2,rf1,rf2):
    #     #which_one=1 # 1 means changing z1 & rf1, so z2 and rf2 are the references
    #     #              0 means changing z2 & rf2 so z1 and rf1 are the references
    #     #print(len(z1),len(z2))
    #     if (np.min(z1)<=np.min(z2))&(np.max(z1)>=np.max(z2)):
    #         changing_this = 0
    #         ref_grid, other_pts, other_data = z1,z2,rf2
    #     elif (np.min(z2)<np.min(z1))&(np.max(z2)>np.max(z1)):
    #         changing_this = 1
    #         ref_grid, other_pts, other_data = z2,z1,rf1
    #     elif len(z1)<=len(z2):
    #         changing_this = 0
    #         ref_grid, other_pts, other_data = z1,z2,rf2
    #     elif len(z1)>len(z2):
    #         changing_this = 1
    #         ref_grid, other_pts, other_data = z2,z1,rf1
    #     else:
    #         changing_this = 0
    #         ref_grid, other_pts, other_data = z1,z2,rf2
    #
    #     new_rf = griddata(other_pts,other_data,ref_grid,method='linear')
    #     return ref_grid,new_rf,changing_this

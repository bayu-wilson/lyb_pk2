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
    # all_dla_NHI = None
    # all_dla_names = None
    # all_dla_redshifts = None
    dla_table = None
    nqso = None
    """Container for quasar spectra"""
    def __init__(self,wavelength=None, flux=None, err_flux=None,
                 resolution=None, dloglambda=None, name=None, redshift=None,
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
        #self.mask_dla = mask_dla

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
            ax.set_ylabel(r"Flux [$10^{-15}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$]",fontsize=fontsize)
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
        if tag.startswith('obs'): #What matters is what the tag starts/ends with
            path_to_cat = "../data/obs/XQ-100_catalogue.txt"
            catalog_table = pd.read_csv(path_to_cat,delim_whitespace=True,usecols=(3,6),
                                names = ("qso_name","redshifts"))
            path_to_dlas = "../data/obs/XQ-100_DLA_catalogue.txt"
            dla_colnames = ['idx','name','z','NHI','err_NHI']
            ####
            dla_table = pd.read_csv(path_to_dlas,delim_whitespace=True,names=dla_colnames)
            qname_array,z_array = catalog_table['qso_name'].values,catalog_table['redshifts'].values
            nqso = 100
            cls.dla_table = dla_table#.iloc[26:27]
            #print()
            # cls.dla_table = pd.DataFrame([dla_table.iloc[i] for i in [38,39,40]])
            #[2,3,7,8,9,14,18,35]])
            #dla_table.iloc[]
            #.iloc[8:9] #dla_table#.iloc[:1] # CHANGE NOV 19 !!!!!!!!!!!!!!!!!!!!!!!!
            #print(dla_table.iloc[8:9])
            #cls.all_dla_NHI = dla_table.NHI.values
            #cls.all_dla_names = dla_table.name.values
            #cls.all_dla_redshifts = dla_table.z.values
        elif tag.startswith('mocks'): #What matters is what the tag starts/ends with
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
        #glob.glob("../data/")

        if name.startswith("mock"):
            all_paths = glob.glob("../data/mocks/XQ-100_{0}/ly*/ref/spectra/mock-*.txt.gz".format(
                                    tag))
            all_paths.sort()
            path = all_paths[cat_index]
            #path = "../data/mocks/XQ-100_{0}/lyb/ref/spectra/{1}_xq{2}_{0}.txt.gz".format(tag,name,nqso)
            qso_table = pd.read_csv(path, delim_whitespace=True,compression='gzip',skiprows = 1,usecols=(0,1,2,3,4),
                                    names = ("wav", "flx", "ferr", "res","dloglam"),float_precision='high')
#             print(qso_table.flx>opt.min_trans)
            mask_wave = ((qso_table.wav>opt.lyb_min*(1+redshift))&
                        (qso_table.wav<opt.lya_rest/opt.lyb_rest*opt.lyb_max*(1+redshift))&
                        (qso_table.flx>opt.min_trans))
            qso_table=qso_table[mask_wave]
            wavelength = qso_table.wav.values
            flux = qso_table.flx.values
            err_flux = qso_table.ferr.values
            unnormalized_flux=np.ones_like(qso_table.flx.values)*np.nan
            resolution = qso_table.res.values # june 24
        else:
            if inis.redside_avg == True: #feb13
                released_path = "../data/obs/XQ-100/released/{0}_uvb-vis.txt".format(name) #feb13
                qso_table = pd.read_csv(released_path, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4), #feb13
                                        names = ("wav", "flx", "ferr", "res","dloglam")) #feb13
                mask_wave = (qso_table.wav.values>opt.lyb_min*(1+redshift))&( #feb13
                           qso_table.wav.values<opt.lya_rest/opt.lyb_rest*opt.lyb_max*(1+redshift))&( #feb13
                           qso_table.flx.values>opt.min_trans) #feb13
                qso_table = qso_table[mask_wave] #feb13
                wavelength = qso_table.wav.values #feb13
                unnormalized_flux = qso_table.flx.values #feb13
                flux = qso_table.flx.values*1
                err_flux = qso_table.ferr.values #feb13
                resolution = np.ones_like(wavelength)
                m1 = wavelength<=5599.14
                resolution[m1] = 41.52075368/(2*np.sqrt(2*np.log(2))) #converting FWHM to sigma_R
                resolution[~m1] = 23.60134842/(2*np.sqrt(2*np.log(2)))

                ### normalized flux with redside_avg ###
                zpix = wavelength/opt.lya_rest - 1 # 3.6,4.8
                redside_avg = np.loadtxt("../output/redside_avg.txt")[:-2] #only using first 4
                temp_edges = np.array([2.0,3.7,3.9,4.1,5.0])
                for i in range(len(temp_edges)-1):
                    mask = (zpix>temp_edges[i])&(zpix<temp_edges[i+1])
                    flux[mask] = unnormalized_flux[mask]/redside_avg[i]
                    err_flux[mask] = err_flux[mask]/redside_avg[i]

            else:
                released_path = "../data/obs/XQ-100/released/{0}_uvb-vis.txt".format(name)
                continuum_path = "../data/obs/XQ-100/continuum/{0}_cont.txt".format(name)
                qso_table = pd.read_csv(released_path, delim_whitespace=True, skiprows = 1,usecols=(0,1,2,3,4),
                                        names = ("wav", "flx", "ferr", "res","dloglam"))
                #print(qso_table)
                cont_table = pd.read_csv(continuum_path, delim_whitespace=True, skiprows = 1,usecols=(0,1),
                                        names = ("wav",'flx'))
                mask_wave = (qso_table.wav.values>opt.lyb_min*(1+redshift))&(
                           qso_table.wav.values<opt.lya_rest/opt.lyb_rest*opt.lyb_max*(1+redshift))&(
                           qso_table.flx.values>opt.min_trans)
                qso_table = qso_table[mask_wave]
                cont_table = cont_table[mask_wave]
                wavelength = qso_table.wav.values
                unnormalized_flux = qso_table.flx.values
                flux = unnormalized_flux/cont_table.flx.values
                err_flux = qso_table.ferr.values/cont_table.flx.values

                # resolution = qso_table.res.values #jan7 #RESOLUTION FROM FILE #RMS RESOLUTION (sigma_R) = FWHM/(2*sqrt(2ln(2)))
                # #June 25th #Carswell
                resolution = np.ones_like(wavelength)
                m1 = wavelength<=5599.14
                resolution[m1] = 41.52075368/(2*np.sqrt(2*np.log(2))) #converting FWHM to sigma_R
                resolution[~m1] = 23.60134842/(2*np.sqrt(2*np.log(2)))
        if "wR2" in tag and not inis.wR2:
            # Using new column
            resolution = np.ones_like(qso_table.res.values)*11 * 0.2 #/ 2.3#(2*np.sqrt(2*np.log(2)))
        dloglambda = qso_table.dloglam.values
        if ("noB" in tag)&(inis.add_beta): #july21
            flux = flux**rescale_flux #np.exp(rescale_flux*np.log(flux)) #flux**rescale_flux #july20
        return cls(name=name,redshift=redshift,
                   wavelength=wavelength, flux=flux, err_flux=err_flux,
                   resolution=resolution,dloglambda=dloglambda,
                   unnormalized_flux=unnormalized_flux)

    def get_chunks(self,which_chunk): # askdfjaskdjf # am I even using this!?!?!?!?
        increases,decreases = 0,0
        inc_list, dec_list = [],[]
        m = self.mask_dla
        m_idx = np.arange(0,len(m),1)
        for i in range(len(m)-1):
            if m[i+1]>m[i]:
                increases += 1
                inc_list.append(i+1) # choose higher (True) one
            if m[i]>m[i+1]:
                decreases += 1
                dec_list.append(i) # choose higher (True) one
        #print("\ninc_list:",inc_list,"\ndec_list:",dec_list,'\n')

        if (len(inc_list) == 0)&(len(dec_list) == 0):
            highest_idx = -1
            lowest_idx = 0
            second_lowest = highest_idx

        elif len(inc_list) == 0: # dla on right end only
            lowest_idx = 0 # start on first index
            highest_idx = dec_list[-1]
            second_lowest = highest_idx #dec_list[0]
        elif len(dec_list) == 0: #dla on left end only
            lowest_idx = inc_list[0]
            highest_idx = -1
            second_lowest = highest_idx
        else: # if lowest increase is greater than lowest decrese ... most typical
            second_lowest = dec_list[0]
            if inc_list[0]<dec_list[0]: #lowest increase less than lowest decrease ... a dla on left edge
                lowest_idx = inc_list[0]
                highest_idx = -1 # !
            if inc_list[-1]<dec_list[-1]: # highest increase is less than highest decrease ... a dla on left edge
                lowest_idx = 0 # !
                highest_idx = dec_list[-1]
            else:
                lowest_idx = 0
                highest_idx = -1 # !


        nsplits = np.min([increases,decreases])
        mask_list = [m_idx[lowest_idx:second_lowest]]
        #print(nsplits,"\ninc_list:",inc_list,"\ndec_list:",dec_list,'\n')
        for i in range(nsplits-1):
            mask_list.append(m_idx[inc_list[i]:dec_list[i+1]])
        #if len(inc_list) >=2:
        if len(inc_list) > 0:
            mask_list.append(m_idx[inc_list[-1]:highest_idx])
        #else:

        # print(0,dec_list[0],inc_list[0],dec_list[1],dec_list[1],-1)
        # mask_list = [m[0:dec_list[0]],
        #          m[inc_list[0]:dec_list[1]],
        #          m[dec_list[1]:-1]]
        # print('a:',inc_list)
        # print('b:',dec_list)
        # print('\n')
        mask = mask_list[which_chunk]

        #mask = np.where(m[dec_list[which_chunk]:inc_list[which_chunk]])[0]

        #mask = mask[which_chunk] # asdifhasiodjfhasdf
        self.wavelength = self.wavelength[mask]
        self.flux = self.flux[mask]
        self.err_flux = self.err_flux[mask]
        self.resolution = self.resolution[mask]
        self.dloglambda = self.dloglambda[mask]
        self.unnormalized_flux = self.unnormalized_flux[mask]
        #self.mask_dla = mask_dla[mask]


    def get_new_forest(self,rescale_flux,wrange): # line 26 in main.py
        #print("1")
        wmin,wmax,wrest,wxs = wrange
        obmax = wmax*(1+self.redshift)
        obmin = wmin*(1+self.redshift)
        # lyb max converted to lya wavelength
        oabmax = obmax*opt.lya_rest/wrest
        # alpha wavelength (large mask)
        mask_abf = (self.wavelength>obmin)&(self.wavelength<oabmax)&(self.flux>opt.min_trans)
        mask_bf = (self.wavelength>obmin)&(self.wavelength<obmax)&(self.flux>opt.min_trans)
        wave_abf = self.wavelength[mask_abf]
        # beta wavelength but incorrect grid
        tmp_wave_abf = wave_abf*wrest/opt.lya_rest
        flux_abf = self.flux[mask_abf]  #alpha fluxes (larger mask)
        tmp_flux_abf = flux_abf**(wxs/opt.xs_alpha) #np.exp(wxs/opt.xs_alpha*np.log(flux_abf)) #july20
        #flux_abf**(opt.xs_beta/opt.xs_alpha)
        tmp_err_flux_abf = self.err_flux[mask_abf]
        # alpha flux (incorrect grid)
        wave_bf = self.wavelength[mask_bf]
        ### LOOK HERE ###
        flux_bf = griddata(points=tmp_wave_abf, values=tmp_flux_abf, xi=wave_bf, method='linear')
        err_flux_bf = griddata(points=tmp_wave_abf, values=tmp_err_flux_abf, xi=wave_bf, method='linear')
        #flux_bf_err = griddata(points=tmp_wave_abf, values=tmp_flux_abf, xi=wave_bf, method='linear')
        #**(1/rescale_lyb)
        flux_bf = flux_bf**(1/rescale_flux)#np.exp((1/rescale_flux)*np.log(flux_bf)) #july20
        err_flux_bf = err_flux_bf**(1/rescale_flux) #np.exp((1/rescale_flux)*np.log(err_flux_bf)) #july20
        #now fluxes are on correct grid
        #print(flux_bf)
        #flux_bf = flux_bf*self.flux[mask_bf]
        #zpix = wave_bf/opt.lyb_rest-1
        #zmask = (zpix>=zedges[zidx])&(zpix<zedges[zidx+1])
        #flux_bf = flux_bf**(rescale_lyb)
        ### LOOK HERE ###
        self.flux[mask_bf] = flux_bf*self.flux[mask_bf]
    def get_zpix(self,wave_rest):
        return self.wavelength/wave_rest - 1
    def get_zmask(self,forest,zpix,zidx,zedges,name):
        (rwave_min,rwave_max,rwave_rest) = forest #
        owave_min = rwave_min*(1+self.redshift)
        owave_max = rwave_max*(1+self.redshift)
        owave_rest = rwave_rest*(1+self.redshift)
        mask = ((self.wavelength>owave_min)&(self.wavelength<owave_max)&
                (self.flux>opt.min_trans)&(zpix>=zedges[zidx])&(zpix<zedges[zidx+1]))
        # if zidx==3: #feb6 #feb11 commented again
        #     only_res = 20 #41.52075368/(2*np.sqrt(2*np.log(2))) # 20 means UV arm. 11 means VIS arm.
        #     #mask = self.resolution==only_res
        #     mask = ((self.wavelength>owave_min)&(self.wavelength<owave_max)&
        #         (self.flux>opt.min_trans)&(zpix>=zedges[zidx])&(zpix<zedges[zidx+1])&
        #         (self.resolution==only_res))



        if inis.cat_name.startswith('obs'):
            if (name in self.dla_table.name.values)&(inis.remove_dla):
                # print(name) #dec2
                sub_table = self.dla_table[self.dla_table.name == name]
                num_dlas = len(sub_table)
                # if num_dlas>1:
                #     print(name)
                mask_dla = np.ones_like(self.wavelength,'bool')
                for n in range(num_dlas):
                    row = sub_table.iloc[n]
                    z_dla,NHI = row.z, row.NHI
                    wave_dla = (1+z_dla)*opt.lya_rest
                    wave_dlb = (1+z_dla)*opt.lyb_rest
                    lo_wave_dla = wave_dla-opt.find_EW(NHI,'alpha')*1.5/2
                    hi_wave_dla = wave_dla+opt.find_EW(NHI,'alpha')*1.5/2
                    lo_wave_dlb = wave_dlb-opt.find_EW(NHI,'beta')*1.5/2
                    hi_wave_dlb = wave_dlb+opt.find_EW(NHI,'beta')*1.5/2
                    mask_dla*=(self.wavelength<lo_wave_dla)|(self.wavelength>hi_wave_dla)
                    mask_dla*=(self.wavelength<lo_wave_dlb)|(self.wavelength>hi_wave_dlb)
                    # lo_wave_lyaf = (1+self.redshift)*opt.lya_min
                    # hi_wave_lyaf = (1+self.redshift)*opt.lya_max
                    # lo_wave_lybf = (1+self.redshift)*opt.lyb_min
                    # hi_wave_lybf = (1+self.redshift)*opt.lyb_max
                    # mask_dla *= (
                    #  (self.wavelength<lo_wave_dla)&(self.wavelength>hi_wave_dla))|(
                    # ((self.wavelength<lo_wave_dlb)&(self.wavelength>hi_wave_dlb)))
                # print(sum(mask_dla))
                #print(lo_wave_dla,lo_wave_dlb,hi_wave_dla,hi_wave_dlb,min(self.wavelength), max(self.wavelength))
                #print(opt.how_many_chunks(mask*mask_dla))
                return mask*mask_dla #july 16
        return mask

    @staticmethod
    def get_npow(mf,nvar,dloglambda): #line 121 main.py
        mf,nvar,dloglambda = np.array(mf),np.array(nvar),np.array(dloglambda)
        dv = opt.c_kms*np.log(10)*dloglambda
        return mf**(-2) * nvar * dv #np.pi / k_nyquist_lya
    # def get_rf(self,mf):
    #     return self.flux/mf - 1
    def get_autopower(self,mf,mask,tag = "None"):
        #if 'noB' in tag.split('_'):
        #    self.flux[mask]
        #    #self.get_lybforest(zidx,zedges, rescale_lyb = 1/0.1)
        rf_x = self.flux[mask]/mf - 1
        dv = opt.c_kms*self.dloglambda[mask]*np.log(10)
        # print(self.dloglambda[mask])
        # import sys
        # sys.exit()
        rf_k = np.fft.fft(rf_x)*dv
        N = len(rf_k)
        V = N*dv
        dk = 2*np.pi/V
        k = dk*np.arange(0,N,1)
        pk = np.abs(rf_k)**2/V
        return k,pk
    @staticmethod
    def get_kmask(kpix,kidx,kedges): # line 190 main.py
        mask = ((kpix>=kedges[kidx])&(kpix<kedges[kidx+1]))#&(kpix!=0))
         #(self.flux>opt.min_trans)&(zpix>zedges[zidx])&(zpix<=zedges[zidx+1]))
        return mask
    def get_pk_subsets(self,kpix,pk,zmask,kmask,corr_tag,npow):
        kpix_sub,pk_sub = kpix[kmask], pk[kmask]
        R = self.resolution[zmask][kmask]
        dv_array = opt.c_kms*np.log(10)*self.dloglambda[zmask][kmask]
        if (('corrNR' in corr_tag) or
            ('corrN' in corr_tag) or
            ('wN' in corr_tag) or
            ('wNR' in corr_tag)):
            pk_sub = pk_sub - npow
        if (('corrNR' in corr_tag) or
            ('corrR' in corr_tag) or
            ('wR' in corr_tag) or
            ('wNR' in corr_tag)):
            pk_sub = pk_sub / opt.window(k=kpix_sub,p=dv_array,R=R)**2
        return pk_sub
    @staticmethod
    def get_xpk_subsets(kpix,pab,qab,dlam,resa,resb,corr_tag,npow,kmask): # 247 main.py
        kpix_sub,pab_sub,qab_sub = kpix[kmask], pab[kmask], qab[kmask]
        dv_array = opt.c_kms*np.log(10)*dlam[kmask]
        if (('corrNR' in corr_tag) or
            ('corrN' in corr_tag) or
            ('wN' in corr_tag) or
            ('wNR' in corr_tag)):
            pab_sub = pab_sub - npow
            qab_sub = qab_sub - npow

        if (('corrNR' in corr_tag) or
            ('corrR' in corr_tag) or
            ('wR' in corr_tag) or
            ('wNR' in corr_tag)):
            pab_sub = pab_sub / opt.window(k=kpix_sub,p=dv_array,R=resa[kmask]
                               )/opt.window(k=kpix_sub,p=dv_array,R=resb[kmask])
            qab_sub = qab_sub / opt.window(k=kpix_sub,p=dv_array,R=resa[kmask]
                              )/opt.window(k=kpix_sub,p=dv_array,R=resb[kmask])
            # import sys
            # print(len(resa),len(resb),np.sum(kmask))
            # print(resa[kmask])
            # print(resb[kmask])
            # sys.exit()

        return pab_sub,qab_sub
    def cross_pk_fft(self,mask_lya,mask_lyb,mf_lya,mf_lyb):
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

        ### DEBUG #aug5
        # x,y,which = self.grid_interp(za[new_af_mask],zb[new_bf_mask],
        #                              rfa[new_af_mask],rfb[new_bf_mask]) #reference, changed, which one
        #
        # finite_mask = np.isfinite(y)
        # if which == 1: #lyb is reference -- confirmed
        #
        #     zb,rfb,ferrb,resb = [i[new_bf_mask][finite_mask] for i in [zb,rfb,ferrb,resb]]
        #     resa = griddata(za,resa,zb)
        #     za,rfa,dlam = (x[finite_mask],y[finite_mask],
        #                         dlamb[new_bf_mask][finite_mask])
        # elif which == 0: #lya is reference -- confirmed
        #     za,rfa,ferra,resa = [i[new_af_mask][finite_mask] for i in [za,rfa,ferra,resa]]
        #     resb = griddata(zb,resb,za)
        #     zb,rfb,dlam = (x[finite_mask],
        #                        y[finite_mask],
        #                        dlama[new_af_mask][finite_mask])
        ### DEBUG #aug5


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

            #rfb = rfb[new_bf_mask]
            #rfa = rfa[:N_b]#aug5
            #dlam = np.copy(dlamb[new_bf_mask])#aug5
        # elif N_a<=N_b: # alpha is reference
        #     rfa = rfa[new_af_mask]
        #     resa = resa[new_af_mask]
        #     rfb = rfb[:N_a]#aug5
        #     resb = resb[:N_a]
        #     dlam = np.copy(dlama[new_af_mask])#aug5
        #dv = opt.c_kms * dlam * np.log(10)
        #rfa_k = np.fft.fft(rfa)*dv #dloglam
        #rfb_k = np.fft.fft(rfb)*dv #dloglam

        N=len(rfa_k)
        V = N*dv
        dk = 2*np.pi/(N*dv)
        k=dk*np.arange(0,N,1)
        P_cross = np.conj(rfa_k)*rfb_k/V
        P_real,P_imag = P_cross.real,P_cross.imag
        return k,P_real,P_imag,dlam,resa,resb
        # else:
        #     return [np.nan],[np.nan],[np.nan],[np.nan],[np.nan],[np.nan]

    @staticmethod
    def grid_interp(z1,z2,rf1,rf2):
        #which_one=1 # 1 means changing z1 & rf1, so z2 and rf2 are the references
        #              0 means changing z2 & rf2 so z1 and rf1 are the references
        #print(len(z1),len(z2))
        if (np.min(z1)<=np.min(z2))&(np.max(z1)>=np.max(z2)):
            changing_this = 0
            ref_grid, other_pts, other_data = z1,z2,rf2
        elif (np.min(z2)<np.min(z1))&(np.max(z2)>np.max(z1)):
            changing_this = 1
            ref_grid, other_pts, other_data = z2,z1,rf1
        elif len(z1)<=len(z2):
            changing_this = 0
            ref_grid, other_pts, other_data = z1,z2,rf2
        elif len(z1)>len(z2):
            changing_this = 1
            ref_grid, other_pts, other_data = z2,z1,rf1
        else:
            changing_this = 0
            ref_grid, other_pts, other_data = z1,z2,rf2

        new_rf = griddata(other_pts,other_data,ref_grid,method='linear')
        return ref_grid,new_rf,changing_this

                # z_ewa = opt.find_EW(NHI,'alpha')/2/opt.lya_rest
                # z_ewb = opt.find_EW(NHI,'beta')/2/opt.lyb_rest
                # lo_z_dla = z_dla - z_ewa
                # hi_z_dla = z_dla + z_ewa
                # lo_z_dlb = z_dla - z_ewb
                # hi_z_dlb = z_dla + z_ewb
            #print(num_dlas)
        #z_dla = "asdf"
        # lo_wave_dla = wave_dla-opt.find_EW(NHI_dla,'alpha')*1.5/2
        # hi_wave_dla = wave_dla+opt.find_EW(NHI_dla,'alpha')*1.5/2
        # lo_wave_dlb = wave_dlb-opt.find_EW(NHI_dla,'beta')*1.5/2
        # hi_wave_dlb = wave_dlb+opt.find_EW(NHI_dla,'beta')*1.5/2
        ### CHECK FOR DLAs ###
        #dla_table = cls.dla_table


        #mask = np.logical_and.reduce(((self.wavelength>owave_min),(self.wavelength<owave_max),
        # (self.flux>opt.min_trans),(zpix>=zedges[zidx]),(zpix<zedges[zidx+1])))
        #return mask
    # def get_forest_mask(self,forest):
        # (rwave_min,rwave_max,rwave_rest) = forest #
        # owave_min = rwave_min*(1+self.redshift)
        # owave_max = rwave_max*(1+self.redshift)
        # owave_rest = rwave_rest*(1+self.redshift)
        # mask = ((self.wavelength>owave_min)&(self.wavelength<owave_max)&
        #  (self.flux>opt.min_trans))
        # return mask
    # def get_nvar(self):
        # return self.err_flux**2

        # if len(dec_list) == 0: #dla on left end only
        #     lowest_idx = inc_list[0]
        #     highest_idx = -1
        # elif inc_list[-1]<dec_list[-1]: # highest increase is less than highest decrease ... a dla on left edge
        #     highest_idx = dec_list[-1]
        # else:
        #     highest_idx = -1

        #if
            # #if inc_list[0]<dec_list[0]: #
            # #    lowest_idx = inc_list[0]
            #
            # if len(dec_list) == 0: # dla on left end
            #     highest_idx = -1
            # else:
            #     highest_idx = dec_list[-1]


        # elif len(dec_list) == 0:
        #     lowest_idx = inc_list[0]
        #     highest_idx = -1
        # elif inc_list[0]<dec_list[0]:
        #     lowest_idx = inc_list[0]
        # else:
        #     lowest_idx = 0

        # if len(inc_list) == 0:
        #     highest_idx = dec_list[-1]
        # elif len(dec_list) == 0:
        #     highest_idx = -1
        # elif inc_list[-1]<dec_list[-1]:
        #     highest_idx = dec_list[-1]
        # else:
        #     highest_idx = -1#np.where(m)[0][-1]

    # def get_autopower2(self,flux,dloglambda,mf,tag = "None"):
    #     rf_x = flux/mf - 1
    #     dv = opt.c_kms*dloglambda*np.log(10)
    #     rf_k = np.fft.fft(rf_x)*dv
    #     N = len(rf_k)
    #     V = N*dv
    #     dk = 2*np.pi/V
    #     k = dk*np.arange(0,N,1)
    #     pk = np.abs(rf_k)**2/V
    #     return k,pk

        # resolution = qso_table.res.values
        #
        # dla_table = cls.dla_table
        #
        #
        # # if zidx == 6:
        # #     print(z,id,z_dla)
        # if (name in dla_table.name.values)&(inis.remove_dla)&(name!="J0034+1639"):
        #     dla_per_qso = dla_table[dla_table.name.values == name]
        #     mask_dla = np.ones_like(wavelength,dtype=bool)
        #     #mask_dla = []
        #     #print()
        #
        #     for i in range(len(dla_per_qso)):
        #         NHI_dla = dla_per_qso.iloc[i].NHI
        #         z_dla = dla_per_qso.iloc[i].z
        #         wave_dla = (1+z_dla)*opt.lya_rest
        #         wave_dlb = (1+z_dla)*opt.lyb_rest
        #         lo_wave_dla = wave_dla-opt.find_EW(NHI_dla,'alpha')*1.5/2
        #         hi_wave_dla = wave_dla+opt.find_EW(NHI_dla,'alpha')*1.5/2
        #         lo_wave_dlb = wave_dlb-opt.find_EW(NHI_dla,'beta')*1.5/2
        #         hi_wave_dlb = wave_dlb+opt.find_EW(NHI_dla,'beta')*1.5/2
        #         lo_wave_lyaf = (1+redshift)*opt.lya_min
        #         hi_wave_lyaf = (1+redshift)*opt.lya_max
        #         lo_wave_lybf = (1+redshift)*opt.lyb_min
        #         hi_wave_lybf = (1+redshift)*opt.lyb_max
        #         if inis.lya_dlas_in_lyaf:
        #             #if (wave_dla>lo_wave_lyaf)&(wave_dla<hi_wave_lyaf):
        #             #print(z_dla,wave_dla,lo_wave_lyaf,hi_wave_lyaf)
        #             if (wave_dla>lo_wave_lyaf)&(wave_dla<hi_wave_lyaf)|(
        #                 lo_wave_dla<hi_wave_lyaf)&(hi_wave_dla>hi_wave_lyaf)|(
        #                 hi_wave_dla>lo_wave_lyaf)&(lo_wave_dla<lo_wave_lyaf):
        #                 mask_dla*=~((wavelength>lo_wave_dla)&(wavelength<hi_wave_dla))
        #                 #print(name,np.sum((wavelength>lo_wave_dla)&(wavelength<hi_wave_dla)))
        #         if inis.lyb_dlas_in_lybf:
        #             # if (wave_dlb>lo_wave_lybf)&(wave_dlb<hi_wave_lybf)|(
        #             #     lo_wave_dlb<hi_wave_lybf)&(hi_wave_dlb>hi_wave_lybf)|(
        #             #     hi_wave_dlb>lo_wave_lybf)&(lo_wave_dlb<lo_wave_lybf):
        #             #     mask_dla*=~((wavelength>lo_wave_dlb)&(wavelength<hi_wave_dlb))
        #             if (wave_dlb>lo_wave_lybf)&(wave_dlb<hi_wave_lybf):
        #
        #                 mask_dla*=~((wavelength>lo_wave_dlb)&(wavelength<hi_wave_dlb))
        #         else:
        #             mask_dla*=~((wavelength>lo_wave_dla)&(wavelength<hi_wave_dla))
        #         #mask_dla*=~((wavelength>lo_wave_dlb)&(wavelength<hi_wave_dlb))
        #
        # else:
        #
        #     #print(name,redshift,'no dla')
        #     mask_dla = []

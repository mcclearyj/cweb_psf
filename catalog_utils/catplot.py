import unittest
test = unittest.TestCase
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import ipdb, pdb
from astropy.io import fits
from scipy.stats import chi2
import galsim
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error, normalized_root_mse
from catalogutils import AttrDict, set_rc_params
from catalog_hsm_fitter import do_hsm_fit
import catalogaugmenter as ca
import glob
import re
from astropy.table import Table, vstack, hstack, Column
import pandas as pd

class plot():
    def __init__(self, catalog_object, psf_object): #mirage == false
        self.catalog_obj = catalog_object
        self.psf_obj = psf_object
        
        substamps = []
        j = 10
        for starstamp in catalog_object.data['VIGNET']:
            starstamp[starstamp<= -999] = np.nan
            substamps.append(starstamp[-j:,-j:])
            substamps.append(starstamp[0:j,0:j])

        self.sky_level = np.nanmean(substamps)
        self.sky_std = np.nanstd(substamps)
        
        #catalog_object.add_noise_flux([psf_object], sky_level=self.sky_level, sky_std=self.sky_std)
        self.stars = catalog_object.data['VIGNET'] #formerly starmaker
        self.psfs  = catalog_object.data[psf_object.nameColumn()] #formerly psfmaker #try except with name for PSFEX_VIGNET
        self.err_stamps = catalog_object.data['ERR_VIGNET']
        catalog_object.crop([psf_object], vignet_size=75)
        catalog_object.save_new(outname=catalog_object.catalog)
        try:
            self.cropped_stars = catalog_object.data['VIGNET_CROPPED']
        except:
            self.cropped_stars = catalog_object.data['VIGNET']
        try:
            self.cropped_psfs = catalog_object.data[psf_object.nameColumn()+'_CROPPED']
        except:
            self.cropped_psfs = catalog_object.data[psf_object.nameColumn()]
        try:
            self.cropped_err_stamps = catalog_object.data['ERR_VIGNET_CROPPED']
        except:
            self.cropped_err_stamps = catalog_object.data['ERR_VIGNET']
        

class resid_plot(plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)
        self.avg_star = np.nanmean(self.cropped_stars, axis=0)
        
        psfs = self.cropped_psfs
        for i in range(len(psfs)):
            psfs[i] *= self.catalog_obj.data['FLUX_AUTO'][i]
            psfs[i] += np.random.normal(loc = self.sky_level, scale = self.sky_std, size = psfs[i].shape)

        self.cropped_psfs = psfs
        #self.avg_psf = np.nanmean(self.cropped_psfs, axis=0)
        self.avg_psf = np.nanmean(psfs, axis=0)
        #self.psf_std = np.nanstd(psfs, axis=0)
        self.sem = 0.0
        self.avg_residual = np.zeros(self.avg_star.shape)
        set_rc_params(fontsize=16)
        self.titles = []
        self.sum_residuals = 0.0
        self.fwhm_model = []
        self.fwhm_data = []
        self.cropped_size = self.cropped_stars[0].shape
        self.vignet_hsm = []
        self.psf_hsm = []
        self.wg = []


    def preprocessing(self, hsm_fit=False):
        stars = self.cropped_stars
        for i in range(len(stars)):
            #count = 0
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    if stars[i][j, k] < -1000:
                        stars[i][j, k] = np.nan
                        #count += 1
            #print('NaNs: ', count)
      #  for i in range(len(stars)):
       #     stars[i] = stars[i]/np.nansum(stars[i])
            #print(np.nansum(stars[i]))
        self.cropped_stars = stars
        self.avg_star = np.nanmean(self.cropped_stars, axis=0)
        
        psfs = self.cropped_psfs
        for i in range(len(psfs)):
            for j in range(psfs[0].shape[0]):
                for k in range(psfs[0].shape[1]):
                    if psfs[i][j, k] < -10:
                        psfs[i][j, k] = np.nan
        #for i in range(len(psfs)):
         #   psfs[i] = psfs[i]/np.nansum(psfs[i] + 1e-10)
            #print(np.nansum(psfs[i]))
        self.cropped_psfs = psfs
        self.avg_psf = np.nanmean(self.cropped_psfs, axis=0)
        if hsm_fit:
            for starstamp in stars:
                starstamp[starstamp<= -1000] = np.nan
                gs_object = galsim.Image(starstamp, wcs=galsim.PixelScale(0.3), xmin=0, ymin=0)
                try:
                    HSM_fit=gs_object.FindAdaptiveMom()
                    self.vignet_hsm.append(1)
                except:
                    try:
                        gs_object = galsim.Image(starstamp+abs(np.min(starstamp)), wcs=galsim.PixelScale(0.3), xmin=0, ymin=0)
                        HSM_fit=gs_object.FindAdaptiveMom()
                        self.vignet_hsm.append(1)
                    except:
                        self.vignet_hsm.append(0)
                        #print("Vignet HSM failed")

            for psfstamp in self.cropped_psfs:
                gs_object = galsim.Image(psfstamp, wcs=galsim.PixelScale(0.3), xmin=0, ymin=0)
                try:
                    HSM_fit=gs_object.FindAdaptiveMom()
                    self.psf_hsm.append(1)
                except:
                    try:
                        gs_object = galsim.Image(psfstamp+abs(np.min(psfstamp)), wcs=galsim.PixelScale(0.3), xmin=0, ymin=0)
                        HSM_fit=gs_object.FindAdaptiveMom()
                        self.psf_hsm.append(1)
                    except:
                        self.psf_hsm.append(0)
                        #print("PSF HSM failed")

            self.wg = [a == 1 and b == 1 for a, b in zip(self.vignet_hsm, list(self.psf_hsm))]
        else:
            self.wg = [True for i in range(len(self.cropped_stars))]
            #self.wg = [True if num == 1 else False for num in np.array(self.vignet_hsm)]
        '''
        if hsm_fit_webb:
            for starstamp in self.cropped_stars:
                gs_object = galsim.Image(starstamp, wcs=galsim.PixelScale(0.3), xmin=0, ymin=0)
                try:
                    HSM_fit=gs_object.FindAdaptiveMom()
                    self.vignet_hsm.append(1)
                except:
                    try:
                        gs_object = galsim.Image(starstamp+abs(np.min(starstamp)), wcs=galsim.PixelScale(0.3), xmin=0, ymin=0)
                        HSM_fit=gs_object.FindAdaptiveMom()
                        self.vignet_hsm.append(1)
                    except:
                        self.vignet_hsm.append(0)
                        #print("Vignet HSM failed")


            self.wg = [True if num == 1 else False for num in self.vignet_hsm]
        else:
            self.wg = [True for i in range(len(self.cropped_stars))]
        '''
        
    def return_residuals_sum(self):
        return self.sum_residuals/(self.cropped_size[0] * self.cropped_size[1])

    def set_residuals_sum(self):
        self.sum_residuals = np.nansum(self.avg_residual)

    def calc_fwhm(self, pix_scale=0.03):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        fwhm_star = []
        fwhm_psf = []
        
        for i in range(len(stars)):
            if type(stars[i]) is galsim.image.Image:
                gs_object = stars[i]
            else:
                # HSM fits fail if there are too many negative pixels
                gs_object = galsim.Image(stars[i], wcs=galsim.PixelScale(pix_scale))    
            fwhm_star.append(gs_object.calculateFWHM())

        for i in range(len(psfs)):
            if type(psfs[i]) is galsim.image.Image:
                gs_object = psfs[i]
            else:
                gs_object = galsim.Image(psfs[i], wcs=galsim.PixelScale(pix_scale))
            fwhm_psf.append(gs_object.calculateFWHM())
        
        self.fwhm_model = fwhm_star
        self.fwhm_data = fwhm_psf

    def return_fwhm_star(self):
        return self.fwhm_model

    def return_fwhm_psf(self):
        return self.fwhm_data

    def set_titles(self, titles):
        self.titles = titles

    def save_figure(self, outname):
        #vmin = -0.1 #np.nanmin(self.avg_psf)
        vmin = np.minimum(np.nanmin(self.avg_star), np.nanmin(self.avg_psf))
        #vmax = 1
        vmax = np.maximum(np.nanmax(self.avg_star), np.nanmax(self.avg_psf))
        
        vmin2 = np.nanmin(self.avg_residual)
        vmax2 = np.nanmax(self.avg_residual)
        if np.abs(vmin2) > np.abs(vmax2):
            vmax2 = np.abs(vmin2)
        else:
            vmin2 = -1*vmax2

        norm = colors.SymLogNorm(vmin=vmin, vmax=vmax, linthresh=1e-4)
        #norm2 = colors.SymLogNorm(vmin=vmin2, vmax=vmax2, linthresh=1e-4)
        norm2 = colors.TwoSlopeNorm(vcenter=0, vmin=vmin2, vmax=vmax2)
        #norm2 = colors.TwoSlopeNorm(0)
        cmap=plt.cm.turbo
        #cmap=plt.cm.BrBG
        #cmap = plt.cm.inferno
        cmap2=plt.cm.bwr_r
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=[15,7], tight_layout=True)
        
        im = axs[0].imshow(self.avg_star, norm=norm, cmap=cmap)
        axs[0].set_title(self.titles[0])
        divider = make_axes_locatable(axs[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[0].axvline((self.avg_star.shape[0]-1)*0.5,color='black')
        axs[0].axhline((self.avg_star.shape[1]-1)*0.5,color='black')
        
        im = axs[1].imshow(self.avg_psf, norm=norm, cmap=cmap)
        axs[1].set_title(self.titles[1])
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[1].axvline((self.avg_psf.shape[0]-1)*0.5,color='black')
        axs[1].axhline((self.avg_psf.shape[1]-1)*0.5,color='black')
        
        im = axs[2].imshow(self.avg_residual, norm=norm2, cmap=cmap2)
        axs[2].set_title(self.titles[2])
        divider = make_axes_locatable(axs[2])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[2].axvline((self.avg_residual.shape[0]-1)*0.5,color='black')
        axs[2].axhline((self.avg_residual.shape[1]-1)*0.5,color='black')

        fig.tight_layout()

        plt.savefig(outname)

class mean_relative_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object, Mirage=None): #special_arg=None | if special arg is not none, make self.avg_stars middle 75 by 75 of mirage cutouts
        super().__init__(catalog_object, psf_object)
        if Mirage is not None:
            def get_middle_pixels(array_2d, n):    
                row_start = (len(array_2d) - n) // 2
                row_end = row_start + n
                col_start = (len(array_2d[0]) - n) // 2
                col_end = col_start + n
                middle_pixels = [row[col_start:col_end] for row in array_2d[row_start:row_end]]
                return middle_pixels 
            def list_of_arrays_to_3d_array(list_of_2d_arrays):
                return np.array(list_of_2d_arrays)
            print(self.cropped_stars[0].shape)
            #self.cropped_stars = [self.catalog_obj.data['MIRAGE_VIGNET'][i] for i in range(len(self.catalog_obj.data['MIRAGE_VIGNET']))]
            self.cropped_stars = [get_middle_pixels(self.catalog_obj.data['MIRAGE_VIGNET'][i], 75) for i in range(len(self.catalog_obj.data['MIRAGE_VIGNET']))]
            self.cropped_stars = list_of_arrays_to_3d_array(self.cropped_stars)
            print(self.cropped_stars[0].shape)
        #print(self.sky_level, self.sky_std)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        
        resids = []
        for i in range(len(stars)):
            relative_error = np.zeros((stars[0].shape[0], stars[0].shape[1]))
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    relative_error[j,k] = (stars[i][j,k] - psfs[i][j,k])/(stars[i][j,k] + 1.0e-6)
            resids.append(relative_error)
        for i in range(len(resids)):
            for j in range(len(resids[i])):
                for k in range(len(resids[i][j])):
                    if (resids[i][j,k] > 50) or (resids[i][j,k] < -50):
                        resids[i][j,k] = np.nan
        self.avg_residual = np.nanmean(resids, axis=0)
        self.sem = np.nanstd(self.avg_residual) #num pixels or stars?
        self.set_residuals_sum()

    def return_sem(self):
        return self.sem

class mean_absolute_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        resids = []
        for i in range(len(stars)):
            absolute_error = np.zeros((stars[0].shape[0], stars[0].shape[1]))
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    absolute_error[j,k] = np.abs(stars[i][j,k] - psfs[i][j,k])/np.abs(stars[i][j,k] + 1.0e-6)
            resids.append(absolute_error)
        for i in range(len(resids)):
            for j in range(len(resids[i])):
                for k in range(len(resids[i][j])):
                    if (resids[i][j,k] > 50) or (resids[i][j,k] < -50):
                        resids[i][j,k] = np.nan
        self.avg_residual = np.nanmean(resids, axis=0)
        self.sem = np.nanstd(self.avg_residual)
        self.set_residuals_sum()
    
    def return_sem(self):
        return self.sem
    
    def save_figure(self, outname):
        #vmin = -0.1 #np.nanmin(self.avg_psf)
        vmin = np.minimum(np.nanmin(self.avg_star), np.nanmin(self.avg_psf))
        #vmax = 1
        vmax = np.maximum(np.nanmax(self.avg_star), np.nanmax(self.avg_psf))
        
        vmin2 = np.nanmin(self.avg_residual)
        vmax2 = np.nanmax(self.avg_residual)
        if np.abs(vmin2) > np.abs(vmax2):
            vmax2 = np.abs(vmin2)
        else:
            vmin2 = -1*vmax2

        norm = colors.SymLogNorm(vmin=vmin, vmax=vmax, linthresh=1e-4)
        #norm2 = colors.SymLogNorm(vmin=0, vmax=np.nanmax(self.avg_residual), linthresh=1e-1)
        #norm2 = colors.TwoSlopeNorm(0)
        cmap=plt.cm.turbo
        #cmap=plt.cm.BrBG
        #cmap = plt.cm.inferno
        cmap2=plt.cm.Blues
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=[15,7], tight_layout=True)
        
        im = axs[0].imshow(self.avg_star, norm=norm, cmap=cmap)
        axs[0].set_title(self.titles[0])
        divider = make_axes_locatable(axs[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[0].axvline((self.avg_star.shape[0]-1)*0.5,color='black')
        axs[0].axhline((self.avg_star.shape[1]-1)*0.5,color='black')
        
        im = axs[1].imshow(self.avg_psf, norm=norm, cmap=cmap)
        axs[1].set_title(self.titles[1])
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[1].axvline((self.avg_psf.shape[0]-1)*0.5,color='black')
        axs[1].axhline((self.avg_psf.shape[1]-1)*0.5,color='black')
        
        im = axs[2].imshow(self.avg_residual, cmap=cmap2)
        axs[2].set_title(self.titles[2])
        divider = make_axes_locatable(axs[2])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[2].axvline((self.avg_residual.shape[0]-1)*0.5,color='black')
        axs[2].axhline((self.avg_residual.shape[1]-1)*0.5,color='black')

        fig.tight_layout()

        plt.savefig(outname)


class mean_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        resids = []
        for i in range(len(stars)):
            error = np.zeros((stars[0].shape[0], stars[0].shape[1]))
            for j in range(stars[0].shape[0]):
                for k in range(stars[0].shape[1]):
                    error[j,k] = (stars[i][j,k] - psfs[i][j,k]) 
                    #if psfs[i][j,k] > stars[i][j,k]:
                        #print("yabba dabba do")
            resids.append(error)

        self.avg_residual = np.nanmean(resids, axis=0)
        self.set_residuals_sum()

class chi_2_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)
        self.chi2_vals = []

    def set_residuals(self, polydim=2, polydeg=1):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        noise_map = self.cropped_err_stamps

        resids = []
        npix = stars[0].shape[0] * stars[0].shape[1]
        nparams = np.prod([polydeg+i+1 for i in range(polydim)])/polydim
        dof = npix - nparams
        chi2_vals = []

        def count_nan_inf(arr):
            # Convert the input list to a NumPy array
            np_arr = np.array(arr)
            # Count the number of NaNs and Infs in the array
            nan_count = np.count_nonzero(np.isnan(np_arr))
            inf_count = np.count_nonzero(np.isinf(np_arr))
            return nan_count, inf_count


        for i in range(len(stars)):
            if self.wg[i] == True:
                resids.append(np.divide(np.square(stars[i] - psfs[i]), np.square(noise_map[i])))
                nan_count, inf_count = count_nan_inf(resids[i])
                chi2_finite_temp = np.divide(np.square(stars[i] - psfs[i]), np.square(noise_map[i]))
                chi2_finite = chi2_finite_temp[np.isfinite(chi2_finite_temp)]
                ddof = chi2_finite.size - nparams - (nan_count + inf_count)
                #if chi2_finite.sum()/ddof > 75:
                    #print(np.sum(np.isnan(stars[i]))/chi2_finite.size)
                    #print("This star had chi square greater than 75: ", i)
                chi2_vals.append(chi2_finite.sum()/ddof)
            else:
                resids.append(np.full(stars[0].shape, np.nan))
                chi2_vals.append(np.nan)
                #chi2_vals.append(0.0)
        self.avg_residual = np.ma.mean(  np.ma.masked_invalid(resids), axis=0).data
        self.set_residuals_sum()
        self.chi2_vals = chi2_vals

    def return_chi2_vals(self):
        return self.chi2_vals
    
    def save_figure(self, outname):
        #vmin = -0.1 #np.nanmin(self.avg_psf)
        vmin = np.minimum(np.nanmin(self.avg_star), np.nanmin(self.avg_psf))
        #vmax = 1
        vmax = np.maximum(np.nanmax(self.avg_star), np.nanmax(self.avg_psf))
        
        vmin2 = np.nanmin(self.avg_residual)
        vmax2 = np.nanmax(self.avg_residual)
        if np.abs(vmin2) > np.abs(vmax2):
            vmax2 = np.abs(vmin2)
        else:
            vmin2 = -1*vmax2

        norm = colors.SymLogNorm(vmin=vmin, vmax=vmax, linthresh=1e-4)
        norm2 = colors.SymLogNorm(vmin=0, vmax=np.nanmax(self.avg_residual), linthresh=1e-4)
        #norm2 = colors.TwoSlopeNorm(0)
        cmap=plt.cm.turbo
        #cmap=plt.cm.BrBG
        #cmap = plt.cm.inferno
        cmap2=plt.cm.Blues
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=[15,7], tight_layout=True)
        
        im = axs[0].imshow(self.avg_star, norm=norm, cmap=cmap)
        axs[0].set_title(self.titles[0])
        divider = make_axes_locatable(axs[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[0].axvline((self.avg_star.shape[0]-1)*0.5,color='black')
        axs[0].axhline((self.avg_star.shape[1]-1)*0.5,color='black')
        
        im = axs[1].imshow(self.avg_psf, norm=norm, cmap=cmap)
        axs[1].set_title(self.titles[1])
        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[1].axvline((self.avg_psf.shape[0]-1)*0.5,color='black')
        axs[1].axhline((self.avg_psf.shape[1]-1)*0.5,color='black')
        
        im = axs[2].imshow(self.avg_residual, norm=norm2, cmap=cmap2)
        axs[2].set_title(self.titles[2])
        divider = make_axes_locatable(axs[2])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        axs[2].axvline((self.avg_residual.shape[0]-1)*0.5,color='black')
        axs[2].axhline((self.avg_residual.shape[1]-1)*0.5,color='black')

        fig.tight_layout()

        plt.savefig(outname)


class ssim_error_plot(resid_plot):
    def __init__(self, catalog_object, psf_object):
        super().__init__(catalog_object, psf_object)

    def set_residuals(self):
        stars = self.cropped_stars
        psfs = self.cropped_psfs
        
        ssims = []
        ssim_val = []
        nrmse = []
        
        for i in range(len(stars)):
            ssim_res = ssim(psfs[i], stars[i], full=True, win_size=3, data_range=psfs.max()-psfs.min())
            ssim_val.append(ssim_res[0])
            ssims.append(ssim_res[1] - np.median(ssim_res[1]))
            nrmse.append(normalized_root_mse(stars[i], psfs[i]))
        
        self.avg_residual = np.nanmean(ssims, axis=0)
        #self.avg_residual = np.nanmean(ssims)
        #print(f'Median SSIM: %.4f\nMedian Norm. RMSE: %.4f' % (np.median(ssim_val), np.median(nrmse)))
        self.set_residuals_sum()


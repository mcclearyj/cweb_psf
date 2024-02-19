import unittest
test = unittest.TestCase
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

from .utils import AttrDict, set_rc_params

class ResidPlots:
    '''
    Resid plots were getting messy, so made it a class!
    Ideally should define a plotter class that gets inherited by these methods
    '''

    def __init__(self, starmaker, psfmaker):
        '''
        Attributes
            stars:   instance of StarMaker
            psf:     instance of PSFMaker
            scale:   scaling for quiverplots
            *_dict:  dicts holding average image & summary statistics for Makers
                     and the various residuals plots

        test.assertIsInstance(self, starmaker, src.starmaker.StarMaker,
                        msg='starmaker must be instance of StarMaker class'
                        )
        test.assertIsInstance(self, psfmaker, src.psfmaker.PSFMaker,
                         msg='psfmaker must be instance of PSFMaker class'
                         )
        '''

        self.stars = starmaker
        self.psfs  = psfmaker

        self.star_dict = {}
        self.psf_dict = {}
        self.resid_dict = {}
        self.chi2_dict = {}


    def _make_im_dict(self, maker, stamps, wg):
        '''
        Create the dict holding average image and summary statistics
        for all the residuals plots

        Inputs
                maker:  should be an instance of either StarMaker or PSFMaker
                stamps: either list of arrays or list of Galsim Image instances
                wg:  where star and psf fits both succeeded
        Returns:
                ellip_dict: dict with e1, e2, and theta for plotting,
                            plus some summary stats, cast to a class
        '''

        fwhm = np.nanmean(maker.fwhm[wg])
        sigma = np.nanmean(maker.hsm_sig[wg])

        try:

            if type(stamps[0]) is galsim.image.Image:
                stamp_arr = []
                for stamp in stamps:
                    stamp_arr.append(stamp.array)
                avg_im = np.nanmean(stamp_arr, axis=0)
            else:
                avg_im = np.nanmean(stamps, axis=0)

        except:
            pdb.set_trace()

        im_dict = dict(avg_im = avg_im,
                        fwhm  = fwhm,
                        sigma = sigma
                        )

        return AttrDict(im_dict)


    def _populate_dicts(self):
        """
        Populate the ellipticity dictionaries for plotting
        """

        psfs  = self.psfs
        stars = self.stars
        resids = np.array(psfs.stamps - stars.stamps)/(np.array(stars.stamps) + 1.0e-6)

        wg = (psfs.hsm_g1 > -9999) & (stars.hsm_g1 > -9999)
        self.star_dict  = self._make_im_dict(stars, stars.stamps, wg)
        self.psf_dict   = self._make_im_dict(psfs, psfs.stamps, wg)
        self.resid_dict = self._make_im_dict(psfs, resids, wg)

    def make_chi2(self, polydim=2, polydeg=1, outname='chi2_residuals.png'):
        """
        Compute the chi-squared for each image. Loop over residuals to get
        chi-squared for each stamp, then create a mean (or maybe median)

        polydim: dimension of PSF fit polynomial (default=2 for X, Y fit)
        polydeg: Degree of PSF polynomial fit (default=1)
        """
        psf = self.psfs
        star = self.stars

        wg = (psf.hsm_g1 > -9999) & (star.hsm_g1 > -9999)

        npix = star.vignet_size * star.vignet_size
        nparams = np.prod([polydeg+i+1 for i in range(polydim)])/polydim
        dof = npix - nparams

        chi2_maps = []
        chi2_vals = []

        for i, psf_stamp in enumerate(np.array(psf.stamps)[wg]):
            # Make image-model chi2 map and append it to list
            noise_map = star.err_stamps[i]
            star_stamp = star.stamps[i]
            chi2_map = np.divide(np.square(star_stamp-psf_stamp),
                                 np.square(noise_map))
            chi2_maps.append(chi2_map)

            # Also compute total reduced chi2 for image-model
            chi2_finite = chi2_map[np.isfinite(chi2_map)]
            ddof = chi2_finite.size-nparams
            chi2_vals.append(chi2_finite.sum()/ddof)

        # Mask out NaNs or Infs
        masked_chi2 = np.ma.masked_invalid(chi2_maps)

        # Average (OK, mean) image
        avg_chi2_im = np.ma.mean(masked_chi2, axis=0).data

        # Total chi2
        chi_square = np.ma.sum(masked_chi2)

        # Calculate reduced chi2
        reduced_chi_square = chi_square / dof / len(chi2_maps)
        mean_reduced_chi_sq = np.mean(chi2_vals)

        # Calculate p-value
        p_value = 1 - chi2.cdf(chi_square, dof)

        # get a dict with all those values!
        chi2_dict = dict(
            avg_im = avg_chi2_im,
            reduced_chi_square = reduced_chi_square,
            mean_reduced_chi_sq = mean_reduced_chi_sq
        )
        self.chi2_dict = AttrDict(chi2_dict)

        # Save the chi image to a fits file, too
        im = fits.PrimaryHDU(np.sqrt(avg_chi2_im))

        for key in list(chi2_dict.keys())[1:]:
            im.header.set(key, chi2_dict[key])

        im.writeto(outname.replace('.png', '.fits'), overwrite=True)


    def _make_mpl_dict(self, index, vmin=None, vmax=None, avg_im=None):
        '''
        EXTREMELY SPECIFIC PLOTTING KEYWORDS
        (Plots may not make sense or errors may be thrown). Assumes that passed
        avg_im is a residual plot of some sort.
        '''
        if (avg_im is not None):
            norm = colors.TwoSlopeNorm(0, vmin=-0.2, vmax=1.1)
            cmap = plt.cm.bwr_r

        else:
            vmin = np.min(self.star_dict.avg_im)
            vmax = np.max(self.star_dict.avg_im)

            norm = colors.SymLogNorm(vmin=vmin,
                            vmax=vmax,
                            linthresh=1e-2)
            cmap=plt.cm.turbo

        mpl_dict = dict(cmap=cmap, norm=norm)

        return mpl_dict

    def _get_plot_titles(self):

        sd = self.star_dict
        pd = self.psf_dict
        rd = self.resid_dict
        xd = self.chi2_dict

        star_title = 'mean HSM $\sigma^{*} = %.4f$ pix\ngs.calculateFWHM() = %.4f$^{\prime\prime}$'\
                    % (sd.sigma, sd.fwhm)
        psf_title = 'mean HSM $\sigma^{PSF} = %.4f$ pix\ngs.calculateFWHM() = %.4f$^{\prime\prime}$'\
                    % (pd.sigma, pd.fwhm)

        # Exclude crazy outliers
        wg = (rd.avg_im.ravel() > -50) & (rd.avg_im.ravel() < 50)
        resid_title = 'mean norm. resid: %1.3f std=%1.3f\n'\
                    % (np.nanmean(rd.avg_im.ravel()[wg]),np.nanstd(rd.avg_im.ravel()[wg]))
        try:
            chi2_title = '$\overline{\chi^2_{dof}} = %.2f$\nmean $\chi^2_{dof} = %.2f$'\
                        % (xd.reduced_chi_square, xd.mean_reduced_chi_sq)
            xd.title = chi2_title
        except:
            pdb.set_trace()

        sd.title = star_title; pd.title = psf_title; rd.title = resid_title

    def _make_fig(self, dicts, mpl_dicts):
        '''
        Generic method to make residuals plots
        '''

        # First things first I'm the reallest
        set_rc_params(fontsize=16)

        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True,
                                    figsize=[15,7], tight_layout=True)
        for i, dc in enumerate(dicts):
            im = axs[i].imshow(dc.avg_im, **mpl_dicts[i])
            axs[i].set_title(dc.title)
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax)
            axs[i].axvline((dc.avg_im.shape[0]-1)*0.5,color='black')
            axs[i].axhline((dc.avg_im.shape[1]-1)*0.5,color='black')
        fig.tight_layout()
        return fig

    def make_flux_resid_plot(self, outname='flux_residuals.png'):
        """ Make flux residual image """

        dicts = [self.star_dict, self.psf_dict, self.resid_dict]
        mpl_dicts=[]

        for i, dct in enumerate(dicts):
            if i==2:
                cmap = plt.cm.bwr_r
                mpl_dict = dict(
                    norm=colors.TwoSlopeNorm(0, vmin=-5, vmax=5),
                    cmap=cmap
                )
            else:
                mpl_dict = self._make_mpl_dict(i)

            mpl_dicts.append(mpl_dict)

        # Make actual plot
        fig = self._make_fig(dicts, mpl_dicts)

        # Save plot
        fig.savefig(outname)

        # Also save a FITS file for making histograms later
        im = fits.PrimaryHDU(self.resid_dict.avg_im)
        for key in list(self.resid_dict.keys())[1:]:
            try:
                im.header.set(key, self.resid_dict[key])
            except:
                pass
        im.writeto(outname.replace('.png', '.fits'), overwrite=True)


    def make_chi2_plot(self, outname='chi2_residuals.png'):
        '''
        Make Chi-squared residual image
        '''

        dicts = [self.star_dict, self.psf_dict, self.chi2_dict]

        mpl_dicts=[]
        for i, dct in enumerate(dicts):
            if i==2:
                mpl_dict = dict(norm=colors.LogNorm(), cmap=plt.cm.RdYlBu_r)
            else:
                mpl_dict = self._make_mpl_dict(i)
            mpl_dicts.append(mpl_dict)

        # Make actual plot
        fig = self._make_fig(dicts, mpl_dicts)

        # Save it
        fig.savefig(outname)

    def make_ssim_ims(self, outname='ssim.png'):
        '''
        Calculate the structural similarity index measure (SSIM) between the
        PSF and star images. It comes from video compression, comparing
        uncompressed to compressed data.

        From Wikipedia (https://en.wikipedia.org/wiki/Structural_similarity):
        SSIM index is a decimal value between -1 and 1,
        where 1 indicates perfect similarity, 0 indicates no similarity, and
        -1 indicates perfect anti-correlation. For an image, it is typically
        calculated using a sliding Gaussian window of size 11x11 or a block
        window of size 8×8.
        The default usage is the block window. Window size must be odd.

        Also calculate  normalized root mean-squared error (NRMSE) between the
        PSF and star images. 0 is perfect correspondance, 1 is no correspondance
        '''

        psfs = self.psfs
        stars = self.stars

        wg = (psfs.hsm_g1 > -9999) & (stars.hsm_g1 > -9999)

        ssims = []
        ssim_val = []
        nrmse = []

        for i, star in enumerate(np.array(stars.stamps)[wg]):
            psf_stamp = psfs.stamps[i]
            ssim_res = ssim(psf_stamp, star, full=True, win_size=3,
                            data_range=psf_stamp.max()-psf_stamp.min())
            ssim_val.append(ssim_res[0])
            ssims.append(ssim_res[1] - np.median(ssim_res[1]))
            nrmse.append(normalized_root_mse(star, psf_stamp))

        title = f'Median SSIM: %.4f\nMedian Norm. RMSE: %.4f'\
                    % (np.median(ssim_val), np.median(nrmse))

        ssim_dict = dict(
            avg_im = np.mean(ssims, axis=0),
            ssim  = np.median(ssim_val),
            nrmse = np.median(nrmse),
            title = title
        )

        # OK ssim too, I guess
        dicts = [self.star_dict, self.psf_dict, AttrDict(ssim_dict)]

        mpl_dicts=[]

        for i, dct in enumerate(dicts):
            if i==2:
                norm = colors.CenteredNorm(0)
                cmap = plt.cm.bwr_r
                mpl_dict = dict(cmap=cmap, norm=norm)
            else:
                mpl_dict = self._make_mpl_dict(i)

            mpl_dicts.append(mpl_dict)

        # Make actual plot
        fig = self._make_fig(dicts, mpl_dicts)
        fig.savefig(outname.replace('flux', 'ssim'))

    def run(self, polydeg=1, resid_name=None, chi2_name=None):
        '''
        Make flux and residuals plots
        '''

        # Populate dicts
        self._populate_dicts()

        # Make chi-square residuals;
        self.make_chi2(polydeg=polydeg, outname=chi2_name)

        # Get titles (they're defined here!)
        self._get_plot_titles()

        # Make flux residual plots
        self.make_flux_resid_plot(resid_name)

        # And make the chi squared plot
        self.make_chi2_plot(chi2_name)

        # Bonus, do the SSIM figure
        self.make_ssim_ims(resid_name)

        return

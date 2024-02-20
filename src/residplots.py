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
    """
    Resid plots were getting messy, so made it a class!
    Ideally should define a plotter class that gets inherited by these methods
    """

    def __init__(self, starmaker, psfmaker):
        """
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
        """

        self.stars = starmaker
        self.psfs  = psfmaker

        self.star_dict = {}
        self.psf_dict = {}
        self.resid_dict = {}
        self.chi2_dict = {}


    def return_flux_resid_vals(self):
        """ Convenience function to access flux residuals """
        # Remove first key (an image)
        resid_dict = self.resid_dict
        remaining_keys = list(resid_dict.keys())[1:-1]

        # Construct a new dictionary with the remaining keys
        modified_resid_dict = {
            key: str(resid_dict[key]) for key in remaining_keys
        }
        return modified_resid_dict

    def return_ssim_vals(self):
        """ Convenience function to access SSIM values """
        # Remove first key (an image)
        ssim_dict = self.ssim_dict
        remaining_keys = list(ssim_dict.keys())[1:]

        # Construct a new dictionary with the remaining keys
        modified_ssim_dict = {
            key: ssim_dict[key] for key in remaining_keys
        }
        return modified_ssim_dict

    def _make_im_dict(self, maker, stamps, wg):
        """
        Create the dict holding average image and summary statistics
        for all the residuals plots

        Inputs
                maker:  should be an instance of either StarMaker or PSFMaker
                stamps: either list of arrays or list of Galsim Image instances
                wg:  where star and psf fits both succeeded
        Returns:
                im_dict: dict with average image, fwhm, sigma,
                         plus some summary stats, cast to a class
        """

        fwhm = np.nanmedian(maker.fwhm[wg])
        sigma = np.nanmedian(maker.hsm_sig[wg])

        try:

            if type(stamps[0]) is galsim.image.Image:
                stamp_arr = []
                for stamp in stamps:
                    stamp_arr.append(stamp.array)
                avg_im = np.nanmean(stamp_arr, axis=0)
                mean_flux = np.nanmean(stamp_arr)
                std_flux = np.nanstd(stamp_arr)
            else:
                avg_im = np.nanmean(stamps, axis=0)
                mean_flux = np.nanmean(stamps)
                std_flux = np.nanstd(stamps)

        except:
            pdb.set_trace()

        im_dict = dict(
            avg_im = avg_im,
            fwhm  = fwhm,
            sigma = sigma,
            mean_flux = mean_flux,
            std_flux = std_flux
        )

        return AttrDict(im_dict)

    def _populate_dicts(self):
        """
        Populate the ellipticity dictionaries for plotting
        """

        psfs  = self.psfs
        stars = self.stars
        resids = np.array(psfs.stamps - stars.stamps) / \
                    (np.array(stars.stamps) + 1.0e-6)

        # filter out garbage
        for i, resid in enumerate(resids):
            wb = (np.abs(resid) > 50)
            resid[wb] = np.nan
            resids[i] = resid

        wg = (psfs.hsm_g1 > -9999) & (stars.hsm_g1 > -9999)
        self.star_dict  = self._make_im_dict(stars, stars.stamps, wg)
        self.psf_dict   = self._make_im_dict(psfs, psfs.stamps, wg)
        self.resid_dict = self._make_im_dict(psfs, resids, wg)

    def _make_mpl_dict(self, index, vmin=None, vmax=None, avg_im=None):
        """
        EXTREMELY SPECIFIC PLOTTING KEYWORDS
        (Plots may not make sense or errors may be thrown). Assumes that passed
        avg_im is a residual plot of some sort.
        """
        if (avg_im is not None):
            norm = colors.TwoSlopeNorm(0, vmin=-1, vmax=1)
            cmap = plt.cm.bwr_r

        else:
            vmin = np.min(self.star_dict.avg_im)
            vmax = np.max(self.star_dict.avg_im)

            norm = colors.SymLogNorm(vmin=vmin,
                            vmax=vmax,
                            linthresh=0.05)
            cmap=plt.cm.turbo

        mpl_dict = dict(cmap=cmap, norm=norm)

        return mpl_dict

    def _get_plot_titles(self):

        sd = self.star_dict
        pd = self.psf_dict
        rd = self.resid_dict
        xd = self.chi2_dict

        star_title = \
            'median HSM $\sigma^{*} = $' + f'{sd.sigma:.4f} pix\n' +\
            f'gs.calculateFWHM() = {sd.fwhm:.4f}' + r'$^{\prime\prime}$'

        psf_title = \
            'median HSM $\sigma^{*} = $' + f'{pd.sigma:.4f} pix\n' +\
            f'gs.calculateFWHM() = {pd.fwhm:.4f}' + r'$^{\prime\prime}$'

        # Exclude crazy outliers
        wg = (rd.avg_im.ravel() > -100) & (rd.avg_im.ravel() < 100)
        resid_title = \
            f'mean norm. resid: {np.nanmean(rd.avg_im.ravel()[wg]):1.3f}' + \
            ' $\pm$' + f'{np.nanstd(rd.avg_im.ravel()[wg]): 1.3f}'

        #resid_title = \
        # 'mean norm. resid: %1.3f std=%1.3f\n' % (rd.mean_flux, rd.std_flux)

        sd.title = star_title; pd.title = psf_title; rd.title = resid_title

    def _make_fig(self, dicts, mpl_dicts):
        """
        Generic method to make residuals plots
        """

        # First things first I'm the reallest
        set_rc_params(fontsize=16)

        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True,
                                    figsize=[15, 5], tight_layout=True)

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
                    norm=colors.TwoSlopeNorm(0),
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

    def make_ssim_ims(self, outname='ssim.png'):
        """
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
        """

        psfs = self.psfs
        stars = self.stars
        psf_stamps = np.array(self.psfs.stamps)
        star_stamps = self.stars.stamps

        # Mask out NaNs or Infs
        masks = (1-np.ma.masked_invalid(psf_stamps).mask) * \
                (1-np.ma.masked_invalid(star_stamps).mask)

        all_good = np.full(len(stars.x), True)
        for i, mask in enumerate(masks):
            # This line picks out "good" stamps whose (unraveled) intersection
            # with the sentinel mask values is empty, so the length of the list is 0
            is_masked = np.size(np.intersect1d(mask, 0)) == 0
            all_good[i] *= is_masked
        # Also exclude failed HSM fits
        wg = (self.psfs.hsm_g1 > -9999) & (self.stars.hsm_g1 > -9999)
        wg *= all_good

        ssims = []
        ssim_val = []
        nrmse = []

        for i, star in enumerate(star_stamps[wg]):
            psf_stamp = psf_stamps[wg][i]
            ssim_res = ssim(psf_stamp, star, full=True, win_size=3,
                            data_range=psf_stamp.max()-psf_stamp.min())
            ssim_val.append(ssim_res[0])
            ssims.append(ssim_res[1] - np.nanmedian(ssim_res[1]))
            nrmse.append(normalized_root_mse(star, psf_stamp))

        title = f'Median SSIM: %.3f\nMedian Norm. RMSE: %.3f'\
                    % (np.median(ssim_val), np.median(nrmse))

        ssim_dict = {
            'avg_im': np.nanmedian(ssims, axis=0),
            'ssim' : np.nanmean(ssim_val),
            'median_ssim': np.nanmedian(ssim_val),
            'std_ssim': np.std(ssim_val),
            'nrmse': np.nanmean(nrmse),
            'std_nrmse': np.std(nrmse),
            'title': title
        }

        self.ssim_dict = AttrDict(ssim_dict)

        print('')
        print(f'For PSF type {psfs.psf_type.upper()}:')
        print(
            f'\tmean SSIM = {self.ssim_dict.ssim:.3f}' + \
            f' +/- {self.ssim_dict.std_ssim:.3f}'
        )
        print(
            f'\tmean NRMSE = {self.ssim_dict.nrmse:.3f} ' + \
            f'+/- {self.ssim_dict.std_nrmse:.3f}'
        )

        # OK ssim too, I guess
        dicts = [self.star_dict, self.psf_dict, self.ssim_dict]

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

    def run(self, polydeg=1, resid_name=None):
        """
        Make flux and residuals plots
        """

        # Populate dicts
        self._populate_dicts()

        # Get titles (they're defined here!)
        self._get_plot_titles()

        # Make flux residual plots
        self.make_flux_resid_plot(resid_name)

        # Bonus, do the SSIM figure
        self.make_ssim_ims(resid_name)

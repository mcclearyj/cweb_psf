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

class ChiSqPlots:
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
        self.chi2_dict = {}

    def return_chi2_resid_vals(self):
        """ Convenience function to access chisq residuals """
        # Wanna remove the avg im before returning
        chi2_dict = self.chi2_dict
        remaining_keys = list(chi2_dict.keys())[1:]

        # Construct a new dictionary with the remaining keys
        modified_chi2_dict = {
            key: chi2_dict[key] for key in remaining_keys
        }
        return modified_chi2_dict

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

    def make_chi2(self, polydim=2, polydeg=1, save_fits=False):
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

        # Make image-model chi2 map and append it to list
        for i, psf_stamp in enumerate(np.array(psf.stamps)[wg]):
            noise_map = star.err_stamps[i]
            star_stamp = star.stamps[i]
            chi2_map = np.divide(
                np.square(star_stamp-psf_stamp),
                np.square(noise_map)
            )
            chi2_maps.append(chi2_map)

            # Also compute total reduced chi2 for image-model
            chi2_finite = chi2_map[np.isfinite(chi2_map)]
            ddof = chi2_finite.size-nparams
            chi2_vals.append(chi2_finite.sum()/ddof)

        # Mask out NaNs or Infs
        masked_chi2 = np.ma.masked_invalid(chi2_maps)

        # Average (OK, mean) image
        avg_chi2_im = np.ma.median(masked_chi2, axis=0).data

        # Total chi2
        chi_square = np.ma.sum(masked_chi2)

        # Calculate reduced chi2
        reduced_chi_square = chi_square / dof / len(chi2_maps)

        # Also descriptive stats
        mean_reduced_chi_sq = np.median(chi2_vals)
        std_reduced_chi_sq = np.std(chi2_vals)

        # Calculate p-value
        p_value = 1 - chi2.cdf(chi_square, dof)

        # get a dict with all those values!
        chi2_dict = dict(
            avg_im = avg_chi2_im,
            reduced_chi_square = reduced_chi_square,
            mean_reduced_chi_sq = mean_reduced_chi_sq,
            std_reduced_chi_sq = std_reduced_chi_sq,
            p_value = p_value
        )
        self.chi2_dict = AttrDict(chi2_dict)

        # Save the chi image to a fits file, too
        if save_fits == True:
            im = fits.PrimaryHDU(np.sqrt(avg_chi2_im))
            for key in list(chi2_dict.keys())[1:]:
                im.header.set(key, chi2_dict[key])
            im.writeto(outname.replace('.png', '.fits'), overwrite=True)

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
        xd = self.chi2_dict

        star_title = \
            'median HSM $\sigma^{*} = $' + f'{sd.sigma:.4f} pix\n' + \
            f'gs.calculateFWHM() = {sd.fwhm:.4f}' + r'$^{\prime\prime}$'

        psf_title = \
            'median HSM $\sigma^{*} = $' + f'{pd.sigma:.4f} pix\n' + \
            f'gs.calculateFWHM() = {pd.fwhm:.4f}' + r'$^{\prime\prime}$'

        chi2_title = \
            'reduced $\chi^2_{dof} = $' + f'{xd.reduced_chi_square:.2f}' + \
            '\n$\overline{\chi^2_{dof}} = $' + \
            f'{xd.mean_reduced_chi_sq:.2f}' + \
            ' $\pm$ ' + f'{xd.std_reduced_chi_sq:.2f}'

        sd.title = star_title; pd.title = psf_title; xd.title = chi2_title

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

    def make_chi2_plot(self, outname='chi2_residuals.png'):
        """
        Make Chi-squared residual image
        """

        # Get list of star, psf, and chi2 dicts
        dicts = [self.star_dict, self.psf_dict, self.chi2_dict]

        # Will hold matplotlib plot parameters
        mpl_dicts=[]

        for i, dct in enumerate(dicts):
            if i==2:
                mpl_dict = dict(
                    norm = colors.LogNorm(),
                    cmap=plt.cm.RdYlBu_r
                )
            else:
                mpl_dict = self._make_mpl_dict(i)
            mpl_dicts.append(mpl_dict)

        # Make actual plot
        fig = self._make_fig(dicts, mpl_dicts)

        # Save it
        fig.savefig(outname)

    def run(self, chisq_name=None, polydeg=1):
        """
        Make flux and residuals plots
        """

        # Populate dicts
        self._populate_dicts()

        # Make chi-square residuals; this stays here!
        self.make_chi2(polydeg=polydeg)

        # Get titles (they're defined here!)
        self._get_plot_titles()

        # And make the chi squared plot
        self.make_chi2_plot(chisq_name)

        return

import unittest
test = unittest.TestCase
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import diagnostics
import ipdb,pdb
from astropy.io import fits
from scipy.stats import chi2


from .utils import AttrDict, set_rc_params

class ResidPlots:
    '''
    Resid plots were getting messy, so made it a class!
    Ideally should define a plotter class that gets inherited by
    '''

    def __init__(self, starmaker, psfmaker):
        '''
        Attributes
            stars:  instance of StarMaker
            psf:  instance of PSFMaker
            scale:  scaling for quiverplots
        '''
        test.assertIsInstance(self, starmaker, diagnostics.starmaker.StarMaker,
                        msg='starmaker must be instance of StarMaker class'
                        )
        test.assertIsInstance(self, psfmaker, diagnostics.psfmaker.PSFMaker,
                         msg='psfmaker must be instance of PSFMaker class'
                         )

        self.stars = starmaker
        self.psfs  = psfmaker

        self.star_dict = {}
        self.psf_dict = {}
        self.resid_dict = {}
        self.chi2_dict = {}


    def _make_im_dict(self, maker, stamps, wg):
        '''
        Calculate ellipticity parameters using input HSM moments

        Inputs
                maker:  should be an instance of either StarMaker or PSFMaker
                wg:  where star and psf fits both succeeded
        Returns:
                ellip_dict: dict with e1, e2, and theta for plotting,
                            plus some summary stats, cast to a class
        '''

        fwhm = np.nanmean(maker.fwhm[wg])
        sigma = np.nanmean(maker.hsm_sig[wg])
        avg_im = np.nanmean(stamps, axis=0)

        im_dict = dict(avg_im = avg_im,
                        fwhm  = fwhm,
                        sigma = sigma
                        )

        return AttrDict(im_dict)


    def _populate_dicts(self):
        '''
        Populate the ellipticity dictionaries for plotting
        '''

        psfs  = self.psfs
        stars = self.stars

        wg = (psfs.hsm_g1 > -9999) & (stars.hsm_g1 > -9999)
        self.star_dict  = self._make_im_dict(stars, stars.stamps, wg)
        self.psf_dict   = self._make_im_dict(psfs, psfs.stamps, wg)
        self.resid_dict = self._make_im_dict(psfs, psfs.resids, wg)

        return


    def make_chi2(self, nparams=3, outname='chi2_residuals.png'):
        '''
        Compute the chi-squared for each image. Loop over residuals to get
        chi-squared for each stamp, then create a mean (or maybe median)
        '''
        psf = self.psfs
        star = self.stars
        observed = self.star_dict.avg_im # for DOF

        chi2_maps = []
        for i, resid in enumerate(psf.resids):
            noise_map = np.full(resid.shape, np.std(star.stamps[i]))
            chi2_map = np.square(np.divide(resid, noise_map))
            chi2_maps.append(chi2_map)

        # Total chi2: should it be chip by chip or computed on avg?
        avg_chi2_im = np.nanmean(chi2_maps, axis=0)
        chi_square = np.sum(avg_chi2_im)


        # Calculate degrees of freedom
        # Q: are there 2 b/c of X, Y or is it 3: X, Y, flux?
        dof = observed.size - nparams

        # Calculate reduced chi2
        reduced_chi_square = chi_square / dof

        # Calculate p-value
        p_value = 1 - chi2.cdf(chi_square, dof)

        # get a dict with all those values!
        chi2_dict = dict(avg_im = avg_chi2_im,
                            chi_square = chi_square,
                            reduced_chi_square = reduced_chi_square,
                            p_value = p_value
                            )
        self.chi2_dict = AttrDict(chi2_dict)

        # Save the chi-squared image to a fits file, too
        im = fits.PrimaryHDU(avg_chi2_im)
        for key in list(chi2_dict.keys())[1:]:
            im.header.set(key, chi2_dict[key])
        im.writeto(outname.replace('.png', '.fits'), overwrite=True)

        return


    def _make_mpl_dict(self, index, vmin=None, vmax=None, avg_im=None):
        '''
        EXTREMELY SPECIFIC PLOTTING KEYWORDS -- CHANGE WITH CAUTION
        (Plots may not make sense or errors may be thrown). Assumes that passed
        avg_im is a residual plot of some sort.

            #norm = colors.TwoSlopeNorm(np.median(avg_im),
            #        vmin=0.8*np.min(avg_im),
            #        vmax=0.8*np.max(avg_im))

        I used to use colors.TwoSlopeNorm and the seismic_r color map for the
        flux residuals, but have decided to go with SymLogNorm for now.
        '''

        if (avg_im is not None):

            norm = colors.SymLogNorm(linthresh=0.01,
                            vmin=0.8*np.min(avg_im),
                            vmax=0.8*np.max(avg_im))

        else:
            if vmin == None:
                vmin = 0.0001
            if vmax == None:
             vmax = 75
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)

        cmap = plt.cm.bwr_r

        mpl_dict = dict(cmap=cmap, norm=norm)

        return mpl_dict


    def _get_plot_titles(self):

        sd = self.star_dict
        pd = self.psf_dict
        rd = self.resid_dict
        xd = self.chi2_dict

        star_title = 'avg star HSM sigma = %.4f\ngs.calculateFWHM() = %.4f'\
                    % (sd.sigma,sd.fwhm)
        psf_title = 'avg PSF HSM sigma = %.4f\ngs.calculateFWHM() = %.4f'\
                    % (pd.sigma,pd.fwhm)
        resid_title = 'sum(mean resid)= %.3f\nmean=%1.2e std=%.3f'\
                    % (np.nansum(rd.avg_im), np.nanmean(rd.avg_im), np.nanstd(rd.avg_im))
        chi2_title = '$\chi^2 = %.3f\ \chi^2_{dof} = %.3f$\np-value = %.3f'\
                    % (xd.chi_square, xd.reduced_chi_square, xd.p_value)


        sd.title = star_title; pd.title = psf_title
        rd.title = resid_title; xd.title = chi2_title

        return


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
            axs[i].axvline((dc.avg_im.shape[0]-1)*0.5,color='black')
            axs[i].axhline((dc.avg_im.shape[1]-1)*0.5,color='black')
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax)

        return fig


    def make_flux_resid_plot(self, outname='flux_residuals.png'):
        '''
        Make Flux residual image
        '''

        dicts = [self.star_dict, self.psf_dict, self.resid_dict]

        # Make plot param dicts
        mpl_dicts=[]
        for i, dct in enumerate(dicts):
            if i==2:
                mpl_dict = self._make_mpl_dict(index=i, avg_im=dct.avg_im)
            else:
                mpl_dict = self._make_mpl_dict(index=i)
            mpl_dicts.append(mpl_dict)

        # Make actual plot
        fig = self._make_fig(dicts, mpl_dicts)

        # Save plot
        fig.savefig(outname)


    def make_chi2_plot(self, outname='chi2_residuals.png'):
        '''
        Make Chi-squared residual image
        '''

        dicts = [self.star_dict, self.psf_dict, self.chi2_dict]

        #pdb.set_trace()

        mpl_dicts=[]
        for i, dct in enumerate(dicts):
            star_norm = colors.LogNorm(vmin=np.min(self.star_dict.avg_im),
                                vmax=np.max(self.star_dict.avg_im))
            if i==2:
                mpl_dict = dict(norm=colors.LogNorm(), cmap=plt.cm.gist_ncar)
            else:
                mpl_dict = dict(norm=star_norm, cmap=plt.cm.turbo)
            mpl_dicts.append(mpl_dict)

        # Make actual plot
        fig = self._make_fig(dicts, mpl_dicts)

        # Save it
        fig.savefig(outname)


    def run(self, resid_name=None, chi2_name=None):
        '''
        Make flux and residuals plots
        '''

        # Populate dicts
        self._populate_dicts()

        # Make chi-square residuals;
        self.make_chi2(nparams=3, outname=chi2_name)

        # Get titles (they're defined here!)
        self._get_plot_titles()

        # Make flux residual plots
        self.make_flux_resid_plot(resid_name)

        # And make the chi squared plot
        self.make_chi2_plot(chi2_name)

        return

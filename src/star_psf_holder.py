import psfex
import galsim
import treecorr
import piff
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from astropy.table import Table, vstack
import pdb
from argparse import ArgumentParser
import yaml
import fitsio
from src.plotter import compare_rho_stats
from src.hsm_fitter import do_hsm_fit

###
### Guess I had the right idea the first time!
###

class StarPSFHolder:
    """
    Object to hold position, vignettes, and information for different PSF types.
    As currently envisioned, one of these is instantiated per PSF type and/or star.
    Parameters
        psf_type: what PSF type is being stored (or star)
        vignet_size: box size
        pixel_scale: detector pixel scale (default 0.03 as/pix)
        rho_params: parameters for rho statistics (convenience thing)
        stamps: actual renderings of the stars or PSFs (aka vignettes)
        x/y: positions on image where PSF is being evaluated
        hsm_{sig/g1/g2}: HSM adaptive moment fit parameters
        fwhm: PSF or star FWHM
    """

    def __init__(self, psf_type, vignet_size=None,
                 pixel_scale=None, rho_params=None, vb=False):
        self.psf_type = psf_type
        self.vignet_size = vignet_size
        self.pixel_scale = pixel_scale
        self.rho_params = rho_params
        self.vb = vb
        self.stamps = []
        self.x = []
        self.y = []
        self.hsm_sig = []
        self.hsm_g1 = []
        self.hsm_g2 = []
        self.fwhm = []

        self._load_defaults(rho_params, pixel_scale)

    def _load_defaults(self, rho_params, pixel_scale):
        if rho_params == None:
            self.rho_params={
                'min_sep':200,'max_sep':10000,'nbins':11
            }
            if self.vb == True:
                print(f'Warning: no rho params supplied')
                print(f'Setting rho_params to default {self.rho_params}')
        if pixel_scale == None:
            self.pixel_scale = 0.03
            if self.vb == True:
                print(f'Warning: no pixel scale supplied')
                print(f'Setting pixel scale to default {self.pixel_scale}')

    def copy_from_holder(self, holder):
        """
        Convenience function to set attributes of StarPSFHolder from a
        pre-existing StarPSFHolder
        """
        self.x = holder.x
        self.y = holder.y
        self.sky_level = holder.sky_level
        self.sky_std = holder.sky_std
        self.star_fluxes = holder.star_fluxes
        self.vignet_size = holder.vignet_size
        self.pixel_scale = holder.pixel_scale

    def load_psf_stamps(self, catalog):
        """ Load PSF stamps. Catalog should be astropy Table """
        psf_key = f'{self.psf_type.upper()}_VIGNET'
        if psf_key in catalog.dtype.names:
            psfs = catalog[psf_key]
            self.stamps = psfs
        else:
            raise KeyError(f'No PSF of type {psf_type.upper()}_VIGNET found')

    def _calc_sky_bkg(self, im_stamp):
        """
        This will sample the corners to get median/std dev of sky background
        in the supplied image stamp ('im_stamp'). Returns median and std dev
        """
        # Take outer 10% of stamps or 4x4 pixel box, whichever is larger
        j = min(im_stamp.shape[0]//10, 10)

        # Fill an array with these corners, take median and std to get an
        # sky background and variance for the stamp.
        im_stamp[im_stamp<= -999] = np.nan
        substamps = []
        substamps.append(im_stamp[-j:, -j:])
        substamps.append(im_stamp[0:j, 0:j])
        substamps.append(im_stamp[0:j, -j:])
        substamps.append(im_stamp[-j:, 0:j])
        substamps = np.array(substamps)

        # Save as attribute of class (?)
        sky_med = np.nanmean(substamps)
        sky_std = np.nanstd(substamps)

        return sky_med, sky_std

    def calc_sky_bkg(self):
        """ Call _calc_bkg for a (list of) (star) stamps """
        sky_bkg_stats = np.array(list(map(self._calc_sky_bkg, self.stamps)))
        self.sky_level = sky_bkg_stats[:, 0]
        self.sky_std = sky_bkg_stats[:, 1]

    def _trim_stamps(self, stamp):
        # What shape is the stamp?
        stamp_dim = stamp.shape[0]

        # Does stamp need trimming?
        if stamp_dim > self.vignet_size:
            ns = (stamp_dim - self.vignet_size) // 2
            ks = ns; js = -ns
            if self.vb == True:
                print(f'Trimming {ns} pixels from {self.psf_type} stamp')
        else:
            ks = 0; js = None

        # Trim (or don't) and return
        return stamp[ks:js, ks:js]

    def trim_stamps(self, star_holder):
        """
        Trim stamps to be vignet size; in the absence of a supplied vignet
        size argument, just take the smaller of the PSF and star image.

        Input
            star_holder: should be a StarPSFHolder instance
        Returns
            trimmed star stamps too

        TO DO: allow for a PSF dim < star dim!
        """

        stars = star_holder.stamps; psfs = self.stamps
        star_dim = stars[0].shape[0]; psf_dim = psfs[0].shape[0]

        if self.vignet_size == None:
            self.vignet_size = min(star_dim, psf_dim)

        # Trim stars and PSF images if needed
        trimmed_stars = list(
            map(self._trim_stamps, stars)
        )
        trimmed_psfs = list(
            map(self._trim_stamps, psfs)
        )

        # Replace objects with their trimmed versions
        star_holder.stamps = np.array(trimmed_stars)
        self.stamps = np.array(trimmed_psfs)

    def _add_flux(self, psf_stamp, star_flux,
                    sky_level, sky_std):
        """
        Add appropriate background noise and star flux to PSF renderings.
        There may be a nicer way to do this that ensures that the total flux
        in the PSF stamp matches the total star flux.
        """
        psf_stamp /= np.sum(psf_stamp)
        psf_stamp *= star_flux

        # Add noise; empirically, sky level is better than sex bkg
        # for large vignettes. Well, maybe.
        if self.psf_type not in ['piff']:
            if self.vb == True: print("\nAdding noise to PSF\n")
            noise = np.random.normal(
                loc=sky_level,
                scale=sky_std,
                size=psf_stamp.shape
            )
            psf_stamp += noise

        return psf_stamp

    def add_flux(self):
        """ Add flux and bkg noise to PSF stamps """
        psf_ims = self.stamps
        psf_flux_images = list(map(
            self._add_flux, psf_ims, self.star_fluxes,
            self.sky_level, self.sky_std
            )
        )
        self.stamps = psf_flux_images

    def do_hsm_fit(self):
        """ Convenience function for do_hsm_fit code """
        do_hsm_fit(maker=self, verbose=self.vb)

    def run_stars(self):
        """ Routines to run to load up stars """
        # Calculate sky background
        self.calc_sky_bkg()

        # Do HSM fitting
        self.do_hsm_fit()

    def run_psfs(self, stars):
        """ Maybe do something like this?
        Stars should be instance of StarPSFHolder"""

        # Load PSF stamps of type self.psf_type
        self.load_psf()

        # Trim star and PSF stamps
        stars = self.trim_stamps(stars)

        # Add star flux and background noise to PSF stamp
        self.add_flux()

        pass

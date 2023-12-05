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
from astropy.io import fits

from src.plotter import compare_rho_stats


def concatenate_catalogs(catnames, outname=None):
    """This I want to run on all the teeny single exposure catalogs.
    I am going to operate like I have the run config available, too"""

    if outname=None:
        outname='combined_catalog.fits'

    fits_files = glob.glob(catnames)

    # List to hold tables
    table_list = []

    # Loop over all your FITS files
    for fits_file in fits_files:
        table = Table.read(fits_file)
        table_list.append(table)

    # Vertically stack tables
    combined_table = vstack(table_list)

    # Save the combined table to a new FITS file
    combined_table.write(outname, format='fits')

    return combined_table

class PSFDiagnostics:
    """Calculate star stamp sky bkg stats, add to PSF images, make figures"""

    def __init__(self, catalog, config, psf_type=None, vignet_size=None):
        self.config = config # Should already have been read in
        self.psf_type = psf_type
        self.vignet_size = vignet_size
        self.stars = []
        self.psfs = []
        self.sky_level = []
        self.sky_std = []

        # open if file, else set catalog to catalog
        if type(catalog) == str:
            self.cat_fits = fits.open(catalog, 'rw')
            hdu = config['input_catalog']['hdu']
            self.catalog = self.cat_fits[hdu].read()
        else:
            self.catalog = catalog

        # load stamps, also trim to vignet_size
        self._load_stamps()


    def _load_stamps(self):
        """ Load star and PSF stamps, trimming to VIGNET size if needed """

        # Stars are easy
        stars = self.catalog['VIGNET']

        # Load PSFs too
        psf_key = f'{self.psf_type.upper()}_VIGNET'
        if psf_key in self.catalog.colnames:
            psfs = self.catalog[f'{psf_type.upper()}_VIGNET']
        else:
            raise KeyError(f'No PSF of type {psf_type.upper()}_VIGNET} found')

        # This should be a list/map thing, which would mean it makes MORE
        # sense to have that boxcutter trim boxes method be generic
        self.trim_stamps(stars, psfs)


    def trim_stamps(self, stars, psfs):
        """ Trim stamps to be vignet size; in the absence of a supplied vignet
        size argument, just take the smaller of the PSF and star image """

        star_dim = stars[0].shape[0]; psf_dim = psfs[0].shape[0]

        if self.vignet_size == None:
            self.vignet_size = np.min(star_dim, psf_dim)

        # Trim stars?
        if star_dim > self.vignet_size:
            ns = (star_dim-self.vignet_size) // 2
            ks = ns; js = -ns
        else:
            ks = 0; js = None

        # Trim PSFs?
        if psf_dim > self.vignet_size:
            np = (psf_dim-self.vignet_size) // 2
            kp = np; jp = -np
        else:
            kp = 0; jp = None

        for star_im, psf_im


    def _calc_bkg(self, im_stamp):
        """This will sample the corners to get median/std dev of sky background in
        the supplied image stamp ('im_stamp'). Returns median and std dev"""

        # Take outer 10% of stamps or 4x4 pixel box, whichever is larger
        j = max(im_stamp.shape[0]//10, 5)

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

    def calc_bkg(self, psf_vignets):
        """Call _calc_bkg for a (list of) (star) stamps"""

        sky_bkg_stats = np.array(list(map(_calc_bkg,self.psfs)))
        self.sky_level = sky_bkg_stats[:, 0]
        self.sky_std = sky_bkg_stats[:, 1]

    def _add_flux(self, psf_stamp, star_flux, sky_level, sky_std):
        """Placeholder for a method to add flux!
        If I adopt some kind of a Plotter class, as Eddie does, I may instead
        make all of these attributes.
        """

        # There may be a nicer way to do this to ensure that the total flux
        # in the PSF stamp matches the total star flux. Perhaps multiply by the
        # sum of the stamp? Don't want to add background pre-emptively.
        psf_stamp *= star_flux

        # Add noise; empirically, sky level is better than sex bkg
        # for large vignets. Maybe.
        noise = np.random.normal(loc=sky_level,
                                    scale=sky_std,
                                    size=psf_stamp.shape
                                    )
        psf_stamp+=noise

        return psf_stamp

    def add_flux(self):
        """ Add flux to PSF stamps """
        psf_flux_images = list(map(self._add_flux,
                          self.sky_level, self.sky_std)
                          )
        psf_flux_images = list(map(_add_flux, psf_vignets, sky_level, sky_std)
        self.psfs = psf_flux_images


def run_all(self, stars, vb=False, outdir='./psf_diagnostics'):
    '''
    stars is expected to be an instance of the StarMaker class
    Possible improvements: allow user to supply just X,Y?
    Allow a freestanding bg value?

    Note: it's shady to pass stars in run_all but w/e
    '''

    self.stars = stars
    self.sky_level = stars.sky_level
    self.sky_std = stars.sky_std
    self.x = stars.x
    self.y = stars.y

    if self.vignet_size == None:
        self.vignet_size = stars.vignet_size
        print(f'PSFmaker vignet size is {self.vignet_size}')

    # Render PSF, take residual against stars
    for i in range(len(stars.x)):
        xpos = stars.x[i]; ypos = stars.y[i]
        star_stamp = stars.stamps[i]
        psf_model = self.render_psf(x=xpos, y=ypos, index=i)
        if type(psf_model) is galsim.image.Image:
            psf_stamp = psf_model.array
        else:
            psf_stamp = psf_model
        try:
            self.models.append(psf_model)
            self.stamps.append(psf_stamp)
            self.resids.append(star_stamp - psf_stamp)
        except:
            pdb.set_trace()

    # Do HSM fitting
    do_hsm_fit(maker=self, verbose=vb)

    # Make output quiverplot
    quiv_name = os.path.join(outdir, '_'.join([self.psf_type,'quiverplot.png']))
    quiverplot = QuiverPlot(starmaker=stars, psfmaker=self)
    quiverplot.run(scale=1, outname=quiv_name)

    # Make output star-psf residuals plot
    resid_name = os.path.join(outdir,'_'.join([self.psf_type,'flux_resid.png']))
    chi2_name = os.path.join(outdir,'_'.join([self.psf_type,'chi2.png']))
    resid_plot = ResidPlots(starmaker=stars, psfmaker=self)
    resid_plot.run(resid_name=resid_name, chi2_name=chi2_name)

    # Compute & make output rho-statistics figures
    rho_params = self.rho_params
    if rho_params == None:
        rho_params={'min_sep':200,'max_sep':5000,'nbins':10}

    self.run_rho_stats(stars=stars,rho_params=rho_params,vb=vb,outdir=outdir)

    print("finished running PSFMaker()")
    return

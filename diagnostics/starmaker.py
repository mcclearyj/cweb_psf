import psfex
import galsim, galsim.des
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from astropy.table import Table
import pdb

from diagnostics.hsm_fitter import do_hsm_fit

class StarMaker():
    '''
    Class to store catalog stars and fit information.
     - Read in star entry from some catalog
     - Trim as needed, make a GS object. Get HSM fit.
     - This object will store all star vignets
    '''

    def __init__(self,cat_stars=None, bg_obj=None, pix_scale=0.03):
        '''
        cat_stars is either SExtractor catalog,
        catalog file path, or simply list of np.ndarrays
        '''

        self.cat_stars = cat_stars
        self.pixel_scale = pix_scale
        self.vignet_size = 21
        self.sky_level = 0.0
        self.sky_std = 0.0

        self.x = []
        self.y = []
        self.stamps = []
        self.star_flux = []

        self.hsm_sig = []
        self.hsm_g1 = []
        self.hsm_g2 = []
        self.fwhm = []


    def _read_cat(self,vb=False):
        '''
        This could be expanded to do all the star catalog reading,
        but keep it simple for now & give it a pre-read star catalog
        '''

        #self.cat_stars = cat_stars

        if vb==True:
            print("fitting to %d stars"  % len(self.cat_stars))

        return


    def _set_background(self,bg_obj=None,vb=False):
        '''
        bkg_obj is expected to be an instance of the
        StampBackground class
        '''
        if bg_obj is not None:
            self.sky_level = bg_obj.sky_level
            self.sky_std = bg_obj.sky_std

        else:
            self.sky_level = 0.0
            self.sky_std = 0.0

        if vb==True:
            print("sky level = %.3f +/- %.3f" % (self.sky_level,self.sky_std))

        return


    def _get_star_vignets(self):
        '''
        Make star stamps from SExtractor catalog vignets
        '''

        n = np.floor(0.5*(np.shape(self.cat_stars['VIGNET'])[1]-self.vignet_size))
        n = int(n)

        for i in range(len(self.cat_stars)):
            this_vign = self.cat_stars[i]['VIGNET']
            x_pos = self.cat_stars[i]['X_IMAGE']; y_pos = self.cat_stars[i]['Y_IMAGE']

            this_vign[this_vign <= -999] = np.nan
            this_vign[np.isnan(this_vign)] = self.sky_level
            vign_cutout = this_vign[n:-n,n:-n]
            star_flux = np.nansum(vign_cutout) - np.size(vign_cutout)*self.sky_level

            self.x.append(x_pos); self.y.append(y_pos)
            self.star_flux.append(star_flux)
            self.stamps.append(vign_cutout)

        self.x=np.array(self.x)
        self.y=np.array(self.y)

        return


    def run(self,vb=False,bg_obj=None):

        # Kind of does nothing right now but whatever
        self._read_cat(vb=vb)

        # Set backgrounds using supplied bkg_obj
        self._set_background(bg_obj=bg_obj,vb=vb)

        # Create star vignets
        self._get_star_vignets()

        # Create GS Object & record fit
        do_hsm_fit(self)

        return


class StampBackground():
    '''
    Determining and storing star cutout backgrounds for shape fitting purposes
    Maybe not efficient, but run the first time on stars. If not needed, don't

    : sky_level : median background of star stamp
    : sky_std   : standard deviation of star stamp
    : cat       : filepath, astropy.Table() instance or list of arrays
    : vclip     : side of sub-stamp to sample from later stamp in self.cat

    '''

    def __init__(self,cat=None,sky_level=None,sky_std=None):

        self.cat = cat
        self.vclip = 6
        self.sky_level = 0.0
        self.sky_std = 0.0
        self.substamps = []

        if sky_level is not None:
            print("setting sky level to %d" % sky_level)
            self.sky_level = sky_level
        if sky_std is not None:
            print("setting sky stdev to %d" % sky_std)
            self.sky_std = sky_std

    def _read_star_cat(self):
        pass

        return


    def calc_star_bkg(self, vb=False):
        '''
        Reading in file if needed, compute sky background of either
        SEXtractor VIGNETS or just arrays

        Is there better way to write this than a series of if statements?
        '''

        if self.cat is None:
            print("Catalog can't be 'None' if calculating sky background!")
            print("Please supply `cat` parameter")
            sys.exit()

        if type(self.cat) == 'str':
            obj_cat = Table.read(self.cat)
            cutouts = obj_cat['VIGNET']
        elif type(self.cat) == Table:
            cutouts=self.cat['VIGNET']
        elif (type(self.cat) == list) and (type(self.cat[0]) == np.ndarray):
            cutouts = self.cat
        else:
            # Go with God
            cutouts = self.cats

        self.sky_level, self.sky_std = self._calc_stamp_bkg(cutouts)

        if vb==True:
            print('star bkg = %.3f +/- %.3f' % (self.sky_level, self.sky_std))
        return self.sky_level, self.sky_std


    def _calc_stamp_bkg(self,cutouts):
        '''
        Input:
            cutouts : a list of np.ndarrays representing star images
        Output:
            sky_level : the median sky level in the star cutouts
            sky_std : standard deviation of sky level in star cutouts

        '''

        j = self.vclip

        for cutout in cutouts:
            cutout[cutout<= -999] = np.nan
            self.substamps.append(cutout[-j:,-j:])
            self.substamps.append(cutout[0:j,0:j])

        sky_level = np.nanmedian(self.substamps)
        sky_std = np.nanstd(self.substamps)

        return sky_level, sky_std

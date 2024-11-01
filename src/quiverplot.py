import unittest
test = unittest.TestCase
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import ipdb, pdb

from .utils import AttrDict, set_rc_params

class QuiverPlot:
    """
    Quiverplot plots were getting messy, so made it a class!
    """

    def __init__(self, starmaker, psfmaker):
        """
        Attributes
            stars:  instance of StarMaker
            psf:  instance of PSFMaker
            scale:  scaling for quiverplots
        """

        self.stars = starmaker
        self.psfs  = psfmaker

        self.scale = 0
        self.star_dict  = {}
        self.psf_dict  = {}
        self.resid_dict = {}

    def return_hsm_resid_vals(self):
        """ Convenience function to access HSM residuals """
        return_dict = {
            'median_g': self.resid_dict.median_g,
            'median_e': self.resid_dict.median_e,
            'median_sigma': self.resid_dict.median_sigma,
            'mean_sigma': self.resid_dict.mean_sigma,
            'std_sigma': self.resid_dict.std_sigma
        }

        return return_dict

    def _make_ellip_dict(self, maker, wg):
        """
        Calculate ellipticity parameters using input HSM moments

        Inputs
                maker:  should be an instance of either StarMaker or PSFMaker
                wg:  where star and psf fits both succeeded
        Returns:
                ellip_dict: dict with e1, e2, and theta for plotting,
                            plus some summary stats, cast to a class
        """

        g1 = maker.hsm_g1[wg]
        g2 = maker.hsm_g2[wg]
        sigma = maker.hsm_sig[wg] * maker.pixel_scale

        g = np.sqrt(g1**2 + g2**2)
        theta = 0.5*np.arctan2(g2, g1)
        e1 = g * np.cos(theta)
        e2 = g * np.sin(theta)
        e = np.sqrt(e1**2 + e2**2)
        # Summary stats for plot titles
        median_g = np.median(g)
        median_sigma = np.median(sigma)
        median_e = np.median(e)

        ellip_dict = {
            'g1': g1,
            'g2': g2,
            'e1': e1,
            'e2': e2,
            'e': e,
            'theta': theta,
            'sigma': sigma,
            'median_g': median_g,
            'median_sigma': median_sigma,
            'median_e': median_e
        }

        return AttrDict(ellip_dict)


    def _make_resid_dict(self, maker1, maker2):
        """
        Calculate ellipticity parameters for star-psf model residuals
        Residual quantities follow the PIFF convention.

        Inputs:
                maker1: should be an instance of either StarMaker or PSFMaker
                maker2: should be an instance of either StarMaker or PSFMaker
        Returns:
                resid_dict: dict with e1, e2, and theta for plotting,
                            plus some summary stats
        """

        resid_g1 = maker1.g1 - maker2.g1
        resid_g2 = maker1.g2 - maker2.g2
        resid_sigma = maker1.sigma - maker2.sigma

        resid_g  = np.sqrt(resid_g1**2 + resid_g2**2)
        theta_resid = 0.5 * np.arctan2(resid_g2, resid_g1)
        resid_e1 = resid_g * np.cos(theta_resid)
        resid_e2 = resid_g * np.sin(theta_resid)
        resid_e = np.sqrt(resid_e1**2 + resid_e2**2)

        median_resid_g = np.median(resid_g)
        median_resid_e = np.median(resid_e)
        median_resid_sigma = np.median(resid_sigma)
        mean_resid_sigma = np.mean(resid_sigma)
        std_resid_sigma = np.std(resid_sigma)

        resid_dict = {
            'e1': resid_e1,
            'e2': resid_e2,
            'e': resid_e,
            'sigma': resid_sigma,
            'theta': theta_resid,
            'median_g': median_resid_g,
            'median_e': median_resid_e,
            'median_sigma': median_resid_sigma,
            'mean_sigma': mean_resid_sigma,
            'std_sigma': std_resid_sigma
        }

        return AttrDict(resid_dict)


    def _populate_dicts(self):
        """
        Populate the ellipticity dictionaries for plotting
        """

        psfs  = self.psfs
        stars = self.stars

        wg = (psfs.hsm_g1 > -9999) & (stars.hsm_g1 > -9999)
        self.x = stars.x[wg]
        self.y = stars.y[wg]

        self.star_dict  = self._make_ellip_dict(stars, wg)
        self.psf_dict   = self._make_ellip_dict(psfs, wg)
        self.resid_dict = self._make_resid_dict(self.star_dict, self.psf_dict)


    def _make_quiver_dict(self, vc1, vc2=0, scale=1):
        """
        Set up quiverplot plot dictionary a la PIFF
        Inputs
                sig1: center for star & PSF norms
                sig2: center for residplot norms
        Returns
                q_dict:  dict for quiverplot params
                qkey_dict:  dict for quiverkey params
        """

        scale_units = 'width' # For quiver plots
        norm = colors.CenteredNorm(vcenter=vc1, halfrange=0.06)
        div_norm = colors.CenteredNorm(vcenter=vc2, halfrange=0.05)

        qkey_scale = 0.05
        qkey_label = r'$e_{HSM} = {%.2f}$' % qkey_scale
        fontprops = {'size':14, 'weight':'bold'}

        q_dict = dict(
            cmap='cividis',
            width=90,
            units='xy',
            pivot='mid',
            headaxislength=0,
            headwidth=0,
            headlength=0,
            norm=norm,
            scale=scale,
            scale_units=scale_units
        )

        qkey_dict = dict(
            X=0.2,
            Y=0.02,
            U=qkey_scale,
            labelpos='N',
            label=qkey_label,
            fontproperties=fontprops
        )

        return q_dict, qkey_dict

    def _get_plot_titles(self):

        sd = self.star_dict
        pd = self.psf_dict
        rd = self.resid_dict

        star_title = \
            'median $\sigma^{*}_{HSM} = %.2f$ mas; $e^{*}_{HSM} = %.5f$'\
                        % (sd.median_sigma*1000, sd.median_e)
        psf_title = \
            'median $\sigma^{PSF}_{HSM} = %.2f$ mas; $e^{PSF}_{HSM} = %.5f$'\
                        % (pd.median_sigma*1000, pd.median_e)
        resid_title = \
            'median $\sigma^{resid}_{HSM} = %.2f$ mas; $e^{resid}_{HSM} = %.5f$'\
                        % (rd.median_sigma*1000, rd.median_e)

        return [star_title, psf_title, resid_title]


    def _make_plot(self, dicts, quiver_dict, qkey_dict, titles, scale):
        """
        Make quiverplots, first making quiver_dicts
        """
        # First things first
        set_rc_params(fontsize=14)

        # Custom settings for quiverplot
        plt.rcParams.update({'xtick.direction': 'out'})
        plt.rcParams.update({'ytick.direction': 'out'})
        plt.rcParams.update({'legend.fontsize': 14})

        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True,
                                    figsize=[15, 8], tight_layout=True)

        for i, dc in enumerate(dicts):
            # Set up masks for really circular objects
            min_ellip = 0.01

            # Define the minimum ellipse size and create a mask of circles
            mask = dc.e <= min_ellip
            round = axs[i].scatter(self.x[mask], self.y[mask], s=9,
                        facecolor='black', edgecolors='black')

            mask = dc.e > min_ellip
            q = axs[i].quiver(
                self.x[mask], self.y[mask],
                dc.e1[mask], dc.e2[mask], dc.sigma[mask],
                angles=np.rad2deg(dc.theta[mask]), **quiver_dict
            )
            # adjust x, y limits
            lx, rx = axs[i].get_xlim()
            axs[i].set_xlim(lx-2000, rx+2000)
            ly, ry = axs[i].get_ylim()
            axs[i].set_ylim(ly-500, ry+500)
            key = axs[i].quiverkey(q, **qkey_dict)
            ax_divider = make_axes_locatable(axs[i])
            cax = ax_divider.append_axes("bottom", size="5%", pad="7%")
            cbar = fig.colorbar(q, cax=cax, orientation="horizontal")
            axs[i].set_title(titles[i])

        return fig

    def make_hex_plots(self, outname='hexbins.png'):
        """ Make hex bins for fun and profit """

        # First things first I'm the reallest
        dicts = [self.star_dict, self.psf_dict, self.resid_dict]

        set_rc_params(fontsize=16)

        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True,
                                figsize=[15, 8], tight_layout=True)
        # First do the e1 map
        titles = ['Stars $e_1$',
            'PSF model $e_1$', 'Star-model residuals $e_1$']
        for i, dc in enumerate(dicts):
            if i in [0, 2]:
                vmin = 0.8*np.mean(dc.e1)
                #vmin = np.min(dc.e1)
                vmax = 1.2*np.median(dc.e1)
            else:
                vmin = np.min(dc.e1)
                vmax = np.max(dc.e1)
            im = axs[i].hexbin(self.x, self.y, C=dc.e1, gridsize=(6, 6),
                    cmap=plt.cm.RdYlBu_r)
            axs[i].set_title(titles[i])
            lx, rx = axs[i].get_xlim()
            axs[i].set_xlim(lx-2000, rx+2000)
            ly, ry = axs[i].get_ylim()
            axs[i].set_ylim(ly-1000, ry+1000)
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax)
        fig.savefig(outname.split('hex')[0]+'e1_hex.'+outname.split('hex')[1])

        # Then do e2 map
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True,
                                figsize=[15, 8], tight_layout=True)
        titles = ['Stars $e_2$',
            'PSF model $e_2$', 'Star-model residuals $e_2$']
        for i, dc in enumerate(dicts):
            '''
            if i==0:
                #vmin = 0.8*np.mean(dc.e2)
                #vmin = np.min(dc.e2)
                vmax = 1.2*np.mean(dc.e2)
            else:
                vmin = np.min(dc.e2)
                vmax = np.max(dc.e2)
            '''
            im = axs[i].hexbin(self.x, self.y, C=dc.e2, gridsize=(6, 6),
                    cmap=plt.cm.RdYlBu_r)
            lx, rx = axs[i].get_xlim()
            axs[i].set_xlim(lx-2000, rx+2000)
            ly, ry = axs[i].get_ylim()
            axs[i].set_ylim(ly-1000, ry+1000)
            axs[i].set_title(titles[i])
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax)
        fig.savefig(outname.split('hex')[0]+'e2_hex.'+outname.split('hex')[1])

        # Then sigma map
        fig, axs = plt.subplots(nrows=1, ncols=3, sharey=True,
                                    figsize=[15, 8], tight_layout=True)
        titles = ['Stars $\sigma_{HSM}$',
            'PSF model $\sigma_{HSM}$', 'Star-model residuals $\sigma_{HSM}$']
        vdict = [
            [0.8*np.median(dicts[0].sigma), 1.2*np.median(dicts[0].sigma)],
            [0.8*np.median(dicts[1].sigma), 1.1*np.median(dicts[1].sigma)],
            [None, None],
        ]
        vdict = np.array(vdict)
        for i, dc in enumerate(dicts):
            im = axs[i].hexbin(self.x, self.y, C=dc.sigma, gridsize=(6, 6),
                    cmap=plt.cm.RdYlBu_r, vmin=vdict[i,0], vmax=vdict[i,1])
            lx, rx = axs[i].get_xlim()
            axs[i].set_xlim(lx-2000, rx+2000)
            ly, ry = axs[i].get_ylim()
            axs[i].set_ylim(ly-1000, ry+1000)
            axs[i].set_title(titles[i])
            divider = make_axes_locatable(axs[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax)
        fig.savefig(outname.split('hex')[0]+'sigma_hex.'+outname.split('hex')[1])


    def run(self, scale=1, outname='quiverplot.png'):
        """
        Take a table or whatever, make quiverplot with it
        Filter out failed fits first!
        """

        # Populate the dicts
        self._populate_dicts()

        # Get list of them
        dicts = [self.star_dict, self.psf_dict, self.resid_dict]

        # Get titles (they're defined here!)
        titles = self._get_plot_titles()

        # Get quiver dicts, which set plot params
        vc1 = self.star_dict.median_sigma
        quiver_dict, qkey_dict = self._make_quiver_dict(vc1, scale)

        # Make plot
        fig = self._make_plot(dicts, quiver_dict, qkey_dict, titles, scale)

        # Save to file
        fig.savefig(outname)

        # Bonus round make hexplots!
        self.make_hex_plots(outname.replace('quiver', 'hexgrid'))


        return

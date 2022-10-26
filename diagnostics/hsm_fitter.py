import numpy as np
import galsim,galsim.des



def do_hsm_fit(maker,verbose=False):
    '''
    Get Hirata-Seljak fit moments for PSF and star images. From the GalSim docs:

    This method estimates the best-fit elliptical Gaussian to the object (see Hirata & Seljak 2003
    for more discussion of adaptive moments).  This elliptical Gaussian is computed iteratively
    by initially guessing a circular Gaussian that is used as a weight function, computing the
    weighted moments, recomputing the moments using the result of the previous step as the weight
    function, and so on until the moments that are measured are the same as those used for the
    weight function.

    Might seem strange to take galsim.Image() output of des_psfex and
    piff.draw, take array part of galsim.Image and redraw with a pixel scale
    wcs since the GSObjects already *had* a WCS. However, it can be shown
    that the result of HSM fit is the same in both cases.
    '''

    for i,stamp in enumerate(maker.stamps):
        try:
            gs_object = galsim.Image(stamp, wcs=galsim.PixelScale(maker.pixel_scale))
            HSM_fit=gs_object.FindAdaptiveMom()

            maker.hsm_sig.append(HSM_fit.moments_sigma)
            maker.hsm_g1.append(HSM_fit.observed_shape.g1)
            maker.hsm_g2.append(HSM_fit.observed_shape.g2)
            maker.fwhm.append(gs_object.calculateFWHM())

        except:
            print("HSM fit for stamp #%d failed, skipping" % i)
            maker.hsm_sig.append(-9999)
            maker.hsm_g1.append(-9999)
            maker.hsm_g2.append(-9999)
            maker.fwhm.append(gs_object.calculateFWHM())

    maker.hsm_sig = np.array(maker.hsm_sig)
    maker.hsm_g1 = np.array(maker.hsm_g1)
    maker.hsm_g2 = np.array(maker.hsm_g2)
    maker.fwhm = np.array(maker.fwhm)

    return

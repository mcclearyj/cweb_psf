import webbpsf
from astropy.io import fits
import os
from skimage.transform import rotate

def run_webb_psf(image_file, oversample_lw, outdir=None):
    '''
    Create a WebbPSF object for an image file given its filter and date
    of observation. Save to disk as a FITS file. Note that getting the WebbPSF by
    date of observation requires the astroquery Python package to be installed

    Inputs:
        image_file: filename of image for which WebbPSF should be made.
        oversample_lw: boolean specifying whether or not long-wavelength channel
                        data should be oversampled.
        outdir: directory to which PSF should be saved (default: same as image)

    Outputs:
        FITS image containing appropriate WebbPSF model for image
    '''

    # Get header keywords
    imhead = fits.getheader(image_file, ext=0)
    filter_name = imhead['FILTER']
    date = imhead['DATE']
    detector = imhead['DETECTOR'].replace('LONG', '5')
    roll_angle = fits.getval(image_file, 'ROLL_REF', ext=1)

    # Define output file name
    if outdir is None:
        outdir = os.path.dirname(image_file)
    image_name = os.path.basename(image_file)
    output_file = os.path.join(outdir,
                    image_name.replace('.fits', '_WebbPSF.fits'))
    # Set oversample scale
    if (filter_name in ['F277W', 'F444W']) and (oversample_lw == True):
        oversample = 2
    else:
        oversample = 1

    # create WebbPSF instance for given image parameters
    nc = webbpsf.NIRCam()
    nc.options['parity'] = 'odd'
    nc.filter =  filter_name
    nc.oversample = oversample
    nc.detector = detector

    # Load OSS by date
    nc.load_wss_opd_by_date(date)

    # Calculate PSF & save to file. Make it big so we don't truncate anything.
    # If crop_psf=False, the distorted image will have different dimensions.
    psf = nc.calc_psf(crop_psf=True, fov_pixels=301)
    print(f"nc.oversample is {nc.oversample}")
    psf.writeto(output_file, overwrite=True)

    # NOW! We need to rotate the PSF image by the appropriate roll angle
    for i in range(len(psf)):
        image_data = psf[i].data
        rotated_image = rotate(image_data, (360.0-roll_angle))
        psf[i].data = rotated_image

    psf.writeto(output_file.replace('WebbPSF', 'WebbPSF_rot'), overwrite=True)

    return

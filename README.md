# cweb_psf

PSF characterization code repository for COSMOS-Web lensing working group.
Still under active development! Comments and issues welcome. 

`environment.yml` is @mcclearyj's Conda environment for this project; it may
serve as an example for other users -- TO BE UPDATED

### get_jwst_psf.py ###

Make SExtractor catalogs, star catalogs, and `piff` and `psfex` PSF models for
COSMOS-Web observations in `*i2d` format. Code also adds an `ERR` cutout for
the source into  the star catalog.

Can call `get_jwst_psf.py` with:
```
# Astromatic & PIFF config files
export CONFIGDIR='/your/path/name/cweb_psf/astro_config'
export DATADIR='/your/path/name/data/April2023'
export CODEDIR='/your/path/name/cweb_psf'

get_jwst_psf.py $DATADIR/image1 [image2 ...] -config RUN_CONFIG [-outdir OUTDIR]
          -configdir [CONFIGDIR] [--overwrite] [--help] [--vb]

positional arguments:
  image1 [image2 ...]   Images to process (wildcards OK)

options:
  -h, --help            show this help message and exit
  -outdir OUTDIR        Where to save files [default: ./]
  -configdir CONFIGDIR  Location of SE/PSFEx/PIFF configs [default: astro_config]
  -config CONFIGFILE    Config file for running this script
  --vb                  Print detailed messages [does nothing for now]
```
Dependencies:
numpy, re, os, astropy, matplotlib, pdb, glob, argparse, esutil

### master_psf_diagnostic.py ###

Makes number of diagnostic figures and tables for input PSFEx and/or PIFF PSF models.
Imports classes from `psfmaker.py` and `starmaker.py` in `diagnostics`

Options include:
- Test a PSFEx model [`-psfex_model`]
- Test a PSFEx model using esheldon `psfex` package: [`--epsfex`] (default)
- Test a PSFEx model using `galsim` rendering [`--gpsfex`] (needs `-im_name` too)
- Test a PIFF model [`-piff_model`]
- Test a PSFEx model with no spatial variation [`-single_model`]
- Test a WebbPSF model [`-webb_model`]

Can call `master_psf_diagnostic.py` with

```
master_psf_diagnostic.py basedir star_cat [-h] [-min_snr MIN_SNR]
    [-pix_scale PIX_SCALE] [-outdir OUTDIR] [-psfex_name PSFEX_NAME]
    [-im_name IM_NAME] [-piff_model PIFF_NAME] [--epsfex]
    [--gpsfex] [-single_model SIMPLE_PSFEX_MODEL]
    [-webb_model WEBBPSF_MODEL] [-vignet_size] [--verbose]

positional arguments:
  basedir               Directory containing star catalogs & images
  star_cat              Cleaned-up star catalog to use for PSF diagnostic

options:
  -h, --help            show this help message and exit
  -min_snr MIN_SNR      Optional S/N cut for star catalog [default=None]
  -pix_scale PIX_SCALE  Image/PSF pixel scale [default=0.03]
  -outdir OUTDIR        Output directory for diagnostics [default=./psf_diagnostics]
  -psfex_name PSFEX_NAME
                        PSFEx model filename
  -im_name IM_NAME      FITS image filename for GalSim PSFEx diagnostic
  -piff_model PIFF_NAME PIFF psf model filename
  -single_model PSFEX_MODEL
                        PSFEx model with no spatial variation
  -webb_model WEBBPSF_MODEL WebbPSF model
  -vignet_size          Sub-section of star & PSF stamps to consider
  --epsfex              Run esheldon psfex diagnostic (default)
  --gpsfex              Run galsim.des_psfex diagnostic
  --verbose, -v         Verbosity

```
Dependencies:
  numpy, re, os, astropy, matplotlib, pdb, argparse, galsim, galsim.des, piff, glob, treecorr, psfex

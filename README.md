# cweb_psf

PSF characterization code repository for COSMOS-Web lensing working group.

`environment.yml` is @mcclearyj's Conda environment for this project; it may
serve as an example for other users -- TO BE UPDATED

### get_jwst_psf.py ###

Quick-and-dirty script to make `piff` PSF models for
COSMOS-Web simulated observations. Currently runs `piff` in single-image mode.
Still under active development.  

Can call `get_jwst_psf.py` with:
```
export CONFIGDIR='/your/path/name/cweb_psf/astro_config'
export DATADIR='/your/path/name/mock_data/DEC2022_v.0.0.1'
export OUTDIR='/your/path/name/mock_data/DEC2022_v.0.0.1/DEC2022_v.0.0.1/working'

get_jwst_psf.py images [images ...] [-outdir OUTDIR] [-configdir CONFIGDIR]
    [-truthstars TRUTHSTARS] [--overwrite] [--help] [--vb]

positional arguments:
  images                Images to process (wildcards OK)

options:
  -h, --help            show this help message and exit
  -outdir OUTDIR        Where to save files
  -configdir CONFIGDIR  Location of SEx/PIFF config files
  -truthstars TRUTHSTARS
                        Star catalog to use for PSF fits
  --overwrite           Overwrite sci/weight files
  --vb                  Print detailed messages [does nothing for now]
```
Dependencies:
numpy, re, os, astropy, matplotlib, pdb, glob, argparse

### master_psf_diagnostic.py ###

Makes number of diagnostic figures and tables for input PSFEx and/or PIFF PSF models.
Imports classes from `psfmaker.py` and `starmaker.py` in `diagnostics`

Options include:
- Test a PSFEx model using esheldon `psfex` package: [`--epsfex`]
- Test a PSFEx model using `galsim` rendering [`--gpsfex`]
- Test a PIFF model [`--piff`]

Can call `master_psf_diagnostic.py` with

```
master_psf_diagnostic.py basedir star_cat [-h] [-min_snr MIN_SNR]
    [-pix_scale PIX_SCALE] [-outdir OUTDIR] [-psfex_name PSFEX_NAME]
    [-im_name IM_NAME] [-piff_name PIFF_NAME] [--epsfex]
    [--gpsfex] [--piff] [--noisefree] [--verbose]

positional arguments:
  basedir               Directory containing star catalogs & images
  star_cat              Star catalog to use for PSF diagnostic

options:
  -h, --help            show this help message and exit
  -min_snr MIN_SNR      Optional S/N cut for star catalog [default=None]
  -pix_scale PIX_SCALE  Image/PSF pixel scale [default=0.03]
  -outdir OUTDIR        Output directory for diagnostics [default=./psf_diagnostics]
  -psfex_name PSFEX_NAME
                        PSFEx model filename
  -im_name IM_NAME      FITS image filename for GalSim PSFEx diagnostic
  -piff_name PIFF_NAME  PIFF psf model filename
  --epsfex              Run esheldon psfex diagnostic
  --gpsfex              Run galsim.des_psfex diagnostic
  --piff                Run PIFF diagnostic
  --noisefree           Disable adding noise to PSF stamps
  --verbose, -v         Verbosity

```
Dependencies:
  numpy, re, os, astropy, matplotlib, pdb, argparse, galsim, galsim.des, piff, glob, treecorr, psfex

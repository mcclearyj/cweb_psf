# cweb_psf

PSF characterization code repository for COSMOS-Web lensing working group.

`environment.yml` is @mcclearyj's Conda environment for this project; it may
serve as an example for other users.

####`get_jwst_psf.py`####
Quick-and-dirty script to make `piff` PSF models for
COSMOS-Web simulated observations. Currently runs `piff` in single-image mode.
Still under active development.  

Can call `get_jwst_psf.py` with:
```
export CONFIGDIR='/your/path/name/cweb_psf/astro_config'
export DATADIR='/your/path/name/mock_data/DEC2022_v.0.0.1'
export OUTDIR='/your/path/name/mock_data/DEC2022_v.0.0.1/DEC2022_v.0.0.1/working'

python get_jwst_psf.py -basedir $DATADIR -outdir $OUTDIR -configdir $CONFIGDIR
```
Dependencies:
  numpy, re, os, astropy, matplotlib, pdb, glob, argparse

#### `master_psf_diagnostic.py` ####
Makes number of diagnostic figures and tables for input PSFEx and/or PIFF PSF models.
Calls the `psfmaker.py` and `starmaker.py`.

Options include:
- Test a PSFEx model using esheldon psfex package: [`--epsfex`]
- Test a PSFEx model using `galsim` rendering [`--gpsfex`]
- Test a PIFF model [`--piff`]

Can call `master_psf_diagnostic.py` with

```
  master_psf_diagnostic.py [-h] [-min_snr MIN_SNR] [-outdir OUTDIR]
      [-psfex_name PSFEX_NAME] [-im_name IM_NAME] [-piff_name PIFF_NAME]
      [--gpsfex] [--epsfex] [--piff]
      [--noisefree] [--verbose] imdir star_cat
```
Dependencies:
  numpy, re, os, astropy, matplotlib, pdb, argparse, galsim. galsim.des, piff, glob, treecorr

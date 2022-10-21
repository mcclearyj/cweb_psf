# cweb_psf

Can call `get_jwst_psf.py` with:
```
conda activate your_conda_env

export CONFIGDIR='/your/path/name/cweb_psf/astro_config'
export DATADIR='/your/path/name/mock_data/DEC2022_v.0.0.1'
export OUTDIR='/your/path/name/mock_data/DEC2022_v.0.0.1/DEC2022_v.0.0.1/working'

python get_jwst_psf.py -basedir $DATADIR -outdir $OUTDIR -configdir $CONFIGDIR
```

Can call the `master_psf_diagnostic.py` with

```
  master_psf_diagnostic.py [-h] [--min_snr MIN_SNR] [--outdir OUTDIR]
      [--epsfex] [--gpsfex] [--piff] [--psfex_name PSFEX_NAME]
      [--im_name IM_NAME] [--piff_name PIFF_NAME]
      [--noisefree] [--verbose] imdir star_cat
```

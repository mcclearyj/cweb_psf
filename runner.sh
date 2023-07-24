#!/bin/bash

export CODEDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf'
export CONFIGDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf/astro_config'
export DATADIR='/Users/j.mccleary/Research/jwst_cosmos/real_data/Apr2023/jw01727128001'
export OUTDIR='working'


# Grab all files matching the pattern "jw017*_cal.fits"
files=$(ls $DATADIR/jw01727128001_04101*nrcb*cal.fits 2>/dev/null)

# Check if any files match the pattern
if [ -z "$files" ]; then
  echo "No files found matching the pattern jw017*_cal.fits"
  exit 1
else
  echo "Working on files:"
  echo $files
fi

# Run get_jwst_psf on the files
python $CODEDIR/get_jwst_psf.py -config $CODEDIR/'configs/apr2023_config.yaml'  $files


# Loop through each file
for file in $files; do
  # Get the filename without extension
  filename=$(basename "$file" .fits)

  # Run the command for each file
  python $CODEDIR/master_psf_diagnostic.py . \
    $OUTDIR/"$filename"_pex_stars.fits \
    -min_snr 5 -pix_scale 0.03 -outdir $OUTDIR/psf_diagnostic_plots/"$filename" \
    -psfex_model $OUTDIR/psfex-output/"$filename"/"$filename"_starcat.psf \
    --epsfex -webb_model $DATADIR/"$filename"_WebbPSF.fits -vignet_size 75

done

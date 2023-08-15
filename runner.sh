#!/bin/bash

#export CODEDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf'
#export CONFIGDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf/astro_config'
#export DATADIR='/Users/j.mccleary/Research/jwst_cosmos/real_data/Apr2023/jw01727128001'
export OUTDIR='working'


# Grab all files matching the pattern "jw017*_cal.fits"
files=$(ls $DATADIR/jw01727*cal.fits 2>/dev/null)

# Check if any files match the pattern
if [ -z "$files" ]; then
  echo "No files found matching the pattern jw017*_cal.fits"
  exit 1
else
  echo "Working on files:"
  echo $files
fi

# Run get_jwst_psf on the files
for file in $files; do
  python get_jwst_psf.py -config 'configs/box_cutter.yaml' $file
  if [ $? -ne 0 ]; then
    echo "Error processing file $file. Continuing with the next file."
  fi
done
#python get_jwst_psf.py -config 'configs/box_cutter.yaml'  $files


done

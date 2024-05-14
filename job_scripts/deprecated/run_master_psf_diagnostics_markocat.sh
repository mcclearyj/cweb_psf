#!/bin/bash

# Set up some environmental variables

export USERDIR="/Users/j.mccleary/Research"
export BASEDIR="${USERDIR}/jwst_cosmos/real_data/Jan2024/marko_psfex"

export CODEDIR="$USERDIR/jwst_cosmos/cweb_psf"
export DATADIR="$BASEDIR"
export OUTDIR="$BASEDIR/psf_diagnostic_plots"


## Define parent directory

if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

##
## The parentheses define a bash array -- necessary for expanding the wildcard.
## The braces around TARGET show we are accessing that variable (disambiguation)

IM_CATS=(${DATADIR}/${TARGET}*_cal_starcat.fits)
COADD_IM_NAME="${COADDDIR}/${TARGET}_coadd_${band}_starcat.fits"

#IM_CATS+=(${COADD_IM_NAME})

## Flexing

ls ${IM_CATS[@]}
ls $COADD_IM_NAME

##
## OK, let's walk through and make a separate output directory for every
## exposure, then run our master PSF diagnostic code

export IM_CAT="COSMOS-Web-PSFEx-SE-Cat_f277w_B4_starcat.fits"
export DATA_DR="/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2024/marko_psfex"

for IM_CAT in "${IM_CATS[@]}"; do

    echo ""
    echo "Working on ${IM_CAT}..."

    # isolate base name
    IMCAT_BASE=$(basename "$IM_CAT" .fits)
    echo $IMCAT_BASE

    # create an output directory
    CAT_OUTDIR=${OUTDIR}/${IMCAT_BASE}


    python $CODEDIR/master_psf_diagnostic.py \
           $DATADIR $IM_CAT \
           -min_snr 10 -pix_scale 0.03 -outdir $OUTDIR \
           -single JWST-PSFEx_out_f277w_B4_psf_v4.0.psf \
           -vignet_size 101
done

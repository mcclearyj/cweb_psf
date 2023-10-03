#!/bin/bash

# Set up some environmental variables


export TARGET='MACSJ1931'

#export TARGET='Abell3827'
export band='b'
export USERDIR="/Users/j.mccleary/Research"
export BASEDIR="${USERDIR}/SuperBIT/real_data"

export CODEDIR="$USERDIR/jwst_cosmos/cweb_psf"
export DATADIR="$BASEDIR/$TARGET/$band/cat"
export COADDDIR="$BASEDIR/$TARGET/$band/coadd"
export OUTDIR="$DATADIR/psf_diagnostic_plots"


## Define parent directory

if [ ! -d "$OUTDIR" ]; then
    mkdir -p "$OUTDIR"
fi

##
## The parentheses define a bash array -- necessary for expanding the wildcard.
## The braces around TARGET show we are accessing that variable (disambiguation)

IM_CATS=(${DATADIR}/${TARGET}*_clean_starcat.fits)
COADD_IM_NAME="${COADDDIR}/${TARGET}_coadd_${band}_starcat.fits"

#IM_CATS+=(${COADD_IM_NAME})

## Flexing

ls ${IM_CATS[@]}
ls $COADD_IM_NAME

## 
## OK, let's walk through and make a separate output directory for every
## exposure, then run our master PSF diagnostic code

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
           -min_snr 10 -pix_scale 0.141 -outdir $CAT_OUTDIR \
           -psfex_model "${DATADIR}/${IMCAT_BASE}.psf"  --epsfex \
           -vignet_size 51
done


##
## Then also do the coadd catalog

echo ""
echo "Working on ${COADD_IM_NAME}..."

# isolate base name
IMCAT_BASE=$(basename "$COADD_IM_NAME" .fits)

# create an output directory
CAT_OUTDIR=${OUTDIR}/${IMCAT_BASE}

python $CODEDIR/master_psf_diagnostic.py \
       $COADDDIR $COADD_IM_NAME \
       -min_snr 10 -pix_scale 0.141 -outdir $CAT_OUTDIR \
       -psfex_model "${COADDDIR}/${IMCAT_BASE}.psf"  --epsfex \
       -vignet_size 51

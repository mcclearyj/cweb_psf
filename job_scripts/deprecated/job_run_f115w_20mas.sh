#!/bin/bash
#PBS -S /bin/sh
#PBS -N lw_mosaic_psfs
#PBS -o lw_mosaic_psfs.out 
#PBS -j oe 
#PBS -l nodes=1,walltime=23:59:59
#PBS -M j.mccleary@northeastern.edu

source /n23data1/mccleary/miniconda3/bin/activate 
conda activate cweb_psf

export DATADIR="/n23data1/mccleary/real_images/20mas_resamp"
export BANDPASSES="f115w f150w"

export CODEDIR="/n23data1/mccleary/cweb_psf"
export CONFIG=${CODEDIR}/configs/20mas_mosaics_config.yaml

for bandpass in $BANDPASSES; do
    IMAGES=(${DATADIR}/${bandpass}/visit_*_nircam_${bandpass}_COSMOS-Web_20mas_*_i2d.fits)

    if [ ${#IMAGES[@]} -gt 0 ]; then
        echo ""
        echo "Input images: ${IMAGES[*]}"
        echo "Using config: $CONFIG"
        echo ""

        python ${CODEDIR}/get_galaxy_cutouts.py -c ${CONFIG} "${IMAGES[@]}"

    else
        echo "No valid files found for ${bandpass}"
    fi

    echo "done with ${bandpass}"
done

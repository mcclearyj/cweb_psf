#!/bin/bash
#PBS -S /bin/sh
#PBS -N sw_mosaic_psfs
#PBS -o sw_mosaic_psfs.out 
#PBS -j oe 
#PBS -l nodes=1,walltime=6:00:00
#PBS -M j.mccleary@northeastern.edu

source /n23data1/mccleary/miniconda3/bin/activate 
conda activate cweb_psf

export DATADIR="/n23data2/mfranco/jwst/COSMOS-Web_april_2023/mosaics/v0_002"
export BANDPASSES="f115w f150w"
TILES=("A1" "A5" "A6" "A10")
export CODEDIR="/n23data1/mccleary/cweb_psf"
export CONFIG=${CODEDIR}/configs/apr_mosaics_cutout_config.yaml

for bandpass in $BANDPASSES; do
    IMAGES=()
    for tile in "${TILES[@]}"; do
        IMAGE=${DATADIR}/mosaic_nircam_${bandpass}_COSMOS-Web_30mas_${tile}_v0_002_i2d.fits
        if [ -f "$IMAGE" ]; then
            IMAGES+=("$IMAGE")
        else
            echo "File not found: $IMAGE"
        fi
    done

    if [ ${#IMAGES[@]} -gt 0 ]; then
        echo ""
        echo "Input images: ${IMAGES[*]}"
        echo "Using config: $CONFIG"
        echo ""

        python ${CODEDIR}/get_galaxy_cutouts.py --config ${CONFIG} "${IMAGES[@]}"
    else
        echo "No valid files found for ${bandpass}"
    fi

    echo "done with ${bandpass}"
done

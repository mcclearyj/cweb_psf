#!/bin/bash
#PBS -S /bin/sh
#PBS -N sim_mosaic_psfs
#PBS -o sim_mosaic_psfs.out 
#PBS -j oe 
#PBS -l nodes=1,walltime=23:00:00
#PBS -M j.mccleary@northeastern.edu

source /n23data1/mccleary/miniconda3/bin/activate 
conda activate cweb_psf

export DATADIR="/n23data1/mfranco/jwst/COSMOS-Web/DEC2022_v.0.0.4/mosaics"
export BANDPASSES="f277w f444w"
export CODEDIR="/n23data1/mccleary/cweb_psf"
export CONFIG=${CODEDIR}/configs/sim_mosaics_cutout_config.yaml 

# Initialize image array
IMAGES=()

# Populate it with paths to mosaics
for bandpass in $BANDPASSES; do
    IMAGE=${DATADIR}/mosaic_nircam_${bandpass}_COSMOS-Web_i2d.fits
    if [ -f "$IMAGE" ]; then
        IMAGES+=("$IMAGE")
    else
        echo "File not found: $IMAGE"
    fi
done


# If at least one image is found, run code
if [ ${#IMAGES[@]} -gt 0 ]; then
    echo ""
    echo "Input images: ${IMAGES[*]}"
    echo "Using config: $CONFIG"
    echo ""

    python ${CODEDIR}/get_galaxy_cutouts.py --config ${CONFIG} "${IMAGES[@]}"
    
else
    echo "No valid files found!"
fi

echo "All done"


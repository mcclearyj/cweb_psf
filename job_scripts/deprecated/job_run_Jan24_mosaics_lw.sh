#!/bin/bash
#PBS -S /bin/sh
#PBS -N starcats_lw
#PBS -o ./starcats_lw_ims.out 
#PBS -j oe 
#PBS -l nodes=1,walltime=3:00:00

module () {
  eval $(/usr/bin/modulecmd bash $*)
}

source /n23data1/mccleary/miniconda3/bin/activate 
conda activate cweb_psf

export DATADIR="/n23data2/mfranco/jwst/COSMOS-Web_jan_2024/mosaics/12_jan_24"
export CODEDIR="/n23data1/mccleary/cweb_psf"
export CONFIG=${CODEDIR}/"configs/real_mosaics_config.yaml"
export OUTDIR="/n23data1/mccleary/real_data/Jan2024/"

export BANDPASSES="f277w f444w"

for bandpass in $BANDPASSES; do

    export IMAGES=(${DATADIR}/${bandpass}/"*_COSMOS-Web_30mas_B*_i2d.fits")
    export 
    if [ ${#IMAGES[@]} -gt 0 ]; then
        echo ""
        echo "Input images: ${IMAGES[*]}"
        echo "Using config: $CONFIG"
        echo ""
	
	python ${CODEDIR}/get_galaxy_cutouts.py -config ${CONFIG} -outdir ${OUTDIR}/${bandpass} ${IMAGES[@]}
    else
        echo "No valid files found for ${bandpass}"
    fi

    echo "done with ${bandpass}"
done


echo "done"

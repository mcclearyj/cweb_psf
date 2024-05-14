export DATADIR="/n23data1/mfranco/jwst/COSMOS-Web/DEC2022_v.0.0.4/mosaics"
export BANDPASSES="f115w f150w f277w f444w"
export CODEDIR="/n23data1/mccleary/cweb_psf"
export CONFIG=${CODEDIR}/configs/sim_mosaics_cutout_config.yaml 

# Initialize image array
IMAGES=()

# Fill image array with images
for bandpass in $BANDPASSES; do
    IMAGE=${DATADIR}/mosaic_nircam_${bandpass}_COSMOS-Web_i2d.fits
    if [ -f "$IMAGE" ]; then
        IMAGES+=("$IMAGE")
    else
        echo "File not found: $IMAGE"
    fi
done 

# Run code if there is at least one image
if [ ${#IMAGES[@]} -gt 0 ]; then
    echo ""
    echo "Input images: ${IMAGES[*]}"
    echo "Using config: $CONFIG"
    echo ""

    echo "python ${CODEDIR}/get_galaxy_cutouts.py --config ${CONFIG} ${IMAGES[@]}"
    
else
    echo "No valid files found for ${bandpass}"
fi

echo "done with ${bandpass}"


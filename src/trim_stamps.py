###
### Functions to trimp star and PSF stamps to be same size
###
import numpy as np
from src.star_psf_holder import StarPSFHolder

def _trim_stamps(stamp):
    # What shape is the stamp?
    stamp_dim = stamp.shape[0]

    # Does stamp need trimming?
    if stamp_dim > self.vignet_size:
        ns = (stamp_dim - self.vignet_size) // 2
        ks = ns; js = -ns
        if self.vb == True:
            print(f'Trimming {ns} pixels from {self.psf_type} stamp')
    else:
        ks = 0; js = None

    # Trim (or don't) and return
    return stamp[ks:js, ks:js]

def trim_stamps(star_holder, psf_holder):
    """
    Trim stamps to be vignet size; in the absence of a supplied vignet
    size argument, just take the smaller of the PSF and star image.

    Input
        star_holder: should be a StarPSFHolder instance
    Returns
        trimmed star stamps too

    TO DO: allow for a PSF dim < star dim!
    """
    stars = star_holder.stamps; psfs = psf_holder.stamps
    star_dim = stars[0].shape[0]; psf_dim = psfs[0].shape[0]
    star_vs = star_holder.vignet_size; psf_vs = psf_holder.vignet_size
    if psf_vs == None:
        if star_vs == None:
            psf_holder.vignet_size = min(star_dim, psf_dim)
            star_holder.vignet_size = min(star_dim, psf_dim)
        else:
            psf_holder.vignet_size = star_vs

    # Trim stars and PSF images if needed
    trimmed_stars = list(
        map(_trim_stamps(), stars)
    )
    trimmed_psfs = list(
        map(_trim_stamps(), psfs)
    )

    # Replace objects with their trimmed versions
    star_holder.stamps = np.array(trimmed_stars)
    psf_holder.stamps = np.array(trimmed_psfs)

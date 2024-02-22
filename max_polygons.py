import numpy as np
from astropy.io import ascii as at
from scipy.spatial import ConvexHull
from astropy import wcs
import matplotlib.path as mpltPath
import os
import re
import functools

def which_visit_OG(ra, dec):
    """
    Determine which COSMOS-Web visit covers a particular coordinate by
    constructing a convex hull.

    From YouTube (https://www.youtube.com/watch?v=B2AJoQSZf4M): given a set of
    points on a 2 dimensional plane, a convex cull is a geometric object, i.e.,
    a polygon, that encloses all of those points. The vertices of this polygon
    maximize the polygon's area while minimizing its circumference. If
    extraneous points were to be included in the polygon, then its covering area
    is reduced while the circumference is increased. Likewise, if some of
    the vertices of the convex hull are omitted, then the resulting polygon
    won’t cover all of the points in the set.
    """

    points = np.zeros((1, 2))
    points[0, 0] = ra
    points[0, 1] = dec

    # This file contains the vertices of the footprints of each visit
    coords_file = at.read('./coords_f277w_CW_JAN2024.txt')

    # Loop over each entry in coordinate file...
    for ii in range(np.size(coords_file)):

        # Initialize array that will store 4 x 2 vertex coords.
        coords_pointings = np.zeros((4, 2))

        # Set up RA pointings
        coords_pointings[:, 0] = np.array([
            coords_file['x1'][ii], coords_file['x2'][ii],
            coords_file['x4'][ii], coords_file['x3'][ii]
        ])

        # Set up DEC pointings
        coords_pointings[:, 1] = np.array([
            coords_file['y1'][ii], coords_file['y2'][ii],
            coords_file['y4'][ii], coords_file['y3'][ii]
        ])

        # Create the convex hull that covers these four points
        hull_pointings = ConvexHull(coords_pointings)

        ploygon_hull_pointings = np.zeros(
            (np.size(coords_pointings[hull_pointings.vertices, 0]), 2)
        )
        ploygon_hull_pointings[:,0] = \
            coords_pointings[hull_pointings.vertices, 0]
        ploygon_hull_pointings[:,1] = \
            coords_pointings[hull_pointings.vertices, 1]

        path_pointings = mpltPath.Path(ploygon_hull_pointings)
        inside_pointings = path_pointings.contains_points(points)

        if inside_pointings == True:
            print(os.path.basename(coords_file['name'][ii]))

def _which_visit(coords_file, this_ra, this_dec):
    """
    Determine which COSMOS-Web visit covers a particular coordinate specified
    with (this_ra, this_dec) by constructing a convex hull.

    Lightly adapted from Max Franco's original function, included above as
    which_visit_OG().

    From YouTube (https://www.youtube.com/watch?v=B2AJoQSZf4M): given a set of
    points on a 2 dimensional plane, a convex cull is a geometric object, i.e.,
    a polygon, that encloses all of those points. The vertices of this polygon
    maximize the polygon's area while minimizing its circumference. If
    extraneous points were to be included in the polygon, then its covering area
    is reduced while the circumference is increased. Likewise, if some of
    the vertices of the convex hull are omitted, then the resulting polygon
    won’t cover all of the points in the set.
    """

    points = np.zeros((1, 2))
    points[0, 0] = this_ra
    points[0, 1] = this_dec

    # Loop over each entry in coordinate file...
    for ii in range(np.size(coords_file)):

        # Initialize array that will store 4 x 2 vertex coords.
        coords_pointings = np.zeros((4, 2))

        # Set up RA pointings
        coords_pointings[:, 0] = np.array([
            coords_file['x1'][ii], coords_file['x2'][ii],
            coords_file['x4'][ii], coords_file['x3'][ii]
        ])

        # Set up DEC pointings
        coords_pointings[:, 1] = np.array([
            coords_file['y1'][ii], coords_file['y2'][ii],
            coords_file['y4'][ii], coords_file['y3'][ii]
        ])

        # Create the convex hull that covers these four points
        hull_pointings = ConvexHull(coords_pointings)

        ploygon_hull_pointings = np.zeros(
            (np.size(coords_pointings[hull_pointings.vertices, 0]), 2)
        )
        ploygon_hull_pointings[:,0] = \
            coords_pointings[hull_pointings.vertices, 0]
        ploygon_hull_pointings[:,1] = \
            coords_pointings[hull_pointings.vertices, 1]

        path_pointings = mpltPath.Path(ploygon_hull_pointings)
        inside_pointings = path_pointings.contains_points(points)

        if inside_pointings == True:
            break

    return os.path.basename(coords_file['name'][ii])

def which_visit(coords_file, ras, decs):
    """
    For a given Max Franco-format coords_file list of visit names, return
    the visit number to which each RA/Dec pair belongs. Uses map() for
    maximum efficiency!

    Inputs
        coords_file: list of visit names with footprint edges (can be string)
        ras: list of right ascensions
        decs: list of declinations
    Returns
        visits: list of three-digit visits
    """

    ## Read in coords file if needed
    if type(coords_file) == str:
        coords_file = at.read(coords_file)

    ## Grab full visit names from coords_file
    visit_names = list(map(
        functools.partial(_which_visit, coords_file), ras, decs
    ))

    ## Now isolate visit
    visits = []
    for visit_name in visit_names:
        visits.append(re.search(r"jw01727(\d){3}", visit_name).group()[-3:])
    return visits

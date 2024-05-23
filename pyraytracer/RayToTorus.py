# % function [intersection_points surface_normals distance_traveled crossing_into] = ...
# %     RayToTorus(starting_points, incoming_directions, torus_center, torus_axis, torus_r1, torus_r2)

# % RayToTorus specifics:
# %
# % other inputs:
# %       torus_center         -   3-element vector, giving the point at
# %                                     the torus center
# %       torus_axis           -   3-element vector pointing along the
# %                                     torus axis
# %       torus_r1             -   scalar giving torus large radius
# %
# %       torus_r2             -   scalar giving torus small radius (height)
# %
# %
# % 9/3/15, CED

import sys
import numpy as np
import warnings

def RayToTorus(starting_points, indir, torus_center, torus_axis, r1, r2):
    params = len(locals())

    intersection_points = []
    surface_normals = []
    distance_traveled = []
    crossing_into = []

    # check inputs
    if params < 4 or len(torus_center) != 3 or len(torus_axis) != 3 or not isinstance(r1, (int, float)) or not isinstance(r2, (int, float)) or \
        starting_points.shape[1] != 3 or indir.shape[1] != 3 or starting_points.shape[0] != indir.shape[0]:
        raise Exception('Improper input to RayToTorus')
    torus_center = np.transpose(torus_center[:,np.newaxis])
    torus_axis = np.transpose(torus_axis[:,np.newaxis])
    numrays = starting_points.shape[0]

    if np.sum(torus_axis**2) > 0:
        torus_axis = torus_axis / np.abs(np.sqrt(np.sum(torus_axis**2)))
    else:
        raise Exception('Invalid torus for RayToTorus')

    # %% solve quadratic for distance_traveled
    # % ( || u + l*v || - r1 )^2 + ( || y + l*w || )^2 = r2^2
    # % ||u+lv||^2 + r1^2 + ||y+lw||^2 - r2^2 = 2*r1*||u+lv||
    # % [ uu + 2luv + llvv + r1^2 + yy + 2lyw + llww - r2^2 ] ^2 = 4 * r1^2 * (uu + 2luv + llvv)
    # % a4*l^4 + a3*l^3 + a2*l^2 + a1*l + a0 = 0
    x = starting_points - torus_center
    y = (x @ np.transpose(torus_axis) @ torus_axis)
    u = y - x
    w = (indir @ np.transpose(torus_axis) @ torus_axis)
    v = w - indir

    a4 = (np.sum(v**2, axis=1) + np.sum(w**2, axis=1)) ** 2
    a3 = 4 * (np.sum(v**2, axis=1) + np.sum(w**2, axis=1)) * (np.sum(u*v, axis=1) + np.sum(y*w, axis=1))
    a2 = 4 * (np.sum(u*v, axis=1) + np.sum(y*w, axis=1)) ** 2 - (4 * r1**2 * np.sum(v**2, axis=1)) + \
         (2 * (np.sum(v**2, axis=1) + np.sum(w**2, axis=1))) * (np.sum(u**2, axis=1) + r1**2 + np.sum(y**2, axis=1) - r2**2)
    a1 = 4 * (np.sum(u*v, axis=1) + np.sum(y*w, axis=1)) * \
         (np.sum(u**2, axis=1) + r1**2 + np.sum(y**2, axis=1) - r2**2) - (8 * r1**2 * np.sum(u*v, axis=1))
    a0 = (np.sum(u**2, axis=1) + r1**2 + np.sum(y**2, axis=1) - r2**2)**2 - 4 * r1**2 * np.sum(u**2, axis=1)

    a = np.array([a4, a3, a2, a1, a0]) # I believe matlab treats this as a nx5, numpy as a 5xn
    a = np.transpose(a)

    quartic_cut = np.not_equal(a[:, 0], 0)
    cubic_cut = np.logical_and(~quartic_cut, np.not_equal(a[:, 1], 0))
    quad_cut = np.logical_and(~np.logical_or(quartic_cut, cubic_cut), np.not_equal(a[:, 2], 0))
    linear_cut = np.logical_and(~np.logical_or.reduce((quartic_cut, cubic_cut, quad_cut)), np.not_equal(a[:, 3], 0))
    nan_cut = ~np.logical_or.reduce((quartic_cut, cubic_cut, quad_cut, linear_cut))

    distance_traveled = np.zeros((numrays, 4))
    distance_traveled[nan_cut, :] = np.nan

    if np.any(linear_cut):
        distance_traveled[linear_cut, 0:2] = -a[linear_cut, 4] / a[linear_cut, 3]
        distance_traveled[linear_cut, 2:4] = np.nan

    if np.any(quad_cut):
        distance_traveled[quad_cut, 0:2] = np.tile(-0.5 * a[quad_cut, 3] / a[quad_cut, 2], (1, 2)) + \
                                           ((0.5 * np.sqrt(a[quad_cut, 3]**2 - 4 * a[quad_cut, 2] * a[quad_cut, 4]) /
                                             a[quad_cut, 2])) * np.array([1, -1])
        distance_traveled[quad_cut, 2:4] = np.nan

    if np.any(cubic_cut):
        for ix in np.transpose(np.nonzero(cubic_cut)):
            distance_traveled[ix, 0:3] = np.roots(a[ix, 1:5].flatten())
        distance_traveled[cubic_cut, 3] = np.nan

    if np.any(quartic_cut):
        for ix in np.transpose(np.nonzero(quartic_cut)):
            warnings.filterwarnings('ignore')  # ignore ComplexWarning here, want to get rid of imaginary parts
            distance_traveled[ix, :] = np.roots(a[ix, :].flatten())

    # find intersection_points
    intersection_points = starting_points[:, :, np.newaxis] + distance_traveled[:, np.newaxis, :] * indir[:, :, np.newaxis] # check reshaping

    # find surface_normals
    surface_normals = np.zeros(intersection_points.shape)
    for i_l in range(4):
        x = intersection_points[:, :, i_l] - torus_center
        y = (x @ np.transpose(torus_axis) @ torus_axis)
        u = x - y
        surface_normals[:, :, i_l] = (1 / r2) * (u * (1 - (r1 / np.sqrt(np.sum(u**2, axis=1))[:,np.newaxis])) + y)

    crossing_into = np.round_(-np.sign(np.sum(indir[:, :, np.newaxis] * surface_normals, axis=1)))
    surface_normals = surface_normals * crossing_into[:, np.newaxis, :]
    crossing_into = np.reshape(crossing_into, (-1, 4))

    return [intersection_points, surface_normals, distance_traveled, crossing_into]


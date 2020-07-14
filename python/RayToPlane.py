import numpy as np
import numpy.matlib


def RayToPlane(starting_points, indir, plane_pt, plane_norm):

    intersection_points = []
    surface_normals = []
    distance_traveled = []
    crossing_into = []

    # Check inputs
    if len(plane_pt) != 3 or len(plane_norm) != 3 or starting_points.shape[1] != 3 or indir.shape[1] != 3 or starting_points.shape[0] != indir.shape[0]:
        raise Exception('Improper input to RayToPlane')
    plane_pt = np.transpose(plane_pt)
    plane_norm = np.transpose(plane_norm)
    numrays = starting_points.shape[0]

    # Normalize directions | some of this should probably be done in RayToShape as it is common among all RayToXXX
    goodray_cut = np.sum(indir ** 2, 1) > 0
    if any(goodray_cut):
        indir[goodray_cut, :] = indir[goodray_cut, :] / np.transpose(np.matlib.repmat(abs(np.sqrt(np.sum(indir ** 2, 1))), 3, 1))
    if sum(plane_norm ** 2) > 0:
        plane_norm = plane_norm / abs(np.sqrt(np.sum(plane_norm**2)))
    else:
        raise Exception('Invalid plane for RayToPlane')

    # Find intersection points
    # distance_traveled = ((repmat(plane_point, numrays, 1) - starting_points) * plane_normal') ./ (incoming_directions * plane_normal');
    with np.errstate(divide='ignore', invalid='ignore'):
        distance_traveled = np.sum((plane_pt - starting_points) * plane_norm, axis=1, keepdims=True) / np.sum(indir * plane_norm, axis=1, keepdims=True)
        # distance_traveled = ((np.matlib.repmat(plane_pt, numrays, 1) - starting_points) * plane_norm) / (indir * plane_norm)

    # Removed repmat for distance_traveled and np.sign(indir * plane_norm), were already nx3; should it be nx3 or nx1?
    intersection_points = starting_points[:,:,np.newaxis] + distance_traveled[:,np.newaxis,:] * indir[:,:,np.newaxis]
    surface_normals = -plane_norm * np.sign(np.sum(indir * plane_norm))
    # surface_normals = -np.matlib.repmat(plane_norm, numrays, 1) * np.sign(indir * np.transpose(plane_norm))
    crossing_into = np.round_(-np.sign(np.sum(indir * np.transpose(plane_norm), axis=1, keepdims=True)))

    return [intersection_points, surface_normals, distance_traveled, crossing_into]


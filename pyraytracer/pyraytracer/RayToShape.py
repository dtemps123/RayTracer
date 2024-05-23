#Called by IntersectFunction (defined in IntersectFunction.py) with the proper inputs
import numpy as np

## pyraytracer files
from .RayToCylinder import RayToCylinder
from .RayToPlane import RayToPlane
from .RayToSphere import RayToSphere
from .RayToTorus import RayToTorus
from .RayToQuadSurface import RayToQuadSurface
# Will have to import other geometries after writing them

def RayToShape(shape, sp, indir, param_list):

    # normalize directions
    goodray_cut = np.sum(indir ** 2, 1) > 0
    if np.any(goodray_cut):
        indir[goodray_cut, :] = indir[goodray_cut, :] / np.abs(np.sqrt(np.sum(indir[goodray_cut] ** 2, 1)))[:, np.newaxis]

    #maybe turn each of these into try-except statements, in case the specified
    # param_list isn't the right size of elements
    if shape == "cylinder":
        output = RayToCylinder(sp, indir, param_list[0], param_list[1], param_list[2])
    elif shape == "sphere":
        output = RayToSphere(sp, indir, param_list[0], param_list[1])
    elif shape == "torus":
        output = RayToTorus(sp, indir, param_list[0], param_list[1], param_list[2], param_list[3])
    elif shape == "plane":
        output = RayToPlane(sp, indir, param_list[0], param_list[1])
    elif shape == "quadsurface":
        output = RayToQuadSurface(sp, indir, param_list[0], param_list[1], param_list[2])
    else:
        raise Exception('Geometry has unrecognized shape')
        
    return output

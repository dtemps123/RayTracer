"""
.pyraytracer

Optical ray tracing code that propagates rays through a surface based geometry.
"""

__version__ = "0.1.0"
__author__ = 'C Eric Dahl & Dylan J Temples'
__credits__ = 'Fermi National Accelerator Laboratory'

try:
    from .rayInterfaces import *
    from .RayleighScatteringClass import *
    from .RefractionReflectionAtInterface import *
    from .UnifiedReflectorModel import *
    from .IntersectFunction import *
    from .surface import *
    from .RayTracer2 import *
    
    from .RayToShape import *
    from .RayToPlane import *
    from .RayToCylinder import *
    from .RayToSphere import *
    from .RayToTorus import *
    from .RayToQuadSurface import *
    
    from .RayTracer2_Display import *
    
except ImportError as err:
    print("\033[1;31mERROR\033[0m: Import error from pyraytracer lib.")
    print(err)

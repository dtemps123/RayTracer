import numpy as np
import math
import random
import RayTracer2
import surface


def ParabolicReflector(foc_z: float, n_in=10):
    surface_list = []
    """Make work both ways; start light at focus or start light incoming"""
    # Create a circular paraboloidal reflector along the z-axis with focus coordinate as given and vertex at origin
    # Focus is located at (0, 0, foc_z), foc_z = 1/(4a^2)
    # Equation of circular paraboloid: z = a(x^2 + y^2)
    a = .5 * math.sqrt(1 / foc_z)
    print("a: " + str(a))

    # First, create ray starting points and isotropic rays
    # Coordinates of startingpoint -- same for all rays
    x = 0
    y = 0

    n = n_in  # number of rays

    ray_startpoints = np.empty((n, 3))
    ray_startpoints[..., 0] = x
    ray_startpoints[..., 1] = y
    ray_startpoints[..., 2] = foc_z

    # rays = initial [forward direction (3), s1 polarization axis (3), s0, s1, s2, s3]
    # Let s0 be 1 and rays be unpolarized (other Stokes parameters = 0 + no polarization axis)
    test_rays = np.zeros((n, 10))
    test_rays[..., 3] = 1
    test_rays[..., 6] = 1

    # Assign initial forward directions in spherical coords for easier isotropism
    """These isotropic directions may be biased towards the origin, as they are not shifted for a non-origin starting point"""
    """
    for i in range(n):
        azimuth = random.random() * 2 * math.pi
        ran_a = np.random.uniform(-1, 1)
        polar = np.arccos(ran_a)  # -1 < theta < 1 so random numbers are not biased towards the pole

        test_rays[i, 0] = np.sin(polar) * np.cos(azimuth)  # x
        test_rays[i, 1] = np.sin(polar) * np.sin(azimuth)  # y
        test_rays[i, 2] = np.cos(polar)  # z
    """
    test_rays[0, 0] = 1
    test_rays[0, 1] = 0
    test_rays[0, 2] = -1
    test_rays[1, 0] = -1
    test_rays[1, 1] = 0
    test_rays[1, 2] = -1

    # Paraboloidal reflector consists of a single quadric surface
    reflector = surface.surface()
    reflector.description = "paraboloidal reflector with focus at (0, 0, " + str(foc_z) + ")"
    reflector.shape = 'quadsurface'
    reflector.param_list = [np.vstack(([a, 0, 0], [0, a, 0], [0, 0, 0])), np.array([0, 0, -1]), 0]
    reflector.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= 0) * (p[:, 2, :] < 10), (p.shape[0], -1)) # cut off paraboloid at z = 10cm
    reflector.n_outside = np.inf
    reflector.n_inside = 1
    reflector.surface_type = 'normal'
    reflector.absorption = 0
    surface_list.append(reflector)

    return ray_startpoints, test_rays, surface_list

def donut(n_in):
    surface_list = []

    center = np.array([0, 0, 0])
    ax = np.array([0, 0, 1])
    r1 = 5
    r2 = 8

    # First, create ray starting points and isotropic rays
    # Coordinates of startingpoint -- same for all rays
    x = 6
    y = 0
    z = 1

    n = n_in  # number of rays

    ray_startpoints = np.empty((n, 3))
    ray_startpoints[..., 0] = x
    ray_startpoints[..., 1] = y
    ray_startpoints[..., 2] = z

    # rays = initial [forward direction (3), s1 polarization axis (3), s0, s1, s2, s3]
    # Let s0 be 1 and rays be unpolarized (other Stokes parameters = 0 + no polarization axis)
    test_rays = np.zeros((n, 10))
    test_rays[..., 3] = 1
    test_rays[..., 6] = 1

    # Assign initial forward directions in spherical coords for easier isotropism
    """These isotropic directions may be biased towards the origin, as they are not shifted for a non-origin starting point"""
    for i in range(n):
        azimuth = random.random() * 2 * math.pi
        a = np.random.uniform(-1, 1)
        polar = np.arccos(a)  # -1 < theta < 1 so random numbers are not biased towards the pole

        test_rays[i, 0] = np.sin(polar) * np.cos(azimuth)  # x
        test_rays[i, 1] = np.sin(polar) * np.sin(azimuth)  # y
        test_rays[i, 2] = np.cos(polar)  # z

    # Donut consists of a single torus and outer sphere for collecting escaped rays
    reflector = surface.surface()
    reflector.description = "fully enclosed donut aligned on z-axis with center at origin"
    reflector.shape = 'torus'
    reflector.param_list = [center, ax, r1, r2]
    reflector.inbounds_function = lambda p: np.reshape(np.abs(p[:,2,:]) <= (r2 - r1),(p.shape[0], -1))
    reflector.n_outside = np.inf
    reflector.n_inside = 1
    reflector.surface_type = 'normal'
    reflector.absorption = 0
    surface_list.append(reflector)

    collector = surface.surface()
    collector.description = "spherical surface for checking for escaped rays"
    collector.shape = 'sphere'
    collector.param_list = [np.array([0, 0, 0]), 20]
    collector.inbounds_function = lambda p: np.reshape((p[:,0,:]**2 + p[:,1,:]**2 + p[:,2,:]**2) == 20, (p.shape[0], -1))
    collector.n_outside = np.inf
    collector.n_inside = np.inf
    collector.surface_type = 'normal'
    collector.absorption = 1
    surface_list.append(collector)

    return ray_startpoints, test_rays, surface_list

def main():
    """
    z = 3   # height of focus
    n = 2  # number of rays

    [starts, rays, surfaces] = ParabolicReflector(z, n)
    [ray_interfaces, absorption_table, raytable] = RayTracer2.RayTracer2(starts, rays, surfaces)

    print("Height of focus: " + str(z))

    for i in range(len(ray_interfaces)):
        print("Direction of incoming " + str(i+1) + ": " + str(ray_interfaces[i].incoming_ray))
        print("Direction of reflected " + str(i + 1) + ": " + str(ray_interfaces[i].reflected_ray))
    """

    n = 4 # number of rays

    [starts, rays, surfaces] = donut(n)
    [ray_interfaces, absorption_table, raytable] = RayTracer2.RayTracer2(starts, rays, surfaces)

    print("Any escaped rays?")
    print(absorption_table)

if __name__ == "__main__":
    main()


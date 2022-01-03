# SBC Geometry from the module itself plus some added shapes

import numpy as np
import math
import random
import RayTracer2
import surface
import matplotlib.pyplot as plt

void = np.inf
n_target = 1.17 # n=1.224 for Ar @ 940nm, 90K
n_jar = 1.4512 # SiO2 @ 940nm
n_hydraulic = 1.22 # 1.21 @ CF4; 940nm, 146K ;  1.237 @ 7eV, 146K, another said 1.515
n_pressurewindow = 1.7569 # Al2O3 @ 940nm
n_pressurewall = void
n_air = 1.00



## Useful equations
# n0 = index of refraction at known density
# r = rho/rho0 (rho = density now, rho0 = density corresponding to n0)
# returns index of refraction at density rho
# clausius_mossotti = @(n0, r)(sqrt(((1 + 2*r).*n0.*n0 + 2 - 2*r)./((1 - r).*n0.*n0 + 2 + r)));
clausius_mossotti = lambda n0, r: math.sqrt(((1 + 2*r) * n0**2 + 2 - 2*r) / ((1 - r) * n0**2 + 2 + r))

## Dimensions in cm
ojar_thick = .5 # thickness of cylinder wall
ojar_cylrad = 12 # outer radius of cylinder
ojar_axrad = 24 # outer radius of sphere (along cylinder axis)
ojar_knucklerad = 4
ojar_cyllength = 30
ojar_elevation = 0

ijar_thick = .5 # thickness of cylinder wall
ijar_cylrad = 10.5 # outer radius of cylinder
ijar_axrad = 21 # outer radius of sphere (along cylinder axis)
ijar_knucklerad = 3
ijar_cyllength = 10.0
ijar_elevation = -2.54*7.74

z0 = 2.54*(20.73+7.74) # +distance from seal to z=0
vp_focuselev = 65.40 - z0
vp_focuslen = 29.23 #25;%33.5; % distance from convergence point to CF seal
vp_phi = np.pi/3

### the below variables (vp_{variable}) are for an outdated viewport ###
### Updated viewport variables are listed with the viewport surface code ###

vp_s = 10.0 # radial position of air-side center of viewport
vp_elev = 60.0 # vertical position of air-side center of viewport
vp_win_rad = 1.375*.5*2.54 # radius of glass (black on circumference)
vp_air_rad = 1.25*.5*2.54 # radius of air-side can (black on circumference)
vp_can_rad = 2.54
vp_can_wall = .0625*2.54
vp_flange_rad = .5*5.5*2.54
vp_nip_rad = 2.54*4.75*.5 # radius of hydraulic-side nipple (black on circumference)
vp_win_thick = .2*2.54 #%mpf is .23*2.54;
vp_nip_top = -.254*2.54 #% mpf is .08*2.54% 
vp_theta = .3737 #20*pi/180;%21.8*pi/180;
vp_can_OAL = 6.43*2.54
vp_flange_thick = np.array([2.88, .69, .625, .625, .625])* 2.54


rdcone_gap = 1.5 # this is used to estimate the amount of light that may get through the top cone and the rdcone, since Solidworks shows
                   # a decent seperation between those general shapes.  This is a pretty rough estimate, we'll play with it to see how
                    # much of an effect it has.
rd_rad = 12.5 # reflector - diffuser radius
rd_top = 0
rd_bot = -25
rdcone_top = 8.96
rdcone_toprad = 15.13
rdtopcone_rad = rdcone_toprad
rdtopcone_bot = rdcone_top + rdcone_gap
rdtopcone_apex = rdtopcone_bot+6.5
rdbotcone_apex = -15.2-10
rdbotcone_rad = 10.5
rdbotcone_bot = -20.0-10


pv_top = vp_focuselev + vp_focuslen*np.cos(vp_theta)-2.54*(5.49-1.5)
pv_bot = pv_top - 2.54*31.17
pv_rad = 8*2.54
pv_thick = .375*2.54
pv_axrad = (3.07)*2.54

## Derived dimensions
t_o = np.array([0, ojar_thick])
t_i = np.array([0, ijar_thick])

r1 = np.concatenate((ojar_cylrad - t_o, ijar_cylrad - t_i))
r2 = np.concatenate((ojar_knucklerad - t_o, ijar_knucklerad - t_i))
r3 = np.concatenate((ojar_axrad - t_o, ijar_axrad - t_i))

s = r3 * (r1-r2) / (r3-r2) # axis to knuckle-dome transition

z = r2 * np.sqrt(1 - (s / r3) ** 2) #  equator to knuckle-dome transition

d = r3 * z * ((1 / r3) - (1 / r2)) # equator to dome sphere center

vp_axis = np.array([0, -math.sin(vp_theta), math.cos(vp_theta)])
vp_center = np.array([0, -vp_s, vp_elev])


head_out_Q = np.vstack(([pv_rad ** (-2), 0, 0], [0, pv_rad ** (-2), 0],  [0, 0, pv_axrad ** (-2)]))
head_in_Q = np.vstack(([pv_rad-pv_thick ** (-2), 0, 0], [0, pv_rad-pv_thick ** (-2), 0], [0, 0, pv_axrad-pv_thick ** (-2)]))
head_out_P = np.array([0, 0, ((-2)*pv_top) * (pv_axrad ** (-2))]) # .* is element-wise, * is matrix multiplication
head_in_P = np.array([0, 0, ((-2)*pv_top) * (pv_axrad-pv_thick ** (-2))])
head_out_R = (pv_top / pv_axrad)**2 - 1
head_in_R = (pv_top / (pv_axrad - pv_thick))**2 - 1

rd_cone_b = (rdcone_toprad - rd_rad) / (rdcone_top - rd_top)
rd_cone_z0 = rd_top - (rd_rad / rd_cone_b)
rd_cone_Q = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, -rd_cone_b**2])) # - (rd_cone_b ** 2) or (-rd_cone_b) ** 2  ?
rd_cone_P = np.array([0, 0, 2* rd_cone_b**2 * rd_cone_z0])
rd_cone_R = -(rd_cone_b * rd_cone_z0)**2

rd_topcone_b = rdtopcone_rad / (rdtopcone_apex - rdtopcone_bot)
rd_topcone_Q = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, -rd_topcone_b**2]))
rd_topcone_P = np.array([0, 0, 2 * rd_topcone_b**2 * rdtopcone_apex])
rd_topcone_R = -(rd_topcone_b * rdtopcone_apex)**2

rd_botcone_b = rdbotcone_rad / (rdbotcone_apex - rdbotcone_bot)
rd_botcone_Q = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, -rd_botcone_b**2]))
rd_botcone_P = np.array([0, 0, 2 * rd_botcone_b**2 * rdbotcone_apex])
rd_botcone_R = -(rd_botcone_b * rdbotcone_apex)**2

## The surfaces themselves

## Useful equations
# n0 = index of refraction at known density
# r = rho/rho0 (rho = density now, rho0 = density corresponding to n0)
# returns index of refraction at density rho
# clausius_mossotti = @(n0, r)(sqrt(((1 + 2*r).*n0.*n0 + 2 - 2*r)./((1 - r).*n0.*n0 + 2 + r)));
clausius_mossotti = lambda n0, r: math.sqrt(((1 + 2*r) * n0**2 + 2 - 2*r) / ((1 - r) * n0**2 + 2 + r))

surface_list = []

## Inner Jar
#inner_jar_density = -2.203

inner_jar_inner_cyl = surface.surface()
inner_jar_inner_cyl.description = 'inside surface of inner quartz jar cylinder'
inner_jar_inner_cyl.shape = 'cylinder'
inner_jar_inner_cyl.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), r1[3]] # shouldn't r be 10?
inner_jar_inner_cyl.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= (ijar_elevation-ijar_cyllength)) * (p[:, 2, :] < ijar_elevation), (p.shape[0], -1)) # -66.40278, -21.61143
inner_jar_inner_cyl.n_outside = n_jar
inner_jar_inner_cyl.n_inside = n_hydraulic
inner_jar_inner_cyl.surface_type = 'normal'
inner_jar_inner_cyl.absorption = 0
surface_list.append(inner_jar_inner_cyl)

inner_jar_outer_cyl = surface.surface()
inner_jar_outer_cyl.description = 'outside surface of inner jar cylinder'
inner_jar_outer_cyl.shape = 'cylinder'
inner_jar_outer_cyl.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), r1[2]] # 10.5
inner_jar_outer_cyl.inbounds_function = lambda p: np.reshape((p[:,2,:] >= (ijar_elevation-ijar_cyllength)) * (p[:,2,:] < ijar_elevation), (p.shape[0], -1)) # -66.40278, -21.61143
inner_jar_outer_cyl.n_outside = n_target
inner_jar_outer_cyl.n_inside = n_jar
inner_jar_outer_cyl.surface_type = 'normal'
inner_jar_outer_cyl.absorption = 0
surface_list.append(inner_jar_outer_cyl)

# will have to split up above into different sections with insulation rather than CF4 as the previous/next material (?)

inner_jar_inner_dome = surface.surface()
inner_jar_inner_dome.description = 'inside surface of inner jar dome'
inner_jar_inner_dome.shape = 'sphere'
inner_jar_inner_dome.param_list = [np.array([0, 0, ijar_elevation + d[3]]), r3[3]] # -37.9745, 20.5
inner_jar_inner_dome.inbounds_function = lambda p: np.reshape((p[:,2,:] > (z[3] + ijar_elevation)), (p.shape[0], -1)) # -21.61143
inner_jar_inner_dome.n_outside = n_jar
inner_jar_inner_dome.n_inside = n_hydraulic
inner_jar_inner_dome.surface_type = 'normal'
inner_jar_inner_dome.absorption = 0
surface_list.append(inner_jar_inner_dome)

inner_jar_outer_dome = surface.surface()
inner_jar_outer_dome.description = 'outside surface of inner jar dome'
inner_jar_outer_dome.shape = 'sphere'
inner_jar_outer_dome.param_list = [np.array([0, 0, ijar_elevation + d[2]]), r3[2]] # -37.9745, 21.0
inner_jar_outer_dome.inbounds_function = lambda p: np.reshape((p[:,2,:] > (z[2] + ijar_elevation)), (p.shape[0], -1)) # -21.61143
inner_jar_outer_dome.n_outside = n_target
inner_jar_outer_dome.n_inside = n_jar
inner_jar_outer_dome.surface_type = 'normal'
inner_jar_outer_dome.absorption = 0
surface_list.append(inner_jar_outer_dome)


## Outer Jar
outer_jar_inner_cyl = surface.surface()
outer_jar_inner_cyl.description = 'inner surface of outer jar cylinder'
outer_jar_inner_cyl.shape = 'cylinder'
outer_jar_inner_cyl.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), r1[1]] # 11.5
outer_jar_inner_cyl.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= (ojar_elevation-ojar_cyllength)) * (p[:, 2, :] < ojar_elevation), (p.shape[0], -1)) # -66.40278, 0
outer_jar_inner_cyl.n_outside = n_jar
outer_jar_inner_cyl.n_inside = n_target
outer_jar_inner_cyl.surface_type = 'normal'
outer_jar_inner_cyl.absorption = 0
surface_list.append(outer_jar_inner_cyl)

outer_jar_outer_cyl = surface.surface()
outer_jar_outer_cyl.description = 'outside surface of outer jar cylinder'
outer_jar_outer_cyl.shape = 'cylinder'
outer_jar_outer_cyl.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), r1[0]] # 12
outer_jar_outer_cyl.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= (ojar_elevation-ojar_cyllength)) * (p[:, 2, :] < ojar_elevation), (p.shape[0], -1)) # -66.40278, 0
outer_jar_outer_cyl.n_outside = n_hydraulic
outer_jar_outer_cyl.n_inside = n_jar
outer_jar_outer_cyl.surface_type = 'normal'
outer_jar_outer_cyl.absorption = 0
surface_list.append(outer_jar_outer_cyl)

outer_jar_inner_dome = surface.surface()
outer_jar_inner_dome.description = 'inner surface of outer jar dome'
outer_jar_inner_dome.shape = 'sphere'
outer_jar_inner_dome.param_list = [np.array([0, 0, ojar_elevation + d[1]]), r3[1]] # -18.3303, 23.5
outer_jar_inner_dome.inbounds_function = lambda p: np.reshape((p[:,2,:] > 0), (p.shape[0], -1)) # < needed too ?
outer_jar_inner_dome.n_outside = n_jar
outer_jar_inner_dome.n_inside = n_target
outer_jar_inner_dome.surface_type = 'normal'
outer_jar_inner_dome.absorption = 0
surface_list.append(outer_jar_inner_dome)

outer_jar_outer_dome = surface.surface()
outer_jar_outer_dome.description = 'outer surface of outer jar domee'
outer_jar_outer_dome.shape = 'sphere'
outer_jar_outer_dome.param_list = [np.array([0, 0, ojar_elevation + d[0]]), r3[1]] # -18.3303, 24
outer_jar_outer_dome.inbounds_function = lambda p: np.reshape((p[:,2,:] > (z[0] + ojar_elevation)), (p.shape[0], -1)) # 0
outer_jar_outer_dome.n_outside = n_hydraulic
outer_jar_outer_dome.n_inside = n_jar
outer_jar_outer_dome.surface_type = 'normal'
outer_jar_outer_dome.absorption = 0
surface_list.append(outer_jar_outer_dome)

## Knuckles
iiknuckle = surface.surface()
iiknuckle.description = 'inner surface of inner knuckle'
iiknuckle.shape = 'torus'
iiknuckle.param_list = [np.array([0, 0, ijar_elevation]), np.array([0, 0, 1]), r1[3] - r2[3], r2[3]]
iiknuckle.inbounds_function = lambda p: np.reshape((p[:,2,:] > ijar_elevation) * (p[:,2,:] <= (z[3]+ijar_elevation)) *((p[:,0,:]**2 + p[:,1,:]** 2) > ((r1[3] - r2[3])**2)), (p.shape[0], -1))
iiknuckle.n_outside = n_jar
iiknuckle.n_inside = n_hydraulic
iiknuckle.surface_type = 'normal'
iiknuckle.absorption = 0
surface_list.append(iiknuckle)

oiknuckle = surface.surface()
oiknuckle.description = 'outer surface of inner knuckle'
oiknuckle.shape = 'torus'
oiknuckle.param_list = [np.array([0, 0, ijar_elevation]), np.array([0, 0, 1]), r1[2] - r2[2], r2[3]]
oiknuckle.inbounds_function = lambda p: np.reshape((p[:,2,:] > ijar_elevation) * (p[:,2,:] <= (z[2] + ijar_elevation))
                                                   * ((p[:,0,:]**2 + p[:,1,:]**2) > ((r1[2] - r2[2])**2)), (p.shape[0], -1))
oiknuckle.n_outside = n_hydraulic
oiknuckle.n_inside = n_jar
oiknuckle.surface_type = 'normal'
oiknuckle.absorption = 0
surface_list.append(oiknuckle)

ioknuckle = surface.surface()
ioknuckle.description = 'inner surface of outer knuckle'
ioknuckle.shape = 'torus'
ioknuckle.param_list = [np.array([0, 0, ojar_elevation]), np.array([0, 0, 1]), r1[1] - r2[1], r2[1]]
ioknuckle.inbounds_function = lambda p: np.reshape((p[:,2,:] > ojar_elevation) * (p[:,2,:] <= (z[1] + ojar_elevation))
                                                   * ((p[:,0,:]**2 + p[:,1,:]**2) > ((r1[1] - r2[1])**2)), (p.shape[0], -1))
ioknuckle.n_outside = n_jar
ioknuckle.n_inside = n_target
ioknuckle.surface_type = 'normal'
ioknuckle.absorption = 0
surface_list.append(ioknuckle)

ooknuckle = surface.surface()
ooknuckle.description = 'outer surface of outer knuckle'
ooknuckle.shape = 'torus'
ooknuckle.param_list = [np.array([0, 0, ojar_elevation]), np.array([0, 0, 1]), r1[0] - r2[0], r2[0]]
ooknuckle.inbounds_function = lambda p: np.reshape((p[:,2,:] > ojar_elevation) * (p[:,2,:] <= (z[0] + ojar_elevation))
                                                   * ((p[:,0,:]**2 + p[:,1,:]**2) > ((r1[0] - r2[0])**2)), (p.shape[0], -1))
ooknuckle.n_outside = n_hydraulic
ooknuckle.n_inside = n_jar
ooknuckle.surface_type = 'normal'
ooknuckle.absorption = 0
surface_list.append(ooknuckle)

## other black surfaces to trap rays
rd = surface.surface()
rd.description = 'reflector/diffuser'
rd.shape = 'cylinder'
rd.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
rd.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_bot) * (p[:, 2, :] <= rd_top), (p.shape[0], -1)) \
& ~(SiPM_1.inbounds_function(p))& ~(SiPM_2.inbounds_function(p))& ~(SiPM_3.inbounds_function(p))\
& ~(SiPM_4.inbounds_function(p))& ~(SiPM_5.inbounds_function(p))& ~(SiPM_6.inbounds_function(p))\
& ~(SiPM_7.inbounds_function(p))& ~(SiPM_8.inbounds_function(p))& ~(SiPM_9.inbounds_function(p))\
& ~(SiPM_10.inbounds_function(p))& ~(SiPM_11.inbounds_function(p))& ~(SiPM_12.inbounds_function(p))\
& ~(SiPM_13.inbounds_function(p))& ~(SiPM_14.inbounds_function(p))& ~(SiPM_15.inbounds_function(p))\
& ~(SiPM_16.inbounds_function(p))& ~(SiPM_17.inbounds_function(p))& ~(SiPM_18.inbounds_function(p))\
& ~(SiPM_19.inbounds_function(p))& ~(SiPM_20.inbounds_function(p))& ~(SiPM_21.inbounds_function(p))\
& ~(SiPM_22.inbounds_function(p))& ~(SiPM_23.inbounds_function(p))& ~(SiPM_24.inbounds_function(p))\
& ~(SiPM_25.inbounds_function(p))& ~(SiPM_26.inbounds_function(p))& ~(SiPM_27.inbounds_function(p))\
& ~(SiPM_28.inbounds_function(p))& ~(SiPM_29.inbounds_function(p))& ~(SiPM_30.inbounds_function(p))\
& ~(SiPM_31.inbounds_function(p))& ~(SiPM_32.inbounds_function(p))
rd.n_outside = n_hydraulic
rd.n_inside = n_hydraulic
rd.surface_type = 'normal'
rd.absorption = 1
surface_list.append(rd)

rd_cone = surface.surface()
rd_cone.description = 'reflector/diffuser cone'
rd_cone.shape = 'quadsurface'                       # QuadSurface needs to be tested
rd_cone.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
rd_cone.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top), (p.shape[0], -1)) \
& ~(SiPM_33.inbounds_function(p))& ~(SiPM_34.inbounds_function(p))& ~(SiPM_35.inbounds_function(p))\
& ~(SiPM_36.inbounds_function(p))& ~(SiPM_37.inbounds_function(p))& ~(SiPM_38.inbounds_function(p))\
& ~(SiPM_39.inbounds_function(p))& ~(SiPM_40.inbounds_function(p))
rd_cone.n_outside = n_hydraulic
rd_cone.n_inside = n_hydraulic
rd_cone.surface_type = 'normal'
rd_cone.absorption = 1
surface_list.append(rd_cone)

rd_topcone = surface.surface()
rd_topcone.description = 'reflector/diffuser topcone'
rd_topcone.shape = 'quadsurface'                       # QuadSurface needs to be tested
rd_topcone.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
rd_topcone.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex)\
& ~(viewport_1.inbounds_function(p)) & ~(viewport_2.inbounds_function(p))\
& ~(viewport_3.inbounds_function(p)) & ~(viewport_flange_1.inbounds_function(p))\
& ~(viewport_flange_2.inbounds_function(p)) & ~(viewport_flange_3.inbounds_function(p)), (p.shape[0], -1))
rd_topcone.n_outside = n_hydraulic
rd_topcone.n_inside = n_hydraulic
rd_topcone.surface_type = 'normal'
rd_topcone.absorption = 1
surface_list.append(rd_topcone)
 
### The viewports for the cameras ###
## Both the viewports' and their flanges are just inbounds functions that pick out pieces of the rd topcone.
#useful viewport variables:
vp_radius = 6.25                # viewport radius
vp_flange_radius = 7            # viewport flange radius (annulus bordering the viewport, mounts the veiwports to the upper cone)
topcone_slope = (rdtopcone_apex - rdtopcone_bot)/(rdcone_toprad)   # the slope of the cone's sides
topcone_theta = np.arctan(topcone_slope) #the angle made between the xy plane and the top cone's sides.

viewport_1 = surface.surface()
viewport_1.description = 'camera viewport on top cone on +y side'
viewport_1.shape = 'quadsurface'                       # QuadSurface needs to be tested
viewport_1.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
viewport_1.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex)\
& ((p[:,0,:]**2+((p[:,1,:]-0.51*rdcone_toprad)/np.cos(topcone_theta))**2)<vp_radius**2),(p.shape[0], -1))
viewport_1.n_outside = n_hydraulic
viewport_1.n_inside = n_hydraulic
viewport_1.surface_type = 'normal'
viewport_1.absorption = 0
surface_list.append(viewport_1)

viewport_2 = surface.surface()
viewport_2.description = 'camera viewport on top cone on -y, +x side'
viewport_2.shape = 'quadsurface'                       # QuadSurface needs to be tested
viewport_2.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
viewport_2.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex)\
& ((((p[:,0,:]-0.51*rdcone_toprad*np.cos(7*np.pi/6))*np.cos(4*np.pi/6)+(p[:,1,:]-0.51*rdcone_toprad*np.sin(7*np.pi/6))\
*np.sin(4*np.pi/6))**2+((((p[:,0,:]-0.51*rdcone_toprad*np.cos(7*np.pi/6))*np.sin(4*np.pi/6))-(p[:,1,:]-0.51*rdcone_toprad\
*np.sin(7*np.pi/6))*np.cos(4*np.pi/6))/np.cos(topcone_theta))**2)<vp_radius**2),(p.shape[0], -1))
viewport_2.n_outside = n_hydraulic
viewport_2.n_inside = n_hydraulic
viewport_2.surface_type = 'normal'
viewport_2.absorption = 0
surface_list.append(viewport_2)

viewport_3 = surface.surface()
viewport_3.description = 'camera viewport on top cone on -y, -x side'
viewport_3.shape = 'quadsurface'                       # QuadSurface needs to be tested
viewport_3.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
viewport_3.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex)\
& ((((p[:,0,:]-0.51*rdcone_toprad*np.cos(-np.pi/6))*np.cos(-4*np.pi/6)+(p[:,1,:]-0.51*rdcone_toprad*np.sin(-np.pi/6))\
*np.sin(-4*np.pi/6))**2+((((p[:,0,:]-0.51*rdcone_toprad*np.cos(-np.pi/6))*np.sin(-4*np.pi/6))-(p[:,1,:]-0.51*rdcone_toprad\
*np.sin(-np.pi/6))*np.cos(-4*np.pi/6))/np.cos(topcone_theta))**2)<vp_radius**2),(p.shape[0], -1))
viewport_3.n_outside = n_hydraulic
viewport_3.n_inside = n_hydraulic
viewport_3.surface_type = 'normal'
viewport_3.absorption = 0
surface_list.append(viewport_2)

#viewport flanges

#flange 1
viewport_flange_1 = surface.surface()
viewport_flange_1.description = 'flange for viewport 1 on top cone on +y side'
viewport_flange_1.shape = 'quadsurface'                       # QuadSurface needs to be tested
viewport_flange_1.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
viewport_flange_1.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex)\
& ((p[:,0,:]**2+((p[:,1,:]-0.51*rdcone_toprad)/np.cos(topcone_theta))**2)<vp_flange_radius**2)\
& ~(viewport_1.inbounds_function(p)),(p.shape[0], -1))
viewport_flange_1.n_outside = n_hydraulic
viewport_flange_1.n_inside = n_hydraulic
viewport_flange_1.surface_type = 'normal'
viewport_flange_1.absorption = 1
surface_list.append(viewport_flange_1)

#flange 2
viewport_flange_2 = surface.surface()
viewport_flange_2.description = 'flange for viewport 2 on top cone on -y, +x side'
viewport_flange_2.shape = 'quadsurface'                       # QuadSurface needs to be tested
viewport_flange_2.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
viewport_flange_2.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex)\
& ((((p[:,0,:]-0.51*rdcone_toprad*np.cos(7*np.pi/6))*np.cos(4*np.pi/6)+(p[:,1,:]-0.51*rdcone_toprad*np.sin(7*np.pi/6))\
*np.sin(4*np.pi/6))**2+((((p[:,0,:]-0.51*rdcone_toprad*np.cos(7*np.pi/6))*np.sin(4*np.pi/6))-(p[:,1,:]-0.51*rdcone_toprad\
*np.sin(7*np.pi/6))*np.cos(4*np.pi/6))/np.cos(topcone_theta))**2)<vp_flange_radius**2)\
& ~(viewport_2.inbounds_function(p)),(p.shape[0], -1))
viewport_flange_2.n_outside = n_hydraulic
viewport_flange_2.n_inside = n_hydraulic
viewport_flange_2.surface_type = 'normal'
viewport_flange_2.absorption = 1
surface_list.append(viewport_flange_2)

#flange 3
viewport_flange_3 = surface.surface()
viewport_flange_3.description = 'flange for viewport 3 on top cone on -y, -x side'
viewport_flange_3.shape = 'quadsurface'                       # QuadSurface needs to be tested
viewport_flange_3.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
viewport_flange_3.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex)\
& ((((p[:,0,:]-0.51*rdcone_toprad*np.cos(-np.pi/6))*np.cos(-4*np.pi/6)+(p[:,1,:]-0.51*rdcone_toprad*np.sin(-np.pi/6))\
*np.sin(-4*np.pi/6))**2+((((p[:,0,:]-0.51*rdcone_toprad*np.cos(-np.pi/6))*np.sin(-4*np.pi/6))-(p[:,1,:]-0.51*rdcone_toprad\
*np.sin(-np.pi/6))*np.cos(-4*np.pi/6))/np.cos(topcone_theta))**2)<vp_flange_radius**2)\
& ~(viewport_3.inbounds_function(p)),(p.shape[0], -1))
viewport_flange_3.n_outside = n_hydraulic
viewport_flange_3.n_inside = n_hydraulic
viewport_flange_3.surface_type = 'normal'
viewport_flange_3.absorption = 1
surface_list.append(viewport_flange_3)

#camera cans
# useful variables:

#this is the z value for all of the center points of the viewports.  This allows us to use those points for our cylinder axis.  
zpoint = -topcone_slope*0.51*rdcone_toprad+rdtopcone_apex 

camcan_h = 7  #height of the camera cans [cm]

# this is the upper z point if you follow the cylinder axis to the top end of the camera can. Important for inbounds function.
zpoint_upper = -topcone_slope*0.51*rdcone_toprad+rdtopcone_apex+(camcan_h*np.cos(topcone_theta))

camcan_1 = surface.surface()
camcan_1.description = 'can leading to camera on the +y side of the detector'
camcan_1.shape = 'cylinder'
camcan_1.param_list = [np.array([0, 0.51*rdcone_toprad, zpoint]), np.array([0, 1, np.cos(topcone_theta)]), vp_radius] 
camcan_1.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= (p[:,1,:]-0.51*rdcone_toprad)*-topcone_slope+zpoint) \
& (p[:, 2, :] <= (p[:,1,:]-(0.51*rdcone_toprad+camcan_h*np.sin(topcone_theta)))*-topcone_slope+zpoint_upper), (p.shape[0], -1))
camcan_1.n_outside = n_hydraulic
camcan_1.n_inside = n_hydraulic
camcan_1.surface_type = 'normal'
camcan_1.absorption = 1
surface_list.append(camcan_1)

camcan_2 = surface.surface()
camcan_2.description = 'can leading to camera on the -x, -y side of the detector'
camcan_2.shape = 'cylinder'
camcan_2.param_list = [np.array([0.51*rdcone_toprad*np.cos(7*np.pi/6), 0.51*rdcone_toprad*np.sin(7*np.pi/6), zpoint]),\
                       np.array([np.cos(7*np.pi/6), np.sin(7*np.pi/6), np.cos(topcone_theta)]), vp_radius] 
camcan_2.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= (p[:,0,:]-0.51*rdcone_toprad*np.cos(7*np.pi/6))\
*topcone_slope+zpoint) & (p[:, 2, :] <= (p[:,0,:]-(0.51*rdcone_toprad*np.cos(7*np.pi/6)-camcan_h*np.sin(topcone_theta)))\
*topcone_slope+zpoint_upper), (p.shape[0], -1))
camcan_2.n_outside = n_hydraulic
camcan_2.n_inside = n_hydraulic
camcan_2.surface_type = 'normal'
camcan_2.absorption = 1
surface_list.append(camcan_2)

camcan_3 = surface.surface()
camcan_3.description = 'can leading to camera on the -x, -y side of the detector'
camcan_3.shape = 'cylinder'
camcan_3.param_list = [np.array([0.51*rdcone_toprad*np.cos(-np.pi/6), 0.51*rdcone_toprad*np.sin(-np.pi/6), zpoint]),\
                       np.array([np.cos(-np.pi/6), np.sin(-np.pi/6), np.cos(topcone_theta)]), vp_radius] 
camcan_3.inbounds_function = lambda p: np.reshape((p[:, 2, :] >= (p[:,0,:]-0.51*rdcone_toprad*np.cos(-np.pi/6))\
*-topcone_slope+zpoint) & (p[:, 2, :] <= (p[:,0,:]-(0.51*rdcone_toprad*np.cos(-np.pi/6)-camcan_h*np.sin(topcone_theta)))\
*-topcone_slope+zpoint_upper), (p.shape[0], -1))
camcan_3.n_outside = n_hydraulic
camcan_3.n_inside = n_hydraulic
camcan_3.surface_type = 'normal'
camcan_3.absorption = 1
surface_list.append(camcan_3)

# upper most cone (only using the part of the cone that caps off the camera cans)
#useful variables for the camera end caps:
camcancaps_rad = rdtopcone_rad+(camcan_h*np.sin(topcone_theta))
camcancaps_apex = rdtopcone_apex + (camcan_h*np.cos(topcone_theta))
cancamcaps_bot = camcancaps_apex - (rdtopcone_apex-rdtopcone_bot)

camcancaps_b = camcancaps_rad / (camcancaps_apex - cancamcaps_bot)
camcancaps_Q = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, -camcancaps_b**2]))
camcancaps_P = np.array([0, 0, 2 * camcancaps_b**2 * camcancaps_apex])
camcancaps_R = -(camcancaps_b * camcancaps_apex)**2

camcancaps = surface.surface()
camcancaps.description = 'end caps for the camera cans'
camcancaps.shape = 'quadsurface'                       # QuadSurface needs to be tested
camcancaps.param_list = [camcancaps_Q, camcancaps_P, camcancaps_R]
camcancaps.inbounds_function = lambda p: np.reshape((p[:, 2, :] > cancamcaps_bot) * (p[:, 2, :] < camcancaps_apex)\
                                , (p.shape[0], -1))
camcancaps.n_outside = n_hydraulic
camcancaps.n_inside = n_hydraulic
camcancaps.surface_type = 'normal'
camcancaps.absorption = 1
surface_list.append(camcancaps)


#SiPMs!!
# Cuts of the cylinder that represent the SiPMs

# Important NOTE: in SBCGeometry: rd_cone is the cone extending outward that has those upward pointed SiPMs.
# add these variables to the top of SBC Geometry:
# Assumes an octagonal setup of the SiPMs (unlikely to change that part of the SBC design)

SiPM_d = 1.269 #dimensions of the square SiPM surface

SiPM_r1 = [-17.88,-16.61]  #SiPMs in the first (bottom most) row's z: bottom z at 0 index, top value at index 1
SiPM_r2 = [-12.71,-11.44]  #SiPMs in the second row's z: bottom z at 0 index, top value at index 1
SiPM_r3 = [-7.54,-6.27]  #SiPMs in the third row's z: bottom z at 0 index, top value at index 1
SiPM_r4 = [-2.47,-1.2]  #SiPMs in the fourth row's z: bottom z at 0 index, top value at index 1
#SiPM_r5 = [b,t]  #SiPMs in the fifth row's (on the rd_cone) z: bottom z at 0 index, top value at index 1

# The convention for SiPM columns is that right is towards positive x and left is towards negative x
### NOTE: columns 1 and 5 use y values, not x. Makes it easier to calculate ###

# column 1 values use y bound
SiPM_c1y = [-1/2*SiPM_d,1/2*SiPM_d] #SiPMs in first column (right most) x: left x at 0 index, right x at 1 index (infinitely thin)
#SiPMs in second column +x angled side x: left x at 0 index, right x at 1 index
SiPM_c2 = [(1/2*(2*rd_rad/(1+np.sqrt(2)))+2*rd_rad/(2*np.sqrt(2)*(1+np.sqrt(2))))-1/2*SiPM_d/np.sqrt(2)\
          ,(1/2*(2*rd_rad/(1+np.sqrt(2)))+2*rd_rad/(2*np.sqrt(2)*(1+np.sqrt(2))))+1/2*SiPM_d/np.sqrt(2)] #only if octagonal
SiPM_c3 = [-1/2*SiPM_d,1/2*SiPM_d] #SiPMs in third column (right most) x: left x at 0 index, right x at 1 index
 #SiPMs in fourth column -x angled side x: left x at 0 index, right x at 1 index
SiPM_c4 = [-(1/2*(2*rd_rad/(1+np.sqrt(2)))+2*rd_rad/(2*np.sqrt(2)*(1+np.sqrt(2))))-1/2*SiPM_d/np.sqrt(2)\
          ,-(1/2*(2*rd_rad/(1+np.sqrt(2)))+2*rd_rad/(2*np.sqrt(2)*(1+np.sqrt(2))))+1/2*SiPM_d/np.sqrt(2)] #only if octagonal
# column 5 values use y bounds 
SiPM_c5y = [-1/2*SiPM_d,1/2*SiPM_d]  #SiPMs in fifth column (left most) x: left x at 0 index, right x at 1 index (infinitely thin)
# heavily based on rd which is the surface they are mounted on.

## Numbering convention for SiPMs:  Start at 1 with the bottom most SiPM in the positive x direction. SiPM 2 continues
## along that same z height row counterclockwise.  Once you get to SiPM 9 you will be starting the next row.

## SiPMs in First Row

# SiPM 1
SiPM_1 = surface.surface()
SiPM_1.description = 'SiPM bottom row at +x, y=0'
SiPM_1.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_1.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_1.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,0,:]>0)), (p.shape[0],-1))
SiPM_1.n_outside = n_hydraulic
SiPM_1.n_inside = n_hydraulic
SiPM_1.surface_type = 'normal'
SiPM_1.absorption = 1
surface_list.append(SiPM_1)

#SiPM 2
SiPM_2 = surface.surface()
SiPM_2.description = 'SiPM bottom row at +x,+y'
SiPM_2.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_2.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_2.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_2.n_outside = n_hydraulic
SiPM_2.n_inside = n_hydraulic
SiPM_2.surface_type = 'normal'
SiPM_2.absorption = 1
surface_list.append(SiPM_2)

#SiPM 3
SiPM_3 = surface.surface()
SiPM_3.description = 'SiPM bottom row at x=0,+y'
SiPM_3.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_3.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_3.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_3.n_outside = n_hydraulic
SiPM_3.n_inside = n_hydraulic
SiPM_3.surface_type = 'normal'
SiPM_3.absorption = 1
surface_list.append(SiPM_3)

#SiPM 4
SiPM_4 = surface.surface()
SiPM_4.description = 'SiPM bottom row at -x,+y'
SiPM_4.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_4.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_4.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_4.n_outside = n_hydraulic
SiPM_4.n_inside = n_hydraulic
SiPM_4.surface_type = 'normal'
SiPM_4.absorption = 1
surface_list.append(SiPM_4)

# SiPM 5
SiPM_5 = surface.surface()
SiPM_5.description = 'SiPM bottom row at -x,y=0'
SiPM_5.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_5.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_5.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,0,:]<0)), (p.shape[0],-1))
SiPM_5.n_outside = n_hydraulic
SiPM_5.n_inside = n_hydraulic
SiPM_5.surface_type = 'normal'
SiPM_5.absorption = 1
surface_list.append(SiPM_5)

#SiPM 6
SiPM_6 = surface.surface()
SiPM_6.description = 'SiPM bottom row at -x,-y'
SiPM_6.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_6.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_6.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_6.n_outside = n_hydraulic
SiPM_6.n_inside = n_hydraulic
SiPM_6.surface_type = 'normal'
SiPM_6.absorption = 1
surface_list.append(SiPM_6)

#SiPM 7
SiPM_7 = surface.surface()
SiPM_7.description = 'SiPM bottom row at x=0,-y'
SiPM_7.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_7.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_7.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_7.n_outside = n_hydraulic
SiPM_7.n_inside = n_hydraulic
SiPM_7.surface_type = 'normal'
SiPM_7.absorption = 1
surface_list.append(SiPM_7)

#SiPM 8
SiPM_8 = surface.surface()
SiPM_8.description = 'SiPM bottom row at +x,-y'
SiPM_8.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_8.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_8.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r1[1]) & \
                                                 (p[:,2,:]>=SiPM_r1[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_8.n_outside = n_hydraulic
SiPM_8.n_inside = n_hydraulic
SiPM_8.surface_type = 'normal'
SiPM_8.absorption = 1
surface_list.append(SiPM_8)

## SiPMs in Second Row

# SiPM 9
SiPM_9 = surface.surface()
SiPM_9.description = 'SiPM second row at +x, y=0'
SiPM_9.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_9.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_9.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,0,:]>0)), (p.shape[0],-1))
SiPM_9.n_outside = n_hydraulic
SiPM_9.n_inside = n_hydraulic
SiPM_9.surface_type = 'normal'
SiPM_9.absorption = 1
surface_list.append(SiPM_9)

#SiPM 10
SiPM_10 = surface.surface()
SiPM_10.description = 'SiPM second row at +x,+y'
SiPM_10.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_10.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_10.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_10.n_outside = n_hydraulic
SiPM_10.n_inside = n_hydraulic
SiPM_10.surface_type = 'normal'
SiPM_10.absorption = 1
surface_list.append(SiPM_10)

#SiPM 11
SiPM_11 = surface.surface()
SiPM_11.description = 'SiPM second row at x=0,+y'
SiPM_11.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_11.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_11.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_11.n_outside = n_hydraulic
SiPM_11.n_inside = n_hydraulic
SiPM_11.surface_type = 'normal'
SiPM_11.absorption = 1
surface_list.append(SiPM_11)

#SiPM 12
SiPM_12 = surface.surface()
SiPM_12.description = 'SiPM second row at -x,+y'
SiPM_12.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_12.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_12.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_12.n_outside = n_hydraulic
SiPM_12.n_inside = n_hydraulic
SiPM_12.surface_type = 'normal'
SiPM_12.absorption = 1
surface_list.append(SiPM_12)

# SiPM 13
SiPM_13 = surface.surface()
SiPM_13.description = 'SiPM second row at -x,y=0'
SiPM_13.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_13.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_13.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,0,:]<0)), (p.shape[0],-1))
SiPM_13.n_outside = n_hydraulic
SiPM_13.n_inside = n_hydraulic
SiPM_13.surface_type = 'normal'
SiPM_13.absorption = 1
surface_list.append(SiPM_13)

#SiPM 14
SiPM_14 = surface.surface()
SiPM_14.description = 'SiPM second row at -x,-y'
SiPM_14.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_14.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_14.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_14.n_outside = n_hydraulic
SiPM_14.n_inside = n_hydraulic
SiPM_14.surface_type = 'normal'
SiPM_14.absorption = 1
surface_list.append(SiPM_14)

#SiPM 15
SiPM_15 = surface.surface()
SiPM_15.description = 'SiPM second row at x=0,-y'
SiPM_15.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_15.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_15.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_15.n_outside = n_hydraulic
SiPM_15.n_inside = n_hydraulic
SiPM_15.surface_type = 'normal'
SiPM_15.absorption = 1
surface_list.append(SiPM_15)

#SiPM 16
SiPM_16 = surface.surface()
SiPM_16.description = 'SiPM second row at +x,-y'
SiPM_16.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_16.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_16.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r2[1]) & \
                                                 (p[:,2,:]>=SiPM_r2[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_16.n_outside = n_hydraulic
SiPM_16.n_inside = n_hydraulic
SiPM_16.surface_type = 'normal'
SiPM_16.absorption = 1
surface_list.append(SiPM_16)

## SiPMs in Third row

# SiPM 17
SiPM_17 = surface.surface()
SiPM_17.description = 'SiPM third row at +x, y=0'
SiPM_17.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_17.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_17.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,0,:]>0)), (p.shape[0],-1))
SiPM_17.n_outside = n_hydraulic
SiPM_17.n_inside = n_hydraulic
SiPM_17.surface_type = 'normal'
SiPM_17.absorption = 1
surface_list.append(SiPM_17)

#SiPM 18
SiPM_18 = surface.surface()
SiPM_18.description = 'SiPM third row at +x,+y'
SiPM_18.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_18.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_18.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_18.n_outside = n_hydraulic
SiPM_18.n_inside = n_hydraulic
SiPM_18.surface_type = 'normal'
SiPM_18.absorption = 1
surface_list.append(SiPM_18)

#SiPM 19
SiPM_19 = surface.surface()
SiPM_19.description = 'SiPM third row at x=0,+y'
SiPM_19.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_19.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_19.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_19.n_outside = n_hydraulic
SiPM_19.n_inside = n_hydraulic
SiPM_19.surface_type = 'normal'
SiPM_19.absorption = 1
surface_list.append(SiPM_19)

#SiPM 20
SiPM_20 = surface.surface()
SiPM_20.description = 'SiPM third row at -x,+y'
SiPM_20.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_20.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_20.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_20.n_outside = n_hydraulic
SiPM_20.n_inside = n_hydraulic
SiPM_20.surface_type = 'normal'
SiPM_20.absorption = 1
surface_list.append(SiPM_20)

# SiPM 21
SiPM_21 = surface.surface()
SiPM_21.description = 'SiPM third row at -x,y=0'
SiPM_21.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_21.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_21.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,0,:]<0)), (p.shape[0],-1))
SiPM_21.n_outside = n_hydraulic
SiPM_21.n_inside = n_hydraulic
SiPM_21.surface_type = 'normal'
SiPM_21.absorption = 1
surface_list.append(SiPM_21)

#SiPM 22
SiPM_22 = surface.surface()
SiPM_22.description = 'SiPM third row at -x,-y'
SiPM_22.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_22.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_22.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_22.n_outside = n_hydraulic
SiPM_22.n_inside = n_hydraulic
SiPM_22.surface_type = 'normal'
SiPM_22.absorption = 1
surface_list.append(SiPM_22)

#SiPM 23
SiPM_23 = surface.surface()
SiPM_23.description = 'SiPM third row at x=0,-y'
SiPM_23.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_23.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_23.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_23.n_outside = n_hydraulic
SiPM_23.n_inside = n_hydraulic
SiPM_23.surface_type = 'normal'
SiPM_23.absorption = 1
surface_list.append(SiPM_23)

#SiPM 24
SiPM_24 = surface.surface()
SiPM_24.description = 'SiPM third row at +x,-y'
SiPM_24.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_24.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_24.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r3[1]) & \
                                                 (p[:,2,:]>=SiPM_r3[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_24.n_outside = n_hydraulic
SiPM_24.n_inside = n_hydraulic
SiPM_24.surface_type = 'normal'
SiPM_24.absorption = 1
surface_list.append(SiPM_24)

## SiPMs in Fourth row

# SiPM 25
SiPM_25 = surface.surface()
SiPM_25.description = 'SiPM fourth row at +x, y=0'
SiPM_25.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_25.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_25.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,0,:]>0)), (p.shape[0],-1))
SiPM_25.n_outside = n_hydraulic
SiPM_25.n_inside = n_hydraulic
SiPM_25.surface_type = 'normal'
SiPM_25.absorption = 1
surface_list.append(SiPM_25)

#SiPM 26
SiPM_26 = surface.surface()
SiPM_26.description = 'SiPM fourth row at +x,+y'
SiPM_26.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_26.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_26.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_26.n_outside = n_hydraulic
SiPM_26.n_inside = n_hydraulic
SiPM_26.surface_type = 'normal'
SiPM_26.absorption = 1
surface_list.append(SiPM_26)

#SiPM 27
SiPM_27 = surface.surface()
SiPM_27.description = 'SiPM fourth row at x=0,+y'
SiPM_27.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_27.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_27.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_27.n_outside = n_hydraulic
SiPM_27.n_inside = n_hydraulic
SiPM_27.surface_type = 'normal'
SiPM_27.absorption = 1
surface_list.append(SiPM_27)

#SiPM 28
SiPM_28 = surface.surface()
SiPM_28.description = 'SiPM fourth row at -x,+y'
SiPM_28.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_28.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_28.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,1,:]>0)), (p.shape[0],-1))
SiPM_28.n_outside = n_hydraulic
SiPM_28.n_inside = n_hydraulic
SiPM_28.surface_type = 'normal'
SiPM_28.absorption = 1
surface_list.append(SiPM_28)

# SiPM 29
SiPM_29 = surface.surface()
SiPM_29.description = 'SiPM fourth row at -x,y=0'
SiPM_29.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_29.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_29.inbounds_function =  lambda p: np.reshape(((p[:,1,:]<=SiPM_c1y[1]) & (p[:,1,:]>= SiPM_c1y[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,0,:]<0)), (p.shape[0],-1))
SiPM_29.n_outside = n_hydraulic
SiPM_29.n_inside = n_hydraulic
SiPM_29.surface_type = 'normal'
SiPM_29.absorption = 1
surface_list.append(SiPM_29)

#SiPM 30
SiPM_30 = surface.surface()
SiPM_30.description = 'SiPM fourth row at -x,-y'
SiPM_30.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_30.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_30.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c4[1]) & (p[:,0,:]>= SiPM_c4[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_30.n_outside = n_hydraulic
SiPM_30.n_inside = n_hydraulic
SiPM_30.surface_type = 'normal'
SiPM_30.absorption = 1
surface_list.append(SiPM_30)

#SiPM 23
SiPM_31 = surface.surface()
SiPM_31.description = 'SiPM fourth row at x=0,-y'
SiPM_31.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_31.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_31.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c3[1]) & (p[:,0,:]>= SiPM_c3[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_31.n_outside = n_hydraulic
SiPM_31.n_inside = n_hydraulic
SiPM_31.surface_type = 'normal'
SiPM_31.absorption = 1
surface_list.append(SiPM_31)

#SiPM 32
SiPM_32 = surface.surface()
SiPM_32.description = 'SiPM fourth row at +x,-y'
SiPM_32.shape = 'cylinder' #SiPMs are actually planar, but this makes the geometry much easier without much of a cost of accuracy.
SiPM_32.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
SiPM_32.inbounds_function =  lambda p: np.reshape(((p[:,0,:]<=SiPM_c2[1]) & (p[:,0,:]>= SiPM_c2[0]) & (p[:,2,:]<=SiPM_r4[1]) & \
                                                 (p[:,2,:]>=SiPM_r4[0]) & (p[:,1,:]<0)), (p.shape[0],-1))
SiPM_32.n_outside = n_hydraulic
SiPM_32.n_inside = n_hydraulic
SiPM_32.surface_type = 'normal'
SiPM_32.absorption = 1
surface_list.append(SiPM_32)

# The next 8 SiPMs are mounted on the cone right above the cylinder.
# some useful variables:
hsl_bottom = 2*rd_rad/(1+np.sqrt(2)) #horizontal side length of octagon SiPMs are mounted to (actually a cone, but this works)
sh = np.sqrt((rdcone_top-rd_top)**2+(rdcone_toprad-rd_rad)**2) #vertical side length of the idealized octagon of the cone
cntheta = np.arctan((rdcone_top-rd_top)/(rdcone_toprad-rd_rad)) # angle made with the sides of the octagon and vertical axis
SiPM_d_angled = np.sin(cntheta)*SiPM_d   # the dimensions of the SiPMs projected onto the z axis
ch = rdcone_top-rd_top # height of the cone
cnslope = (rdcone_top-rd_top)/(rdcone_toprad-rd_rad) #the slope of the line in the yz and xz plane made by the slope of the cone
                                                      #on top of the rd cyclinder
#center of SiPM for +x, +y tilted segment on cone. Can change some signs to make this work for all 4 of the awkward tilted parts
# format = [x-value,z-value]
SiPM_center = [ch/2*(1/cnslope)+rd_rad-np.sqrt(2)/4*hsl_bottom,ch/2]


#SiPM 33
SiPM_33 = surface.surface()
SiPM_33.description = 'SiPM upper cone at +x,y=0'
SiPM_33.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_33.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_33.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,1,:]>=SiPM_c1y[0])\
& (p[:,1,:]<=SiPM_c1y[1]) & (p[:,0,:]>0) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & (p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled)))\
, (p.shape[0], -1))
SiPM_33.n_outside = n_hydraulic
SiPM_33.n_inside = n_hydraulic
SiPM_33.surface_type = 'normal'
SiPM_33.absorption = 1
surface_list.append(SiPM_33)

#SiPM 34
SiPM_34 = surface.surface()
SiPM_34.description = 'SiPM upper cone at +x,+y'
SiPM_34.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_34.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_34.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,0,:]>(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom-np.sin(np.pi/4)*SiPM_d))) & (p[:,0,:]<(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom+np.sin(np.pi/4)*SiPM_d))) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & \
(p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled))) & (p[:,1,:]>0),(p.shape[0], -1))
SiPM_34.n_outside = n_hydraulic
SiPM_34.n_inside = n_hydraulic
SiPM_34.surface_type = 'normal'
SiPM_34.absorption = 1
surface_list.append(SiPM_34)

#SiPM 35
SiPM_35 = surface.surface()
SiPM_35.description = 'SiPM upper cone at x=0,+y'
SiPM_35.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_35.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_35.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,0,:]>=SiPM_c3[0])\
& (p[:,0,:]<=SiPM_c3[1]) & (p[:,1,:]>0) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & (p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled)))\
, (p.shape[0], -1))
SiPM_35.n_outside = n_hydraulic
SiPM_35.n_inside = n_hydraulic
SiPM_35.surface_type = 'normal'
SiPM_35.absorption = 1
surface_list.append(SiPM_35)

#SiPM 36
SiPM_36 = surface.surface()
SiPM_36.description = 'SiPM upper cone at -x,+y'
SiPM_36.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_36.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_36.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,0,:]<-(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom-np.sin(np.pi/4)*SiPM_d))) & (p[:,0,:]>-(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom+np.sin(np.pi/4)*SiPM_d))) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & \
(p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled))) & (p[:,1,:]>0),(p.shape[0], -1))
SiPM_36.n_outside = n_hydraulic
SiPM_36.n_inside = n_hydraulic
SiPM_36.surface_type = 'normal'
SiPM_36.absorption = 1
surface_list.append(SiPM_36)

#SiPM 37
SiPM_37 = surface.surface()
SiPM_37.description = 'SiPM upper cone at -x,y=0'
SiPM_37.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_37.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_37.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,1,:]>=SiPM_c5y[0])\
& (p[:,1,:]<=SiPM_c5y[1]) & (p[:,0,:]<0) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & (p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled)))\
, (p.shape[0], -1))
SiPM_37.n_outside = n_hydraulic
SiPM_37.n_inside = n_hydraulic
SiPM_37.surface_type = 'normal'
SiPM_37.absorption = 1
surface_list.append(SiPM_37)

#SiPM 38
SiPM_38 = surface.surface()
SiPM_38.description = 'SiPM upper cone at -x,-y'
SiPM_38.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_38.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_38.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,0,:]<-(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom-np.sin(np.pi/4)*SiPM_d))) & (p[:,0,:]>-(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom+np.sin(np.pi/4)*SiPM_d))) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & \
(p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled))) & (p[:,1,:]<0),(p.shape[0], -1))
SiPM_38.n_outside = n_hydraulic
SiPM_38.n_inside = n_hydraulic
SiPM_38.surface_type = 'normal'
SiPM_38.absorption = 1
surface_list.append(SiPM_38)

#SiPM 39
SiPM_39 = surface.surface()
SiPM_39.description = 'SiPM upper cone at x=0,-y'
SiPM_39.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_39.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_39.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,0,:]>=SiPM_c3[0])\
& (p[:,0,:]<=SiPM_c3[1]) & (p[:,1,:]<0) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & (p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled)))\
, (p.shape[0], -1))
SiPM_39.n_outside = n_hydraulic
SiPM_39.n_inside = n_hydraulic
SiPM_39.surface_type = 'normal'
SiPM_39.absorption = 1
surface_list.append(SiPM_39)

#SiPM 40
SiPM_40 = surface.surface()
SiPM_40.description = 'SiPM upper cone at +x,-y'
SiPM_40.shape = 'quadsurface'                       # QuadSurface needs to be tested
SiPM_40.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
SiPM_40.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top) & (p[:,0,:]>(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom-np.sin(np.pi/4)*SiPM_d))) & (p[:,0,:]<(p[:,2,:]*\
(1/cnslope)+(rd_rad-np.sqrt(2)/4*hsl_bottom+np.sin(np.pi/4)*SiPM_d))) & (p[:,2,:]>= ((sh/2)-(0.5*SiPM_d_angled))) & \
(p[:,2,:]<=((sh/2)+(0.5*SiPM_d_angled))) & (p[:,1,:]<0),(p.shape[0], -1))
SiPM_40.n_outside = n_hydraulic
SiPM_40.n_inside = n_hydraulic
SiPM_40.surface_type = 'normal'
SiPM_40.absorption = 1
surface_list.append(SiPM_40)
"""SBC Geometry in RayTracer style, adapted from SBC_MCNP_Eu148_vised.i and CreateArBCGeometry.m"""
import numpy as np
import math
import random
import RayTracer2
import surface
import matplotlib.pyplot as plt

# Material Parameters - indices of refraction
# n_argon = 1.45 # arXiv:1502.04213v4
# n_cf4 = 1.0004823 # gaseous!!! can't find liquid
# n_fused_quartz = 1.4424
# n_copper
# n_ss_304L = 2.75681  # typical value, not specific to 304L: https://www.filmetrics.com/refractive-index-database/Stainless-Steel#:~:text=For%20a%20typical%20sample%20of,nm%20are%202.75681%20and%203.792016.
# n_beryl
# n_air
# n_hdpe
"""n_target = 1.17;% n=1.224 for Ar @ 940nm, 90K
n_jar = 1.4512; % SiO2 @ 940nm
n_hydraulic = 1.22;% 1.21 @ 940nm, 146K ;  1.237 @ 7eV, 146K, another said 1.515;
n_pressurewindow = 1.7569; % Al2O3 @ 940nm
n_pressurewall = inf;
n_air = 1.00;"""
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

rd_rad = 12.5 # reflector - diffuser radius
rd_top = 0
rd_bot = -25
rdcone_top = 17.5 + 10
rdcone_toprad = 14.5 + 20
rdtopcone_apex = 21 + 10
rdtopcone_rad = 9.5
rdtopcone_bot = 19.5 + 10
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

rd_stcone_b = (rdcone_toprad - rdtopcone_rad) / (rdtopcone_bot - rdcone_top)
rd_stcone_z0 = rdtopcone_bot + (rdtopcone_rad / rd_stcone_b)
rd_stcone_Q = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, -rd_stcone_b**2]))
rd_stcone_P = np.array([0, 0, 2 * rd_stcone_b**2 * rd_stcone_z0])
rd_stcone_R = -(rd_stcone_b * rd_stcone_z0)**2

rd_topcone_b = rdtopcone_rad / (rdtopcone_apex - rdtopcone_bot)
rd_topcone_Q = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, -rd_topcone_b**2]))
rd_topcone_P = np.array([0, 0, 2 * rd_topcone_b**2 * rdtopcone_apex])
rd_topcone_R = -(rd_topcone_b * rdtopcone_apex)**2

rd_botcone_b = rdbotcone_rad / (rdbotcone_apex - rdbotcone_bot)
rd_botcone_Q = np.vstack(([1, 0, 0], [0, 1, 0], [0, 0, -rd_botcone_b**2]))
rd_botcone_P = np.array([0, 0, 2 * rd_botcone_b**2 * rdbotcone_apex])
rd_botcone_R = -(rd_botcone_b * rdbotcone_apex)**2



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
iiknuckle.inbounds_function = lambda p: np.reshape((p[:,2,:] > ijar_elevation) * (p[:,2,:] <= (z[3]+ijar_elevation)) *
                                                   ((p[:,0,:]**2 + p[:,1,:]** 2) > ((r1[3] - r2[3])**2)), (p.shape[0], -1))
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


## Viewports
sight_glass = surface.surface()
sight_glass.description = 'sight glass wall'
sight_glass.shape = 'cylinder'
sight_glass.param_list = [vp_center, vp_axis, vp_air_rad]
sight_glass.inbounds_function = lambda p: np.reshape(((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) > 0) *
                                                      (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                       (vp_nip_top + vp_flange_thick[1]))), (p.shape[0], -1))
sight_glass.n_outside = n_pressurewall
sight_glass.n_inside = n_air
sight_glass.surface_type = 'normal'
sight_glass.absorption = 1
surface_list.append(sight_glass)

camera_can_i = surface.surface()
camera_can_i.description = 'camera can inner wall'
camera_can_i.shape = 'cylinder'
camera_can_i.param_list = [vp_center, vp_axis, vp_can_rad - vp_can_wall]
camera_can_i.inbounds_function = lambda p: np.reshape(((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) >
                                                        (vp_nip_top + vp_flange_thick[1])) *
                                                       (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                        (vp_can_OAL + vp_flange_thick[1] + vp_nip_top))), (p.shape[0], -1))
camera_can_i.n_outside = n_pressurewall
camera_can_i.n_inside = n_air
camera_can_i.surface_type = 'normal'
camera_can_i.absorption = 1
surface_list.append(camera_can_i)

camera_can_o = surface.surface()
camera_can_o.description = 'camera can outer wall'
camera_can_o.shape = 'cylinder'
camera_can_o.param_list = [vp_center, vp_axis, vp_can_rad]
camera_can_o.inbounds_function = lambda p: np.reshape(((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) >
                                                        (vp_nip_top + vp_flange_thick[1] + vp_flange_thick[2])) *
                                                       (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                        (vp_can_OAL + vp_nip_top + vp_flange_thick[1] - vp_flange_thick[3]))), (p.shape[0], -1))
camera_can_o.n_outside = 1
camera_can_o.n_inside = n_pressurewall
camera_can_o.surface_type = 'normal'
camera_can_o.absorption = 1
surface_list.append(camera_can_o)

flange = surface.surface()
flange.description = 'flange outer edge'
flange.shape = 'cylinder'
flange.param_list = [vp_center, vp_axis, vp_flange_rad]
flange.inbounds_function = lambda p: np.reshape(((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) >
                                                        (-vp_flange_thick[0] + vp_nip_top)) *
                                                       (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                        (vp_nip_top + vp_flange_thick[1] + vp_flange_thick[2])))
                                                / ((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) >
                                                        (vp_nip_top + vp_flange_thick[1] + vp_can_OAL - vp_flange_thick[3])) *
                                                       (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                        (vp_nip_top + vp_flange_thick[1] + vp_can_OAL + vp_flange_thick[4]))), (p.shape[0], -1))
flange.n_outside = 1
flange.n_inside = n_pressurewall
flange.surface_type = 'normal'
flange.absorption = 1
surface_list.append(flange)

win_wall = surface.surface()
win_wall.description = 'window wall'
win_wall.shape = 'cylinder'
win_wall.param_list = [vp_center, vp_axis, vp_win_rad]
win_wall.inbounds_function = lambda p:np.reshape((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) >
                                                        (-vp_win_thick)) *
                                                       (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                        0), (p.shape[0], -1))
win_wall.n_outside = n_hydraulic
win_wall.n_inside = n_pressurewindow
win_wall.surface_type = 'normal'
win_wall.absorption = 1
surface_list.append(win_wall)

win_retainer = surface.surface()
win_retainer.description = 'window retainer outer wall'
win_retainer.shape = 'cylinder'
win_retainer.param_list = [vp_center, vp_axis, vp_win_rad]
win_retainer.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) >
                                                        0) *
                                                       (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                        vp_nip_top), (p.shape[0], -1))
win_retainer.n_outside = n_hydraulic
win_retainer.n_inside = n_pressurewall
win_retainer.surface_type = 'normal'
win_retainer.absorption = 1
surface_list.append(win_retainer)

pv_n_wall = surface.surface()
pv_n_wall.description = 'pressure vessel nipple wall'
pv_n_wall.shape = 'cylinder'
pv_n_wall.param_list = [vp_center, vp_axis, vp_win_rad]
pv_n_wall.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) >
                                                        (-vp_flange_thick[0] + vp_nip_top)) *
                                                       (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) * np.tile(vp_axis, p.shape[0], 1, p.shape[2]), 1) <=
                                                        vp_nip_top), (p.shape[0], -1))
pv_n_wall.n_outside = n_pressurewall
pv_n_wall.n_inside = n_hydraulic
pv_n_wall.surface_type = 'normal'
pv_n_wall.absorption = 1
surface_list.append(pv_n_wall)

vp_air = surface.surface()
vp_air.description = 'air side of viewport'
vp_air.shape = 'plane'
vp_air.param_list = [vp_center, vp_axis]
vp_air.inbounds_function = lambda p: np.reshape(np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2])))**2, 1) <=
                                                vp_air_rad**2, (p.shape[0], -1))
vp_air.n_outside = n_air
vp_air.n_inside = n_pressurewindow
vp_air.surface_type = 'normal'
vp_air.absorption = 0
surface_list.append(vp_air)

vp_hydraulic = surface.surface()
vp_hydraulic.description = 'hydraulic side of viewport'
vp_hydraulic.shape = 'plane'
vp_hydraulic.param_list = [vp_center - (vp_axis * vp_win_thick), vp_axis]
vp_hydraulic.inbounds_function = lambda p: np.reshape(np.sum((p - np.tile(vp_center - (vp_axis * vp_win_thick), (p.shape[0], 1, p.shape[2])))**2, 1) <=
                                                vp_win_rad**2, (p.shape[0], -1))
vp_hydraulic.n_outside = n_pressurewindow
vp_hydraulic.n_inside = n_hydraulic
vp_hydraulic.surface_type = 'normal'
vp_hydraulic.absorption = 0
surface_list.append(vp_hydraulic)

vp_retainer = surface.surface()
vp_retainer.description = 'viewport retainer'
vp_retainer.shape = 'plane'
vp_retainer.param_list = [vp_center, vp_axis]
vp_retainer.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2])))**2, 1) >
                                                vp_air_rad**2) *
                                                (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2])))**2, 1) <=
                                                vp_win_rad**2, (p.shape[0], -1)))
vp_retainer.n_outside = n_pressurewall
vp_retainer.n_inside = n_pressurewindow
vp_retainer.surface_type = 'normal'
vp_retainer.absorption = 1
surface_list.append(vp_retainer)

nip_bottom = surface.surface()
nip_bottom.description = 'nipple bottom'
nip_bottom.shape = 'plane'
nip_bottom.param_list = [vp_center - vp_axis * (vp_flange_thick[0] -vp_nip_top), vp_axis]
nip_bottom.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center - vp_axis*(vp_flange_thick[0] - vp_nip_top), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                > vp_nip_rad**2) *
                                                (np.sum((p - np.tile(vp_center - vp_axis*(vp_flange_thick[0] - vp_nip_top), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                <= vp_flange_rad**2), (p.shape[0], -1))
nip_bottom.n_outside = n_pressurewall
nip_bottom.n_inside = n_hydraulic
nip_bottom.surface_type = 'normal'
nip_bottom.absorption = 1
surface_list.append(nip_bottom)

nip_top = surface.surface()
nip_top.description = 'nipple top'
nip_top.shape = 'plane'
nip_top.param_list = [vp_center + vp_axis * vp_nip_top, vp_axis]
nip_top.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center + vp_axis*vp_nip_top, (p.shape[0], 1, p.shape[2])))**2, 1)
                                                > vp_win_rad**2) *
                                                (np.sum((p - np.tile(vp_center + vp_axis*vp_nip_top, (p.shape[0], 1, p.shape[2])))**2, 1)
                                                <= vp_flange_rad**2), (p.shape[0], -1))
nip_top.n_outside = n_pressurewall
nip_top.n_inside = n_hydraulic
nip_top.surface_type = 'normal'
nip_top.absorption = 1
surface_list.append(nip_top)

can_bot = surface.surface()
can_bot.description = 'can bot'
can_bot.shape = 'plane'
can_bot.param_list = [vp_center + vp_axis * (vp_nip_top + vp_flange_thick[1]), vp_axis]
can_bot.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1]), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                > vp_air_rad**2) *
                                                (np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1]), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                <= vp_flange_rad**2), (p.shape[0], -1))
can_bot.n_outside = n_air
can_bot.n_inside = n_pressurewall
can_bot.surface_type = 'normal'
can_bot.absorption = 1
surface_list.append(can_bot)

can_bot_top = surface.surface()
can_bot_top.description = 'can bottom top'
can_bot_top.shape = 'plane'
can_bot_top.param_list = [vp_center + vp_axis * (vp_nip_top + vp_flange_thick[1] + vp_flange_thick[2]), vp_axis]
can_bot_top.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1] + vp_flange_thick[2]), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                > vp_can_rad**2) *
                                                (np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1] + vp_flange_thick[2]), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                <= vp_flange_rad**2), (p.shape[0], -1))
can_bot_top.n_outside = 1
can_bot_top.n_inside = n_pressurewall
can_bot_top.surface_type = 'normal'
can_bot_top.absorption = 1
surface_list.append(can_bot_top)


can_top_bot = surface.surface()
can_top_bot.description = 'can top bottom'
can_top_bot.shape = 'plane'
can_top_bot.param_list = [vp_center + vp_axis*(vp_nip_top + vp_flange_thick[1] + vp_can_OAL - vp_flange_thick[3]), vp_axis]
can_top_bot.inbounds_function = lambda p: np.reshape((np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1] + vp_can_OAL - vp_flange_thick[3]), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                > vp_can_rad**2) *
                                                (np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1] + vp_can_OAL - vp_flange_thick[3]), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                <= vp_flange_rad**2), (p.shape[0], -1))
can_top_bot.n_outside = n_pressurewall
can_top_bot.n_inside = 1
can_top_bot.surface_type = 'normal'
can_top_bot.absorption = 1
surface_list.append(can_top_bot)

can_top = surface.surface()
can_top.description = 'can top'
can_top.shape = 'plane'
can_top.param_list = [vp_center + vp_axis*(vp_nip_top + vp_flange_thick[1] + vp_can_OAL), vp_axis]
can_top.inbounds_function = lambda p: np.reshape(np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1] + vp_can_OAL), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                <= vp_flange_rad**2, (p.shape[0], -1))
can_top.n_outside = n_pressurewall
can_top.n_inside = n_air
can_top.surface_type = 'normal'
can_top.absorption = 1
surface_list.append(can_top)

can_vtop = surface.surface()
can_vtop.description = 'can very top'
can_vtop.shape = 'plane'
can_vtop.param_list = [vp_center + vp_axis*(vp_nip_top + vp_flange_thick[1] + vp_can_OAL + vp_flange_thick[4]), vp_axis]
can_vtop.inbounds_function = lambda p: np.reshape(np.sum((p - np.tile(vp_center + vp_axis*(vp_nip_top+vp_flange_thick[1] + vp_can_OAL + vp_flange_thick[4]), (p.shape[0], 1, p.shape[2])))**2, 1)
                                                <= vp_flange_rad**2, (p.shape[0], -1))
can_vtop.n_outside = 1
can_vtop.n_inside = n_pressurewall
can_vtop.surface_type = 'normal'
can_vtop.absorption = 1
surface_list.append(can_vtop)

######## DELETE TEMPLATE
can_bot = surface.surface()
can_bot.description = 'can bot'
can_bot.shape = 'plane'
can_bot.param_list = []
can_bot.inbounds_function = lambda p: np.reshape()
can_bot.n_outside = n_air
can_bot.n_inside = n_pressurewall
can_bot.surface_type = 'normal'
can_bot.absorption = 1
surface_list.append(can_bot)


## other black surfaces to trap rays
rd = surface.surface()
rd.description = 'reflector/diffuser'
rd.shape = 'cylinder'
rd.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), rd_rad]
rd.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_bot) * (p[:, 2, :] <= rd_top), (p.shape[0], -1))
rd.n_outside = n_hydraulic
rd.n_inside = n_hydraulic
rd.surface_type = 'normal'
rd.absorption = 1
surface_list.append(rd)

rd_cone = surface.surface()
rd_cone.description = 'reflector/diffuser cone'
rd_cone.shape = 'quadsurface'                       # QuadSurface needs to be tested
rd_cone.param_list = [rd_cone_Q, rd_cone_P, rd_cone_R]
rd_cone.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rd_top) * (p[:, 2, :] < rdcone_top), (p.shape[0], -1))
rd_cone.n_outside = n_hydraulic
rd_cone.n_inside = n_hydraulic
rd_cone.surface_type = 'normal'
rd_cone.absorption = 1
surface_list.append(rd_cone)

rd_scone = surface.surface()
rd_scone.description = 'reflector/diffuser strip cone'
rd_scone.shape = 'quadsurface'                       # QuadSurface needs to be tested
rd_scone.param_list = [rd_stcone_Q, rd_stcone_P, rd_stcone_R]
rd_scone.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdcone_top) * (p[:, 2, :] < rdtopcone_bot) *
                                                  ((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2])))**2, 1) -
                                                    (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) *
                                                    np.tile(vp_axis, (p.shape[0], 1, p.shape[2])), 1)**2)) > vp_nip_rad**2), (p.shape[0], -1))
rd_scone.n_outside = n_hydraulic
rd_scone.n_inside = n_hydraulic
rd_scone.surface_type = 'normal'
rd_scone.absorption = 1
surface_list.append(rd_scone)

rd_topcone = surface.surface()
rd_topcone.description = 'reflector/diffuser topcone'
rd_topcone.shape = 'quadsurface'                       # QuadSurface needs to be tested
rd_topcone.param_list = [rd_topcone_Q, rd_topcone_P, rd_topcone_R]
rd_topcone.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdtopcone_bot) * (p[:, 2, :] < rdtopcone_apex), (p.shape[0], -1))
rd_topcone.n_outside = n_hydraulic
rd_topcone.n_inside = n_hydraulic
rd_topcone.surface_type = 'normal'
rd_topcone.absorption = 1
surface_list.append(rd_topcone)

rd_botcone = surface.surface()
rd_botcone.description = 'reflector/diffuser botcone'
rd_botcone.shape = 'quadsurface'                       # QuadSurface needs to be tested
rd_botcone.param_list = [rd_botcone_Q, rd_botcone_P, rd_botcone_R]
rd_botcone.inbounds_function = lambda p: np.reshape((p[:, 2, :] > rdbotcone_bot) * (p[:, 2, :] < rdbotcone_apex), (p.shape[0], -1))
rd_botcone.n_outside = n_hydraulic
rd_botcone.n_inside = n_hydraulic
rd_botcone.surface_type = 'normal'
rd_botcone.absorption = 1
surface_list.append(rd_botcone)

pv_cyl_out = surface.surface()
pv_cyl_out.description = 'PV - cylinder outer wall'
pv_cyl_out.shape = 'cylinder'
pv_cyl_out.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), pv_rad]
pv_cyl_out.inbounds_function = lambda p: np.reshape((p[:, 2, :] > pv_bot) * (p[:, 2, :] < pv_top), (p.shape[0], -1))
pv_cyl_out.n_outside = 1
pv_cyl_out.n_inside = n_pressurewall
pv_cyl_out.surface_type = 'normal'
pv_cyl_out.absorption = 1
surface_list.append(pv_cyl_out)

pv_cyl_in = surface.surface()
pv_cyl_in.description = 'PV - cylinder inner wall'
pv_cyl_in.shape = 'cylinder'
pv_cyl_in.param_list = [np.array([0, 0, 0]), np.array([0, 0, 1]), pv_rad - pv_thick]
pv_cyl_in.inbounds_function = lambda p: np.reshape((p[:, 2, :] > pv_bot) * (p[:, 2, :] < pv_top), (p.shape[0], -1))
pv_cyl_in.n_outside = n_pressurewall
pv_cyl_in.n_inside = n_hydraulic
pv_cyl_in.surface_type = 'normal'
pv_cyl_in.absorption = 1
surface_list.append(pv_cyl_in)

pv_out_top = surface.surface()
pv_out_top.description = 'PV - outer top'
pv_out_top.shape = 'quadsurface'
pv_out_top.param_list = [head_out_Q, head_out_P, head_out_R]
pv_out_top.inbounds_function = lambda p: np.reshape((p[:, 2, :] > pv_top) *
                                                    ((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2])))**2, 1) -
                                                      (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) *
                                                    np.tile(vp_axis, (p.shape[0], 1, p.shape[2])), 1)**2)) >
                                                     vp_flange_rad**2), (p.shape[0], -1))
pv_out_top.n_outside = 1
pv_out_top.n_inside = n_pressurewall
pv_out_top.surface_type = 'normal'
pv_out_top.absorption = 1
surface_list.append(pv_out_top)

pv_in_top = surface.surface()
pv_in_top.description = 'PV - inner top'
pv_in_top.shape = 'quadsurface'
pv_in_top.param_list = [head_in_Q, head_in_P, head_in_R]
pv_in_top.inbounds_function = lambda p: np.reshape((p[:, 2, :] > pv_top) *
                                                    ((np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2])))**2, 1) -
                                                      (np.sum((p - np.tile(vp_center, (p.shape[0], 1, p.shape[2]))) *
                                                    np.tile(vp_axis, (p.shape[0], 1, p.shape[2])), 1)**2)) >
                                                     vp_flange_rad**2), (p.shape[0], -1))
pv_in_top.n_outside = n_pressurewall
pv_in_top.n_inside = n_hydraulic
pv_in_top.surface_type = 'normal'
pv_in_top.absorption = 1
surface_list.append(pv_in_top)

pv_bottom = surface.surface()
pv_bottom.description = 'PV - bottom'
pv_bottom.shape = 'plane'
pv_bottom.param_list = [np.array([0, 0, pv_bot]), np.array([0, 0, -1])]
pv_bottom.inbounds_function = lambda p: np.reshape((p[:, 0, :]**2 + p[:, 1, :]**2) <= pv_rad**2, (p.shape[0], -1))
pv_bottom.n_outside = n_pressurewall
pv_bottom.n_inside = n_hydraulic
pv_bottom.surface_type = 'normal'
pv_bottom.absorption = 1
surface_list.append(pv_bottom)


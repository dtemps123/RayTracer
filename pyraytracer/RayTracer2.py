'''#function ray_interfaces = RayTracer2(ray_startingpoints, rays, surfacelist, ...
# %     max_scatters, min_travel_length, follow_threshold, tir_handling, full_output, singlechild, output_raytable)
# %
# % RayTracer2 propagates a set of rays through a geometry of reflecting and
# % refracting surfaces.  The function inputs a set of initial rays and a
# % list of surfaces.  It loops over surfaces, finding the intersection point
# % of each ray with each surface, to determine the next scattering point of
# % each ray.  Reflected and refracted rays are generated according to the
# % surface properties.  The function iterates with the new set of rays,
# % until all rays drop below threshold or the maximum number of scatters is
# % reached.  The history of each ray is output in a structure array.
# %
# % Advances since RayTracer include bulk absorption and Rayleigh scattering,
# % advanced surface modeling (following Geant4's Optical Photon "unified"
# % surface physics) and maybe some other things I haven't thought of yet.
# % This *is* backwards compatible -- if you give it the same input as
# % RayTracer it will work just fine -- but it may be slightly slower.
# %
# % inputs:
# %           ray_startingpoints  -  N-by-3 matrix, where N is the number of
# %                                    initial rays to follow, giving the
# %                                    starting point for each ray.
# %           rays                -  N-by-10 matrix giving the initial
# %                                    direction, intensity, and polarization
# %                                    of each ray.  The first 3 columns give
# %                                    the forward direction of the ray
# %                                    (these will be normalized if they
# %                                    aren't already), columns 4-6 give a
# %                                    direction non-parallel to the ray that
# %                                    defines the s1 polarization axis
# %                                    (these will be made normal to the ray
# %                                    direction and normalized, if they
# %                                    aren't already), and columns 7-10 are
# %                                    the stokes parameters s0-s3, giving
# %                                    the intensity and polarization of the
# %                                    ray ( s0 gives the total intensity,
# %                                    and s0^2 >= s1^2 + s2^2 + s3^2, see
# %                                    7.2 in Jackson for more details)
# %           surfacelist         -  A structure array defining the geometry
# %                                    of scattering surfaces -- see
# %                                    Create2LGeometry for details
# %           max_scatters        -  The maximum number of scatters to
# %                                    propagate rays through (the simulation
# %                                    may stop before this point if there
# %                                    are no rays above threshold left to
# %                                    follow).
# %           min_travel_length   -  A minimum travel length between
# %                                    scatters.  This prevents rounding
# %                                    errors from causing a ray to scatter
# %                                    multiple times at the same point.
# %                                    Rays with legitimate travel lengths
# %                                    below this value will be INCORRECTLY
# %                                    RECONSTRUCTED, so keep this small
# %                                    (~1e-5 times your smallest dimension
# %                                    is probably sufficient for most
# %                                    geometries.)
# %           follow_threshold    -  Refracted or reflected rays with an s0
# %                                    below follow_threshold(1) or
# %                                    follow_threshold(2), respectively,
# %                                    will not be followed.  If
# %                                    follow_threshold is a scalar, the same
# %                                    threshold is used for both.
# %           tir_handling        -   This determines what the refracted_rays
# %                                     output is in the case of total internal
# %                                     reflection.  The default (-1) gives a
# %                                     refracted ray tangent to the surface with
# %                                     zero intensity.  Any value >=0 will give
# %                                     a ray with the same direction and
# %                                     polarization as the reflected ray, with
# %                                     intensity equal to the reflected
# %                                     intensity times tir_handling.  This lets
# %                                     you treat tir-rays like refracted rays,
# %                                     which can be handy in geometry sims.
# %                                     NOTE -- if follow_threshold(2) is
# %                                     bigger than max(rays(:,7)) then
# %                                     default tir_handling=1.
# %           full_output         -   If false, the ray_interfaces output is
# %                                     not populated.  Default is true.
# %           singlechild         -   If true, then reflected/refracted rays
# %                                     are never both followed, rather one
# %                                     is chosen by a dice roll, and the ray
# %                                     index is always positive.  If false,
# %                                     this follows the old RayTracer
# %                                     standard of following both reflected
# %                                     and refracted rays, with purely
# %                                     refracted trajectories only getting
# %                                     the positive index.  Default is true.
# %           output_raytable     -   If false, the raytable output is not
# %                                     populated.  Default is false.
# %
# % output:
# %       ray_interfaces  -  a structure array, where the array index is the
# %                            scatter number.  Each element in the array has
# %                            the following fields:
# %           incoming_ray        -  M-by-10 array, where M is the number of
# %                                    rays scattering in this iteration,
# %                                    giving the direction, intensity, and
# %                                    polarization of the incoming rays.
# %                                    Rays that do not scatter are not
# %                                    reported (to report all rays, enclose
# %                                    your geometry in an absorbing box, for
# %                                    example).
# %           reflected_ray       -  M-by-10 array, giving the direction,
# %                                    intensity, and polarization of the
# %                                    reflected rays
# %           refracted_ray       -  M-by-10 array, giving the direction,
# %                                    intensity, and polarization of the
# %                                    refracted rays
# %           intersection_point  -  M-by-3 array, giving the points where
# %                                    the incoming rays scatter
# %           surface_normal      -  M-by-3 array, giving the
# %                                    backward-pointing surface normal at
# %                                    the intersection_point
# %           ray_index           -  M-by-1 vector, giving the index of the
# %                                    incoming ray (the input rays are
# %                                    numbered 1:N) -- a negative ray_index
# %                                    means the ray has undergone at least
# %                                    one reflection in its history, so
# %                                    there will be at most one ray with a
# %                                    given positive index
# %           surface_index       -  M-by-1 vector, giving the index of the
# %                                    scattering surface hit by each
# %                                    incoming_ray, where the index
# %                                    indicates an element of the input
# %                                    surfacelist.  Negative values indicate
# %                                    an outward ray, positive values an
# %                                    inward ray, as defined in the
# %                                    surface geometry.
                                     NOTE: Unlike usual indexing, index 1 refers to the first surface in your list (index 0 does not, issue
                                     with having a -0 required some changes.)
# %           distace_traveled    -  M-by-1 vector, giving the distance
# %                                    traveled by the ray since its last
# %                                    scatter
# %           n_incident          -  M-by-1 vector, giving the index of
# %                                    refraction for the incoming ray
# %           n_transmitted       -  M-by-1 vector, giving the index of
# %                                    refraction for the refracted ray
# %           bulkabs_incident    -  M-by-1 vector, giving the absorption
# %                                    length for the incoming ray
# %           bulkabs_transmitted -  M-by-1 vector, giving the absorbtion
# %                                    length for the refracted ray
# %           rayleigh_incident   -  M-by-1 vector, giving the rayleigh
# %                                    scattering length for the incoming ray
# %           rayleigh_transmitted-  M-by-1 vector, giving the rayleigh
# %                                    scattering length for the refracted ray
# %
# %       absorption_table        -  Array of size [K, 5, S, 2], where
# %                                    K=max_scatters, S=length(surfacelist)
# %                                    Giving total intensity absorbed at
# %                                    each scatter step.  2nd index
# %                                    separates:
# %                                       1 - Surface absorption
# %                                       2 - Bulk absorption
# %                                       3 - Escaped geometry
# %                                       4 - Dropped below threshold
# %                                       5 - Still following
# %                                     3rd index identifies the surface
# %                                     the ray is interacting with in that
# %                                     step, 4th index indicates the ray's
# %                                     orientation with respect to that
# %                                     surface.  For `Escaped geometry'
# %                                     rays, the previously hit surface is
# %                                     indicated.
# %       raytable                -  Array of size [K+1, N, 13] giving the
# %                                     details of each ray's path through
# %                                     the geometry, following positive ray
# %                                     indices only.  First index is scatter
# %                                     index, starting with initial
# %                                     condition.  Second index is ray
# %                                     index.  For third index, columns 1:3
# %                                     give XYZ position, 4:13 give ray
# %                                     details (direction, intensity,
# %                                     polarization, as in rays input).
# %                                 
# % 
# % RayTracer:  12/17/09, CED
# % RayTracer2:  8/17/16, CED'''
import numpy as np
import numpy.matlib
import rayInterfaces
import RayleighScatteringClass
import RefractionReflectionAtInterface
import UnifiedReflectorModel
import IntersectFunction
import surface

# noinspection PyRedundantParentheses
def RayTracer2(ray_startingpoints, rays, surfacelist = [], max_scat = 10, min_travel_len = np.spacing(np.float64(1)), follow_thresh = np.array([0, 0], dtype=np.float64), tir_handling = [], full_output = True, singlechild = True, output_raytable = False):
    #ignore divide by zero / nan warnings
    np.seterr(divide='ignore', invalid='ignore')

    #instantiate return value containers
    ray_interfaces = []
    raytable = []

    #parse inputs into variable
    if (np.size(follow_thresh) == 1 and follow_thresh == 0):
        follow_threshold = follow_thresh
    elif (np.size(follow_thresh) == 1 and follow_thresh != 0):
        #MAKE NEW VARIABLE for follow threshold
        follow_threshold = np.array([follow_thresh,follow_thresh])
    else:
        follow_threshold = follow_thresh[0:2]
        

    if np.size(min_travel_len) == 1 and min_travel_len == np.spacing(np.float64(1)):
        min_travel_length = min_travel_len
    else:
        min_travel_length = min_travel_len[0]


    if isinstance(max_scat, int) and max_scat == 10:
        max_scatters = max_scat
    elif isinstance(max_scat, int):
        max_scatters = max_scat
    else:
        max_scatters = max_scat[0]


#    %% error check inputs
    if np.size(surfacelist) == 0 or np.array(ray_startingpoints).ndim !=2 or np.array(rays).ndim !=2 or np.size(ray_startingpoints,1) !=3 or np.size(rays,1)!=10 or np.size(ray_startingpoints,0)!=np.size(rays,0):
        raise Exception('Improper input to RayTracer2')

    numrays = np.size(rays,0)
    # Normalize rays
    rays[:,0:3] = rays[:,0:3] / np.abs(np.sqrt(np.sum(rays[:,0:3]**2, axis=1, keepdims=True)))
    rays[:,3:6] = rays[:,3:6] / np.abs(np.sqrt(np.sum(rays[:,3:6]**2, axis=1, keepdims=True)))

#    %% initialize absorption table
    absorption_table = np.zeros((max_scatters, 5, len(surfacelist)+1, 2))
    if output_raytable:
        raytable = np.zeros((max_scatters+1, np.size(ray_startingpoints, 0), 13))
        raytable[0, :, 0:3] = ray_startingpoints
        raytable[0, :, 3:13] = rays

#   Adds a null row to the surface list to avoid index -0 issues.
    null_surface = surface.surface()
    surfacelist = np.insert(surfacelist,0,null_surface,axis=0)        
        
#    %% check optional fields in surface list, necessary for backwards compatibility
    bulk_props = np.array(['abslength_outside', 'abslength_inside', 'rayleigh_outside', 'rayleigh_inside'], dtype=object)
    
    for i_s in range(1,np.size(surfacelist)):
        for i_p in range(np.size(bulk_props)):
        
            #this doesnt need to be a try statement because python uses short circuiting when evaluating boolean statements
            if ((~hasattr(surfacelist, bulk_props[i_p])) or (np.size(getattr(surfacelist[i_s],bulk_props[i_p])) == 0)):
                setattr(surfacelist[i_s],bulk_props[i_p],np.inf)

        
        # if (~hasattr(surfacelist, 'unifiedparams') or (np.size(surfacelist[i_s].unifiedparams) == 0)):
        #     surfacelist[i_s].unifiedparams = [0, 1, 0, 1, 0]
        try: # above was replacing ALL unifiedparams with [0, 1, 0, 1, 0]
            getattr(surfacelist[i_s], 'unifiedparams')
            #print(getattr(surfacelist[i_s], 'unifiedparams'))
        except AttributeError:
            surfacelist[i_s].unifiedparams = [0, 1, 0, 1, 0]

#    %% now really set default tir_handling
    if (np.size(tir_handling) == 0):
        if follow_threshold[0] > np.amax(rays[:,6]):
            tir_handling = 1
        else:
            tir_handling = -1

#    %% initialize raylist
    p_start = np.copy(ray_startingpoints) # Does this need to be copied?
    incoming_rays = np.copy(rays)
    ray_index = np.array(range(numrays))
    smix_last = np.ones(numrays, dtype=int)
    six_last = np.zeros((np.size(p_start,0), 1))

#    %%  follow rays
    num_scatters = 0
    while (len(ray_index) != 0):

        if num_scatters >= max_scatters:
            break
        num_scatters += 1

        
#        %% find next scattering surface for each ray
#        % initialize the variables containing scatter properties
        p_next = np.zeros(p_start.shape)         # scatter point
        l_next = np.zeros(p_start.shape[0]) + np.inf  # distance to scatter
        s_next = np.zeros(p_start.shape)            # optical surface normal at scatter
        sm_next = np.zeros(p_start.shape)       # mechanical surface normal at scatter
        n_next = np.zeros((p_start.shape[0],2))        # indices of refraction at scatter
        abs_next = np.zeros(p_start.shape[0])      # absorption at scatter
        six_next = np.zeros(p_start.shape[0])      # scatter surface index
        surfacetype_next = np.zeros(p_start.shape[0])  # surface type
        unifiedsurface_next = np.zeros((p_start.shape[0],5))  # surface parameters at scatter
        rayleigh_next = np.zeros((p_start.shape[0],2)) + np.inf  # Rayleigh scattering length before/after scatter
        abslength_next = np.zeros((p_start.shape[0],2)) + np.inf  # Bulk Absorption length before/after scatter

#        % loop over surfaces
        for n in range(1,len(surfacelist)):
#            % find scatter point, normal, distance, and orientation for
#            % scatters on this surface
            [p_intersect, s_normal, l_ray, s_orientation] = IntersectFunction.IntersectFunction(surfacelist[n],p_start,incoming_rays[:,0:3])

            s_orientation = np.array(s_orientation) #convert into np array
            
#            % save the actual normal vector, because that will be important
            s_normal_mech = np.copy(s_normal) #%N x 3 x M matrix
            
#            % adjust the surface normal as appropriate for the surface_type
            int_surfacetype = 0
            if (surfacelist[n].surface_type == 'diffuse'): #% reflect transmitted portion of incoming ray back into 2pi
                int_surfacetype = 1
                
            if (surfacelist[n].surface_type == 'unified'): #% mimics Geant4 UNIFIED surface for OpticalPhotons
                int_surfacetype = 2
                
            if (surfacelist[n].surface_type == 'retro'): #% reflect incoming ray back where it came from
                s_normal = np.tile(-1 * np.array(incoming_rays[:,0:3]), (1, 1, np.size(p_intersect,2)))
            
#            % find the closest intersection point in the inbounds portion of
#            % this surface (counting only positive, real distances and ignoring
#            % glancing blows)

#            valid_intersections = (surfacelist[n].inbounds_function(p_intersect)) and ( np.imag(l_ray)==0 ) and (s_orientation != 0) and (not np.isnan(l_ray)) and (l_ray < np.inf) and (l_ray > np.matlib.repmat(min_travel_length * int(six_last==n),1,np.size(l_ray,1)))
            # valid_intersection = np.logical_and.reduce(((surfacelist[n].inbounds_function(p_intersect)),(np.equal(np.imag(l_ray), 0)),(np.not_equal(s_orientation, 0)),(~np.isnan(l_ray)),(np.less(l_ray, np.full(l_ray.shape, np.inf))),(np.greater(l_ray, np.matlib.repmat((min_travel_length * np.equal(six_last, n)),1,l_ray.shape[1])))))
            valid_intersection = np.logical_and.reduce(((surfacelist[n].inbounds_function(p_intersect)),(np.equal(np.imag(l_ray), 0)),(np.not_equal(s_orientation, 0)),(~np.isnan(l_ray)),(np.less(l_ray, np.full(l_ray.shape, np.inf))),(np.greater(l_ray, (min_travel_length * np.equal(six_last, n))))))

            l_ray[~valid_intersection] = np.inf
            ix, l_ray = l_ray.argmin(1), l_ray.min(1)
            #l_ray is an array of minimums
            #ix is the indexes of the minimums


#            % find the intersection points, surface normals, and orientations
#            % assosciated with these scatters -- this is all just array
#            % manipulation to extract the rows assosciated with the scatter
#            % points found above
            ixlist = np.tile(np.reshape(range(0,p_intersect.shape[2]),(1,-1,1)),(3, 1, p_intersect.shape[0])) #changed from 1:size(p_intersect,2) to python indexed

                        
            ixcut = (ixlist == np.tile(np.reshape(ix,(1,1,-1)),(3, np.size(p_intersect,2), 1)))
            p_intersect = np.transpose(p_intersect,(1, 2, 0))
            p_intersect = np.transpose(np.reshape(p_intersect[ixcut],(3,-1)))
            s_normal = np.transpose(s_normal,(1, 2, 0))
            
            s_normal = np.transpose(np.reshape(s_normal[ixcut],(3,-1)))
            s_normal_mech = np.transpose(s_normal_mech,(1, 2, 0))
            s_normal_mech = np.transpose(np.reshape(s_normal_mech[ixcut],(3,-1)))

            # ixlist = np.matlib.repmat(np.reshape(range(0,s_orientation.shape[1]),(-1,1)), 1, s_orientation.shape[0])
            ixlist = np.tile(np.reshape(range(0, s_orientation.shape[1]), (-1, 1)), (1, s_orientation.shape[0])) #Check tile
            # ixcut = (ixlist == np.tile(np.transpose(ix[:]),(np.size(s_orientation,1), 1))) # is the : in ix[] necessary? both (n,). same with tranpose
            ixcut = np.equal(ixlist, ix[np.newaxis, :])

            s_orientation = np.transpose(s_orientation)
            s_orientation = np.reshape(s_orientation[ixcut],(-1,1))

            # n_before_after = repmat([surfacelist(n).n_outside surfacelist(n).n_inside],length(s_orientation),1); MATLAB
            # n_before_after(s_orientation<0,:) = repmat([surfacelist(n).n_inside surfacelist(n).n_outside],sum(s_orientation<0),1);
            #n_before_after = np.matlib.repmat(np.array([surfacelist[n].n_outside, surfacelist[n].n_inside]),len(s_orientation),1)
            n_before_after = np.tile(np.array([surfacelist[n].n_outside, surfacelist[n].n_inside]), (len(s_orientation),1))
            #n_before_after[(s_orientation<0).flatten(),:] = np.matlib.repmat(np.array([surfacelist[n].n_inside, surfacelist[n].n_outside]),np.sum((s_orientation<0)),1) #cast bool array to int
            #n_before_after[(s_orientation<0).flatten(),:] = np.tile(np.array([surfacelist[n].n_inside, surfacelist[n].n_outside]),(np.sum((s_orientation<0)),1)) #cast bool array to int
            n_before_after[(s_orientation<0).flatten(),:] = np.array([surfacelist[n].n_inside, surfacelist[n].n_outside])[np.newaxis, :] # is this broadcasting the same as right above?

            abs_before_after = np.tile(np.array([surfacelist[n].abslength_outside, surfacelist[n].abslength_inside]), (np.size(s_orientation),1))
            #abs_before_after[(s_orientation<0).flatten(),:] = np.tile(np.array([surfacelist[n].abslength_inside, surfacelist[n].abslength_outside]),(np.sum((s_orientation<0)),1))
            abs_before_after[(s_orientation<0).flatten(),:] = np.array([surfacelist[n].abslength_inside, surfacelist[n].abslength_outside])[np.newaxis, :]

            scat_before_after = np.tile(np.array([surfacelist[n].rayleigh_outside, surfacelist[n].rayleigh_inside]),(np.size(s_orientation),1))
            #scat_before_after[(s_orientation<0).flatten(),:] = np.tile(np.array([surfacelist[n].rayleigh_inside, surfacelist[n].rayleigh_outside]),(np.sum((s_orientation<0)),1))
            scat_before_after[(s_orientation<0).flatten(),:] = np.array([surfacelist[n].rayleigh_inside, surfacelist[n].rayleigh_outside])[np.newaxis, :]
            
#            % if this is the closest scatter so far, update the scatter
#            % property variables
            scatter_here = np.less(l_ray, l_next)
            l_next[scatter_here] = l_ray[scatter_here]
            s_next[scatter_here,:] = s_normal[scatter_here,:]
            sm_next[scatter_here,:] = s_normal_mech[scatter_here,:]
            p_next[scatter_here,:] = p_intersect[scatter_here,:]
            n_next[scatter_here,:] = n_before_after[scatter_here,:]
            abslength_next[scatter_here,:] = abs_before_after[scatter_here,:]
            rayleigh_next[scatter_here,:] = scat_before_after[scatter_here,:]
            abs_next[scatter_here] = surfacelist[n].absorption
            six_next[scatter_here] = n * s_orientation[scatter_here[:, np.newaxis]] # Turned scatter_here from n,  to n,1
            surfacetype_next[scatter_here] = int_surfacetype
            unifiedsurface_next[scatter_here,:] = np.tile(np.array(surfacelist[n].unifiedparams), (np.sum(scatter_here), 1))
        #end loop through surfaces
        
        

        
#        %% refigure s_normal for rays hitting weird surfaces like diffuse
        diffuse_cut = (surfacetype_next==1)
        if np.any(diffuse_cut): # REMOVE REPMATS
            n_diffuse = np.sum(int(diffuse_cut)) #cast to int

            cos_theta = np.sqrt(np.random.randn(n_diffuse,1)) #% diffuse ~= isotropic
            sin_theta = np.sqrt(1-cos_theta**2)
            phi = np.random.randn(n_diffuse,1) * 2 * np.pi

            x_tmp = np.cross(s_next[diffuse_cut,:], np.matlib.repmat([1, 0, 0], n_diffuse,1), axis=1) #set axis... numpy cross() inputs differ from matlab
            y_tmp = np.cross(s_next[diffuse_cut,:], np.matlib.repmat([0, 1, 0],n_diffuse,1), axis=1) #set axis... numpy cross() inputs differ from matlab
            tmpcut = np.all(x_tmp==0,2)
            x_tmp[tmpcut,:] = y_tmp[tmpcut,:]
            x_tmp = x_tmp / np.matlib.repmat(np.abs(np.sqrt(np.sum(x_tmp**2, 1))), 1, 3)
            y_tmp = np.cross(s_next[diffuse_cut,:], x_tmp)

            outdir = s_next[diffuse_cut,:]*np.matlib.repmat(cos_theta,1,3) + x_tmp*np.matlib.repmat(sin_theta * np.cos(phi),1,3) + y_tmp*np.matlib.repmat(sin_theta * np.sin(phi),1,3)

            s_normal_tmp = outdir - incoming_rays[diffuse_cut,0:3]
            s_next[diffuse_cut,:] = s_normal_tmp / np.matlib.repmat(np.abs(np.sqrt(np.sum(s_normal_tmp**2, 1))), 1, 3)
            

#        %% Now determine which rays get absorbed, rayleigh scattered, or surface scattered
        scatter_cut = np.logical_or((l_next < np.inf), (rayleigh_next[:,0] < np.inf))
        if not np.any(scatter_cut):
            ray_index = []
            continue #move to next scatter
        
        l_to_bulkscatter = -rayleigh_next[:,0] * np.log(np.random.rand(rayleigh_next.shape[0]))

        surface_scatter_cut = np.logical_and(scatter_cut, (np.less_equal(l_next, l_to_bulkscatter)))
        unified_scatter_cut = np.logical_and(surface_scatter_cut, (surfacetype_next == 2))
        normal_scatter_cut = np.logical_and(surface_scatter_cut, ~unified_scatter_cut)
        rayleigh_scatter_cut = np.logical_and(scatter_cut, ~surface_scatter_cut)

#        %% Now implement rayleigh scatters (scattering bit will happen later)
        smix_next = np.copy(six_next)
        if np.any(rayleigh_scatter_cut):
            six_next[rayleigh_scatter_cut] = 0
            l_next[rayleigh_scatter_cut] = l_to_bulkscatter[rayleigh_scatter_cut]
            p_next[rayleigh_scatter_cut, :] = p_start[rayleigh_scatter_cut, :] + (l_to_bulkscatter[rayleigh_scatter_cut, np.newaxis] * incoming_rays[rayleigh_scatter_cut, 0:3]) # removed repmat
        
#        %% Now apply bulk absorption
        trans_frac = np.exp(-l_next / abslength_next[:,0])
        incoming_intensity = np.copy(incoming_rays[:,6])
        bulk_abs = incoming_intensity * (1 - trans_frac)
        incoming_rays[scatter_cut, 6:10] = incoming_rays[scatter_cut, 6:10] * trans_frac[scatter_cut, np.newaxis] # removed (1,4) repmat around trans_frac
        
#        %% now initialize the refracted and reflected ray lists | INCOMING_RAYS IS CHANGED | fixed by copying instead of pointing to
        refracted_rays = np.copy(incoming_rays)
        refracted_rays[:, 6:10] = 0
        reflected_rays = np.copy(incoming_rays)
        reflected_rays[:, 6:10] = 0

#        %% Now handle the scattering
#        % First handle normal, diffuse, and retro surfaces
#        % (all subject to 'normal_scatter_cut')
        if np.any(normal_scatter_cut): # NaN's are happening here
            [refracted_rays[normal_scatter_cut,:], reflected_rays[normal_scatter_cut,:]] = RefractionReflectionAtInterface.RefractionReflectionAtInterface(incoming_rays[normal_scatter_cut,:], s_next[normal_scatter_cut,:], n_next[normal_scatter_cut,0], n_next[normal_scatter_cut,1], tir_handling)
            ######## BOOKMARK ##########

#        % Next handle unified reflecting surfaces
        if np.any(unified_scatter_cut):
            reflected_rays[unified_scatter_cut,:] = UnifiedReflectorModel.UnifiedReflectorModel(incoming_rays[unified_scatter_cut,:], sm_next[unified_scatter_cut,:], n_next[unified_scatter_cut,0], n_next[unified_scatter_cut,1], unifiedsurface_next[unified_scatter_cut,:])
        
#        % apply the absorption coefficient for all surfaces
        if np.any(surface_scatter_cut): # ISSUE IS WITH refracted/reflected_rays[surface_scatter_cut,6:10] TURNING NaN
            refracted_rays[surface_scatter_cut,6:10] = refracted_rays[surface_scatter_cut,6:10] * (1-abs_next[surface_scatter_cut])[:,np.newaxis] # removed (1,4) repmats
            reflected_rays[surface_scatter_cut,6:10] = reflected_rays[surface_scatter_cut,6:10] * (1-abs_next[surface_scatter_cut])[:,np.newaxis]
        
#        % and finally do the Rayleigh-scattered rays
        rsc = RayleighScatteringClass.RayleighScatteringClass()
        if np.any(rayleigh_scatter_cut):
            reflected_rays[rayleigh_scatter_cut, :] = rsc.RayleighScattering(incoming_rays[rayleigh_scatter_cut, :])
            
#        % now cut the raylist in half if singlechild is true --
#        % decide reflection vs refraction by dice roll, and call it refraction
        if singlechild:
            total_amp = reflected_rays[:, 6] + refracted_rays[:, 6]
            # print("refracted: " + str(refracted_rays[:, 6]))
            # print("reflected: " + str(reflected_rays[:, 6]))
            reflection_roll = np.less(np.random.rand(reflected_rays.shape[0]), (reflected_rays[:, 6] / total_amp))
            refracted_rays[reflection_roll, :] = reflected_rays[reflection_roll, :]
            amp_rescale = total_amp / refracted_rays[:, 6]
            amp_rescale[np.isnan(amp_rescale)] = 0
            #print("total amp: " + str(total_amp))
            total_amp[np.isnan(total_amp)] = 0
            refracted_rays[:, 6] = total_amp
            refracted_rays[:, 7:10] = refracted_rays[:, 7:10] * amp_rescale[:,np.newaxis] # removed repmat
            reflected_rays[:, 6:10] = 0

            
        surface_abs = incoming_rays[:, 6] - refracted_rays[:, 6] - reflected_rays[:, 6]
        # print("incoming abs: " + str(incoming_rays[:, 6]))
        # print("refracted abs: " + str(refracted_rays[:, 6]))
        # print("reflected abs: " + str(reflected_rays[:, 6]))
        # print("surface absorption: " + str(surface_abs))

#        %% now fill out the first three pieces of the absorption table
        for i_s in range(1,len(surfacelist)):
            inward_cut = smix_next == i_s
            outward_cut = smix_next == -i_s
            infrom_cut = smix_last == -i_s
            outfrom_cut = smix_last == i_s


            absorption_table[num_scatters-1, 0, i_s, 0] = np.sum(surface_abs[np.logical_and(surface_scatter_cut, inward_cut)])
            absorption_table[num_scatters-1, 0, i_s, 1] = np.sum(surface_abs[np.logical_and(surface_scatter_cut, outward_cut)])
            absorption_table[num_scatters-1, 1, i_s, 0] = np.sum(bulk_abs[np.logical_and(scatter_cut, inward_cut)])
            absorption_table[num_scatters-1, 1, i_s, 1] = np.sum(bulk_abs[np.logical_and(scatter_cut, outward_cut)])
            absorption_table[num_scatters-1, 2, i_s, 0] = np.sum(incoming_intensity[np.logical_and(~scatter_cut, outfrom_cut)])
            absorption_table[num_scatters-1, 2, i_s, 1] = np.sum(incoming_intensity[np.logical_and(~scatter_cut, infrom_cut)])
            
#        %% store output
        if full_output:
#            % we round the indices because the sign command used to find
#            % orientation sometimes gives non-integer surface indices
            curr_ray_interface = rayInterfaces.rayInterfaces()
            curr_ray_interface.incoming_ray = incoming_rays[scatter_cut, :]
            curr_ray_interface.refracted_ray = refracted_rays[scatter_cut, :]
            curr_ray_interface.reflected_ray = reflected_rays[scatter_cut, :]
            curr_ray_interface.intersection_point = p_next[scatter_cut, :]
            curr_ray_interface.surface_normal = sm_next[scatter_cut, :]
            curr_ray_interface.ray_index = np.round_(ray_index[scatter_cut])
            curr_ray_interface.surface_index = np.round_(six_next[scatter_cut])
            curr_ray_interface.distance_traveled = l_next[scatter_cut]
            curr_ray_interface.n_incident = n_next[scatter_cut, 0]
            curr_ray_interface.n_transmitted = n_next[scatter_cut, 1]
            curr_ray_interface.bulkabs_incident = abslength_next[scatter_cut, 0]
            curr_ray_interface.bulkabs_transmitted = abslength_next[scatter_cut, 1]
            curr_ray_interface.rayleigh_incident = rayleigh_next[scatter_cut, 0]
            curr_ray_interface.rayleigh_transmitted = rayleigh_next[scatter_cut, 1]
            ray_interfaces.append(curr_ray_interface)

        
        if output_raytable: # check num_scatters indexing
            raytable_cut = np.logical_and(scatter_cut, (ray_index > 0))
            ray_ix = np.around(ray_index[raytable_cut])
            raytable[num_scatters, ray_ix, 0:3] = p_next[raytable_cut, :]
            raytable[num_scatters, ray_ix, 3:13] = refracted_rays[raytable_cut, :]

#        %% get set for next iteration
#        % follow reflected and refracted rays that are above the follow_threshold
        refracted_rays_to_follow = np.logical_and(scatter_cut, (refracted_rays[:,6] > follow_threshold[0]))
        reflected_rays_to_follow = np.logical_and(scatter_cut, (reflected_rays[:,6] > follow_threshold[1]))
                
        for i_s in range(1,np.size(surfacelist)):
            inward_cut = smix_next == i_s
            outward_cut = smix_next == -i_s

            absorption_table[num_scatters-1, 3, i_s, 0] = np.sum(refracted_rays[np.logical_and.reduce((~refracted_rays_to_follow, scatter_cut, inward_cut)), 6]) + np.sum(reflected_rays[np.logical_and.reduce((~reflected_rays_to_follow, scatter_cut, inward_cut)), 6])
            absorption_table[num_scatters-1, 3, i_s, 1] = np.sum(refracted_rays[np.logical_and.reduce((~refracted_rays_to_follow, scatter_cut, outward_cut)), 6]) + np.sum(reflected_rays[np.logical_and.reduce((~reflected_rays_to_follow, scatter_cut, outward_cut)), 6])
            absorption_table[num_scatters-1, 4, i_s, 0] = np.sum(refracted_rays[np.logical_and(refracted_rays_to_follow, inward_cut), 6]) + np.sum(reflected_rays[np.logical_and(reflected_rays_to_follow, inward_cut), 6])
            absorption_table[num_scatters-1, 4, i_s, 1] = np.sum(refracted_rays[np.logical_and(refracted_rays_to_follow, outward_cut), 6]) + np.sum(reflected_rays[np.logical_and(reflected_rays_to_follow, outward_cut), 6])

        #turn the following into nparrays to preserve 2-d shape from matlab code | must concatenate w/ numpy
        p_start = np.concatenate([p_next[refracted_rays_to_follow,:], p_next[reflected_rays_to_follow,:]])
        incoming_rays = np.concatenate([refracted_rays[refracted_rays_to_follow,:], reflected_rays[reflected_rays_to_follow,:]]) # NO REFLECTED RAYS TO FOLLOW, STILL HAVE SOME NANs
#        % smix_last identifies the volume the ray is traversing next, only used
#        % for rays that escape the geometry in the next step
        smix_last = np.concatenate([-smix_next[refracted_rays_to_follow], smix_next[reflected_rays_to_follow]]).flatten()
        six_last = np.abs(np.concatenate([six_next[refracted_rays_to_follow], six_next[reflected_rays_to_follow]]))[:,np.newaxis]   #Changed from np.abs(np.array([six_next(refracted_rays_to_follow), six_next(reflected_rays_to_follow)]))    ; note parantheses
#        % identify reflected rays with a negative ray_index (refracted rays
#        % also inhered the negative index if they have previously been
#        % reflected)
        ray_index = np.concatenate([ray_index[refracted_rays_to_follow], -np.abs(ray_index[reflected_rays_to_follow])]).flatten()
        
    absorption_table = absorption_table[0:num_scatters, :, :, :]
    absorption_table = np.delete(absorption_table,0,axis=2)  # Deletes the null surface that we made to avoid index -0 issues.

    return [np.array(ray_interfaces), absorption_table, raytable]

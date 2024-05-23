# This cell establishes a class based on the output of RayTracer2 that allows for a bit more of a user friendly display

# Boilderplate
import pandas as pd
import numpy as np
import RayTracer2

class RayTracer2_Display:
    '''This class takes the output data from RayTracer2 and displays it via various methods within this class to make 
    the data a bit easier to follow.
    
    Use the .methods() method to see a list of methods pertaining to displaying the output of RayTracer2
    Use the .descriptions() method to see the more in depth original descriptions from RayTracer2 MatLab version'''
    
    def __init__(self,starting_points,rays,surfacelist,tir_handling=[],max_scatters=10,\
                 min_travel_length=np.spacing(np.float64(1)),follow_thresh = np.array([0, 0], dtype=np.float64),\
                 full_output=True,singlechild=True,output_raytable=False ):
        
        '''starting_points:     N-by-3 array, where N is the number of
                                    initial rays to follow, giving the
                                    starting point for each ray.
        
           rays:                 N-by-10 array giving the initial
                                    direction, intensity, and polarization
                                    of each ray.  The first 3 columns give
                                    the forward direction of the ray
                                    (these will be normalized if they
                                    aren't already), columns 4-6 give a
                                    direction non-parallel to the ray that
                                    defines the s1 polarization axis
                                    (these will be made normal to the ray
                                    direction and normalized, if they
                                    aren't already), and columns 7-10 are
                                    the stokes parameters s0-s3, giving
                                    the intensity and polarization of the
                                    ray ( s0 gives the total intensity,
                                    and s0^2 >= s1^2 + s2^2 + s3^2, see
                                    7.2 in Jackson for more details) 
                                    
            surfacelist:           A class based list defining the geometry
                                    of scattering surfaces
                                    
           max_scatters        -  The maximum number of scatters to
                                   propagate rays through (the simulation
                                   may stop before this point if there
                                   are no rays above threshold left to
                                   follow) Default is 10.
          min_travel_length   -  A minimum travel length between
                                   scatters.  This prevents rounding
                                   errors from causing a ray to scatter
                                   multiple times at the same point.
                                   Rays with legitimate travel lengths
                                   below this value will be INCORRECTLY
                                   RECONSTRUCTED, so keep this small
                                   (~1e-5 times your smallest dimension
                                   is probably sufficient for most
                                   geometries.) Default is 2.220446049250313e-16 
          follow_threshold    -  Refracted or reflected rays with an s0
                                   below follow_threshold(1) or
                                   follow_threshold(2), respectively,
                                   will not be followed.  If
                                   follow_threshold is a scalar, the same
                                   threshold is used for both. Default is 0
          tir_handling        -   This determines what the refracted_rays
                                    output is in the case of total internal
                                    reflection.  The default (-1) gives a
                                    refracted ray tangent to the surface with
                                    zero intensity.  Any value >=0 will give
                                    a ray with the same direction and
                                    polarization as the reflected ray, with
                                    intensity equal to the reflected
                                    intensity times tir_handling.  This lets
                                    you treat tir-rays like refracted rays,
                                    which can be handy in geometry sims.
                                    NOTE -- if follow_threshold(2) is
                                    bigger than max(rays(:,7)) then
                                    default tir_handling=1.
          full_output         -   If false, the ray_interfaces output is
                                    not populated.  Default is true.
          singlechild         -   If true, then reflected/refracted rays
                                    are never both followed, rather one
                                    is chosen by a dice roll, and the ray
                                    index is always positive.  If false,
                                    this follows the old RayTracer
                                    standard of following both reflected
                                    and refracted rays, with purely
                                    refracted trajectories only getting
                                    the positive index.  Default is true.
          output_raytable     -   If false, the raytable output is not
                                    populated.  Default is false.'''
        
        self.output = RayTracer2.RayTracer2(starting_points,rays,surfacelist,output_raytable=True)
    
    # creates a method that lists all available display commands
    def methods(self):
        '''Lists all of the methods, along with short descriptions, that relate to displaying the RayTracer2 output.'''
        
        
        # list of commands
        ray_interfaces_indices = ['incoming_ray','reflected_ray','refracted_ray','intersection_point','surface_normal',\
                                  'ray_index','surface_index','distance_traveled','n_incident','n_transmitted', \
                                  'bulkabs_incident','bulkabs_transmitted','rayleigh_incident','rayleigh_transmitted'\
                                  ,'absorption_table','raytable']
        
        # descriptions for those commands
        descriptions = ['[M x 11 array] M = number of scattered rays','Lists direction, intensity, and polarization of \
        all scattered rays.'],['[M x 11 array]','Lists direction, itensity, and polarization of reflected rays.'],\
        ['[M x 11 array]','Lists direction, itensity, and polarization of refracted rays.'],\
        ['[M x 4 array]','Lists coordinates of the scatter.'],\
        ['[M x 4 array]','Lists the backward pointing surface normal at intersection point.'],\
        ['[M x 2 array]','Ray index (negative if it reflected at least once).'],\
         ['[M x 2 array]','Surface where scattering occurred.'],\
         ['[M x 2 array]','Distance travelled since last scatter '],\
         ['[M x 2 array]','Index of refraction for incoming rays.'],\
         ['[M x 2 array]','Index of refraction for refracted rays.'],\
         ['[M x 2 array]','Lists absorption length of incoming rays.'],\
         ['[M x 2 array]','Lists absorption length of refracted rays.'],\
         ['[M x 2 array]','Lists Rayleigh scattering length for incoming rays.'],\
         ['[M x 2 array]','Lists the Rayleigh scattering length for refracted rays.'],\
         ['[K x 5 x S x 2 array] K=scatters, S=surfaces','Lists various absorption data.'],\
        ['[K+1 x N x 13 array]', 'Follows all scatters of each ray (does not follow reflections).']
        
        # Creates and formats the data frame
        ray_interfaces=self.output[0]
        methods_raw=pd.DataFrame(descriptions, index=pd.Index(ray_interfaces_indices, name='Method'),\
                              columns=['Output','Description'])
        d=dict(selector="th",props=[('text-align', 'center')])
        methods = methods_raw.style.set_properties().set_table_styles([d]).set_caption('Each listed method can be called\
        via the following line of code: self.{method}().  Keep in mind the "self" will be replaced by whatever name you gave\
        to the variable referencing this class.')
        return  methods

        
        
        
        
        
### The display methods ###
      
    
    
    
    def incoming_ray(self,rows=40,fancy=True):
        '''Creates a data frame listing the direction, intensity, polarization, and scatter number of the incoming rays
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
       % M-by-11 data frame, where M is the number of rays scattering in this iteration, giving the direction, intensity, 
       % polarization, and scatter number of the incoming rays. Rays that do not scatter are not reported (to report all rays, 
       % enclose your geometry in an absorbing box, for example). 
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        incoming_ray_list = self.output[0][0].incoming_ray
        incoming_ray_list = np.append(incoming_ray_list,np.transpose([np.zeros(np.shape(self.output[0][0].incoming_ray)[0])]),\
                                      axis=1) 
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].incoming_ray
            # The line below adds the scatter number column
            data = np.append(data,np.transpose([np.ones(np.shape(self.output[0][i].incoming_ray)[0])])*[i],axis=1)
            incoming_ray_list = np.append(incoming_ray_list,data,axis=0)
        
        # Creates the data frame
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            incoming_ray_fancy=pd.DataFrame(incoming_ray_list,columns=pd.MultiIndex((['Direction Vector','Polarization Axis', \
                    'Stokes Parameters',''],['x_d','y_d','z_d','x_p','y_p','z_p','s0','s1','s2','s3','scatter']),\
                    ([0,0,0,1,1,1,2,2,2,2,3],[0,1,2,3,4,5,6,7,8,9,10])))
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            incoming_ray_df = incoming_ray_fancy.head(rows).style.set_table_styles([d]).format({('','scatter') : "{:.0f}"})\
            .set_caption('Incoming Rays Table')
            return incoming_ray_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            incoming_ray_raw=pd.DataFrame(incoming_ray_list,columns=['x_d','y_d','z_d','x_p','y_p','z_p',\
                                                                       's0','s1','s2','s3','scatter'])
            return incoming_ray_raw
    
    
    
    
    
    
    def reflected_ray(self,rows=40,fancy=True):
        '''Creates a data frame listing the direction, intensity, polarization, and scatter number of the reflected rays
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
       % M-by-11 data frame, where M is the number of rays scattering in this iteration, giving the direction, intensity, 
       % polarization, and scatter number of the reflected rays. Rays that do not scatter are not reported (to report all rays, 
       % enclose your geometry in an absorbing box, for example). 
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        reflected_ray_list = self.output[0][0].reflected_ray
        reflected_ray_list = np.append(reflected_ray_list,np.transpose\
                                       ([np.zeros(np.shape(self.output[0][0].reflected_ray)[0])]), axis=1) 
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].reflected_ray
            # The line below adds the scatter number column
            data = np.append(data,np.transpose([np.ones(np.shape(self.output[0][i].reflected_ray)[0])])*[i],axis=1)
            reflected_ray_list = np.append(reflected_ray_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            reflected_ray_fancy=pd.DataFrame(reflected_ray_list,columns=pd.MultiIndex((['Direction Vector',\
            'Polarization Axis','Stokes Parameters',''],['x_d','y_d','z_d','x_p','y_p','z_p','s0','s1','s2','s3','scatter']),\
            ([0,0,0,1,1,1,2,2,2,2,3],[0,1,2,3,4,5,6,7,8,9,10])))
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            reflected_ray_df = reflected_ray_fancy.head(rows).style.set_table_styles([d]).format({('','scatter') : "{:.0f}"})\
            .set_caption('Reflected Rays Table')
            return reflected_ray_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            reflected_ray_raw=pd.DataFrame(reflected_ray_list,columns=['x_d','y_d','z_d','x_p','y_p','z_p',\
                                                                       's0','s1','s2','s3','scatter'])
            return reflected_ray_raw
        
        
        
        
        
    def refracted_ray(self,rows=40,fancy=True):
        '''Creates a data frame listing the direction, intensity, polarization, and scatter number of the refracted rays
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
       % M-by-11 data frame, where M is the number of rays scattering in this iteration, giving the direction, intensity, 
       % polarization, and scatter number of the refracted rays. Rays that do not scatter are not reported (to report all rays, 
       % enclose your geometry in an absorbing box, for example). 
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        refracted_ray_list = self.output[0][0].refracted_ray
        refracted_ray_list = np.append(refracted_ray_list,np.transpose\
                                       ([np.zeros(np.shape(self.output[0][0].refracted_ray)[0])]), axis=1) 
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].refracted_ray
            # The line below adds the scatter number column
            data = np.append(data,np.transpose([np.ones(np.shape(self.output[0][i].refracted_ray)[0])])*[i],axis=1)
            refracted_ray_list = np.append(refracted_ray_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            refracted_ray_fancy=pd.DataFrame(refracted_ray_list,columns=pd.MultiIndex((['Direction Vector',\
            'Polarization Axis','Stokes Parameters',''],['x_d','y_d','z_d','x_p','y_p','z_p','s0','s1','s2','s3','scatter']),\
            ([0,0,0,1,1,1,2,2,2,2,3],[0,1,2,3,4,5,6,7,8,9,10])))
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            refracted_ray_df = refracted_ray_fancy.head(rows).style.set_table_styles([d]).format({('','scatter') : "{:.0f}"})\
            .set_caption('Refracted Rays Table')
            return refracted_ray_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            refracted_ray_raw=pd.DataFrame(refracted_ray_list,columns=['x_d','y_d','z_d','x_p','y_p','z_p',\
                                                                       's0','s1','s2','s3','scatter'])
            return refracted_ray_raw
        
 




    def intersection_point(self,rows=40,fancy=True):
        '''Creates a data frame listing the x,y,z coordinates of each scatter
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
        % M-by-4 data frame, where M is the number of rays scattering in this iteration,giving the points where
        % the incoming rays scatter and the scatter number
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        intersection_point_list = self.output[0][0].intersection_point
        intersection_point_list = np.append(intersection_point_list,np.transpose\
                                       ([np.zeros(np.shape(self.output[0][0].intersection_point)[0])]), axis=1) 
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].intersection_point
            # The line below adds the scatter number column
            data = np.append(data,np.transpose([np.ones(np.shape(self.output[0][i].intersection_point)[0])])*[i],axis=1)
            intersection_point_list = np.append(intersection_point_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            intersection_point_fancy=pd.DataFrame(intersection_point_list,columns=['x','y','z','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            intersection_point_df = intersection_point_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                                    : "{:.0f}"}).set_caption('Intersection Points Table')
            return intersection_point_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            intersection_point_raw=pd.DataFrame(intersection_point_list,columns=['x','y','z','scatter'])
            return intersection_point_raw
        
        
        
        
        
    def surface_normal(self,rows=40,fancy=True):
        '''Creates a data frame listing the x,y,z coordinate vector for the backward facing surface normal at the intersection 
            point.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
        % M-by-4 data frame, where M is the number of rays scattering in this iteration,giving the points where
        % the incoming rays scatter and the scatter number
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        surface_normal_list = self.output[0][0].surface_normal
        surface_normal_list = np.append(surface_normal_list,np.transpose\
                                       ([np.zeros(np.shape(self.output[0][0].surface_normal)[0])]), axis=1) 
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].surface_normal
            # The line below adds the scatter number column
            data = np.append(data,np.transpose([np.ones(np.shape(self.output[0][i].surface_normal)[0])])*[i],axis=1)
            surface_normal_list = np.append(surface_normal_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            surface_normal_fancy=pd.DataFrame(surface_normal_list,columns=['x','y','z','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            surface_normal_df = surface_normal_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                                                    : "{:.0f}"}).set_caption('Surface Normal Table')
            return surface_normal_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            surface_normal_raw=pd.DataFrame(surface_normal_list,columns=['x','y','z','scatter'])
            return surface_normal_raw
        
        
        
        
        
    def ray_index(self,rows=40,fancy=True):
        '''Creates a data frame listing the ray indices, which negative indices representing reflections
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
        % M-by-2 array, giving the index of the incoming ray (the input rays are numbered 1:N) -- a negative 
            ray_index means the ray has undergone at least one reflection in its history, so there will be at most one ray 
            with a given positive index. Also lists scatter number
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        ray_index_list = self.output[0][0].ray_index
        # have to add another step compared to the other methods because the data is only one dimension right now
        ray_index_list = np.append([ray_index_list],[np.zeros(np.shape(self.output[0][0].ray_index)[0])], axis=0) 
        ray_index_list = np.transpose(ray_index_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].ray_index
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].ray_index)[0])*[i]],axis=0)
            data = np.transpose(data)
            ray_index_list = np.append(ray_index_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            ray_index_fancy=pd.DataFrame(ray_index_list,columns=['index','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            ray_index_df = ray_index_fancy.head(rows).style.set_table_styles([d]).format("{:.0f}").set_caption('Ray Index Table')
            return ray_index_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            ray_index_raw=pd.DataFrame(ray_index_list,columns=['index','scatter'])
            return ray_index_raw 
        
        
        
        
        
        
    def surface_index(self,rows=40,fancy=True):
        '''Creates a data frame listing the surface indices of each scatter.
        
        NOTE: Unlike usual indexing, index 1 refers to the first surface in your list (index 0 does not, issue
               with having a -0 required some changes.)
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
       % M-by-2 array, giving the index of the incoming ray (the input rays are numbered 1:N) -- a negative ray_index means 
                the ray has undergone at least one reflection in its history, so there will be at most one ray with a given 
                positive index. Also gives scatter number.
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        surface_index_list = self.output[0][0].surface_index
        # have to add another step compared to the other methods because the data is only one dimension right now
        surface_index_list = np.append([surface_index_list],[np.zeros(np.shape(self.output[0][0].surface_index)[0])], axis=0) 
        surface_index_list = np.transpose(surface_index_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].surface_index
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].surface_index)[0])*[i]],axis=0)
            data = np.transpose(data)
            surface_index_list = np.append(surface_index_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            surface_index_fancy=pd.DataFrame(surface_index_list,columns=['surface','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            surface_index_df = surface_index_fancy.head(rows).style.set_table_styles([d]).format("{:.0f}")\
            .set_caption('Surface Index Table')
            return surface_index_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            surface_index_raw=pd.DataFrame(surface_index_list,columns=['surface','scatter'])
            return surface_index_raw
        
        
        
        
    def distance_traveled(self,rows=40,fancy=True):
        '''Creates a data frame listing the distance traveled of each ray since its last scatter.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
        % M x 2 array. Includes scatter numbers
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        distance_traveled_list = self.output[0][0].distance_traveled
        # have to add another step compared to the other methods because the data is only one dimension right now
        distance_traveled_list = np.append([distance_traveled_list],[np.zeros(np.shape(self.output[0][0]\
                                                                                       .distance_traveled)[0])], axis=0) 
        distance_traveled_list = np.transpose(distance_traveled_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].distance_traveled
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].distance_traveled)[0])*[i]],axis=0)
            data = np.transpose(data)
            distance_traveled_list = np.append(distance_traveled_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            distance_traveled_fancy=pd.DataFrame(distance_traveled_list,columns=['distance','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            distance_traveled_df = distance_traveled_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                                                    : "{:.0f}"}).set_caption('Distance Traveled Table')
            return distance_traveled_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            distance_traveled_raw=pd.DataFrame(distance_traveled_list,columns=['distance','scatter'])
            return distance_traveled_raw
        
        
        
        
        
    def n_incident(self,rows=40,fancy=True):
        '''Creates a data frame listing the index of refraction for the incoming ray.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
                
        % M x 2 array. Includes scatter numbers.
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        n_incident_list = self.output[0][0].n_incident
        # have to add another step compared to the other methods because the data is only one dimension right now
        n_incident_list = np.append([n_incident_list],[np.zeros(np.shape(self.output[0][0].n_incident)[0])], axis=0) 
        n_incident_list = np.transpose(n_incident_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].n_incident
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].n_incident)[0])*[i]],axis=0)
            data = np.transpose(data)
            n_incident_list = np.append(n_incident_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            n_incident_fancy=pd.DataFrame(n_incident_list,columns=['refraction index','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            n_incident_df = n_incident_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                                                    : "{:.0f}"}).set_caption('Index of Refraction (n_incident)\
                                                                       Table')
            return n_incident_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            n_incident_raw=pd.DataFrame(n_incident_list,columns=['refraction index','scatter'])
            return n_incident_raw
        
        
        
    def n_transmitted(self,rows=40,fancy=True):
        '''Creates a data frame listing the index of refraction for the refracted ray.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
                
        % M x 2 array. Includes scatter numbers.
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        n_transmitted_list = self.output[0][0].n_transmitted
        # have to add another step compared to the other methods because the data is only one dimension right now
        n_transmitted_list = np.append([n_transmitted_list],[np.zeros(np.shape(self.output[0][0].n_transmitted)[0])], axis=0) 
        n_transmitted_list = np.transpose(n_transmitted_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].n_transmitted
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].n_transmitted)[0])*[i]],axis=0)
            data = np.transpose(data)
            n_transmitted_list = np.append(n_transmitted_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            n_transmitted_fancy=pd.DataFrame(n_transmitted_list,columns=['refraction index','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            n_transmitted_df = n_transmitted_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                                : "{:.0f}"}).set_caption('Index of Refraction (n_transmitted) Table')
            return n_transmitted_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            n_transmitted_raw=pd.DataFrame(n_transmitted_list,columns=['refraction index','scatter'])
            return n_transmitted_raw
        
        
        
    def bulkabs_incident(self,rows=40,fancy=True):
        '''Creates a data frame listing the absorption length of the incoming ray.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
                
        % M x 2 array. Includes scatter numbers.
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        bulkabs_incident_list = self.output[0][0].bulkabs_incident
        # have to add another step compared to the other methods because the data is only one dimension right now
        bulkabs_incident_list = np.append([bulkabs_incident_list],[np.zeros(np.shape(self.output[0][0].bulkabs_incident)[0])],\
                                          axis=0) 
        bulkabs_incident_list = np.transpose(bulkabs_incident_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].bulkabs_incident
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].bulkabs_incident)[0])*[i]],axis=0)
            data = np.transpose(data)
            bulkabs_incident_list = np.append(bulkabs_incident_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            bulkabs_incident_fancy=pd.DataFrame(bulkabs_incident_list,columns=['absorption length','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            bulkabs_incident_df = bulkabs_incident_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                                : "{:.0f}"}).set_caption('Absorption Length (bulkabs_incident) Table')
            return bulkabs_incident_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            bulkabs_incident_raw=pd.DataFrame(bulkabs_incident_list,columns=['absorption length','scatter'])
            return bulkabs_incident_raw
        
        
        
        
    def bulkabs_transmitted(self,rows=40,fancy=True):
        '''Creates a data frame listing the absorption length of the refracted ray.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
                
        % M x 2 array. Includes scatter numbers.
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        bulkabs_transmitted_list = self.output[0][0].bulkabs_transmitted
        # have to add another step compared to the other methods because the data is only one dimension right now
        bulkabs_transmitted_list = np.append([bulkabs_transmitted_list],[np.zeros(np.shape(self.output[0][0]\
                                                                                           .bulkabs_transmitted)[0])],axis=0) 
        bulkabs_transmitted_list = np.transpose(bulkabs_transmitted_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].bulkabs_transmitted
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].bulkabs_transmitted)[0])*[i]],axis=0)
            data = np.transpose(data)
            bulkabs_transmitted_list = np.append(bulkabs_transmitted_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            bulkabs_transmitted_fancy=pd.DataFrame(bulkabs_transmitted_list,columns=['absorption length','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            bulkabs_transmitted_df = bulkabs_transmitted_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                                : "{:.0f}"}).set_caption('Absorption Length (bulkabs_transmitted) Table')
            return bulkabs_transmitted_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            bulkabs_transmitted_raw=pd.DataFrame(bulkabs_transmitted_list,columns=['absorption length','scatter'])
            return bulkabs_transmitted_raw
        
        

               
    def rayleigh_incident(self,rows=40,fancy=True):
        '''Creates a data frame listing the Rayleigh scattering length of the incoming ray.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
                
        % M x 2 array. Includes scatter numbers.
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        rayleigh_incident_list = self.output[0][0].rayleigh_incident
        # have to add another step compared to the other methods because the data is only one dimension right now
        rayleigh_incident_list = np.append([rayleigh_incident_list],[np.zeros(np.shape(self.output[0][0]\
                                                                                           .rayleigh_incident)[0])],axis=0) 
        rayleigh_incident_list = np.transpose(rayleigh_incident_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].rayleigh_incident
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].rayleigh_incident)[0])*[i]],axis=0)
            data = np.transpose(data)
            rayleigh_incident_list = np.append(rayleigh_incident_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            rayleigh_incident_fancy=pd.DataFrame(rayleigh_incident_list,columns=['Rayleigh length','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            rayleigh_incident_df = rayleigh_incident_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                            : "{:.0f}"}).set_caption('Rayleigh Scattering Length (rayleigh_incident) Table')
            return rayleigh_incident_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            rayleigh_incident_raw=pd.DataFrame(rayleigh_incident_list,columns=['Rayleigh length','scatter'])
            return rayleigh_incident_raw 
        
        
        
        
    def rayleigh_transmitted(self,rows=40,fancy=True):
        '''Creates a data frame listing the Rayleigh scattering length of the refracted ray.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
                
        % M x 2 array. Includes scatter numbers.
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        rayleigh_transmitted_list = self.output[0][0].rayleigh_transmitted
        # have to add another step compared to the other methods because the data is only one dimension right now
        rayleigh_transmitted_list = np.append([rayleigh_transmitted_list],[np.zeros(np.shape(self.output[0][0]\
                                                                                           .rayleigh_transmitted)[0])],axis=0) 
        rayleigh_transmitted_list = np.transpose(rayleigh_transmitted_list)
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[0])[0]):
            data = self.output[0][i].rayleigh_transmitted
            # The line below adds the scatter number column
            data = np.append([data],[np.ones(np.shape(self.output[0][i].rayleigh_transmitted)[0])*[i]],axis=0)
            data = np.transpose(data)
            rayleigh_transmitted_list = np.append(rayleigh_transmitted_list,data,axis=0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            rayleigh_transmitted_fancy=pd.DataFrame(rayleigh_transmitted_list,columns=['Rayleigh length','scatter'])
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            rayleigh_transmitted_df = rayleigh_transmitted_fancy.head(rows).style.set_table_styles([d]).format({('scatter')\
                                            : "{:.0f}"}).set_caption('Rayleigh Scattering Length (rayleigh_transmitted) Table')
            return rayleigh_transmitted_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            rayleigh_transmitted_raw=pd.DataFrame(rayleigh_transmitted_list,columns=['Rayleigh length','scatter'])
            return rayleigh_transmitted_raw 
        
    def absorption_table(self,rows=40,fancy=True):
        '''Creates a data frame listing total intensity absorbed separated by scatter number, absorption type, surface related
        to the absorption
        
        Absorption Types:
            1 - Surface absorption
            2 - Bulk absorption
            3 - Escaped geometry
                Listed surface is the last surface hit before escaping
            4 - Dropped below threshold
            5 - Still following
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
                
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        
        # creates the first row of the array, so appending won't throw us an error
        absorption_list = np.array([np.zeros(12)]) 
        
        # loops through the entire absorption table to recreate it in a 2 dimensional data frame with many columns
        for n in range(np.shape(self.output[1])[2]):
            for i in range(np.shape(self.output[1])[0]):
                data = np.array([n+1,i,self.output[1][i][0][n][0],self.output[1][i][0][n][1],self.output[1][i][1][n][0],\
                                 self.output[1][i][1][n][1],self.output[1][i][2][n][0],self.output[1][i][2][n][1],\
                                 self.output[1][i][3][n][0],self.output[1][i][3][n][1],self.output[1][i][4][n][0],\
                                 self.output[1][i][4][n][1]])
                absorption_list = np.append(absorption_list,[data],axis=0)
        
        # delete the place holder row
        absorption_list = np.delete(absorption_list,0,0)
        
        # Creates the data frame
        
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            absorption_fancy=pd.DataFrame(absorption_list, columns=pd.MultiIndex((['','Surface Absorption','Bulk Absorption',\
        'Escaped Geometry','Dropped Below Threshold','Still Following'],['surface','scatter number',\
        'positive_SA','negative_SA','positive_BA','negative_BA','positive_EG',\
        'negative_EG','positive_BT','negative_BT','positive_SF','negative_SF']),\
        ([0,0,1,1,2,2,3,3,4,4,5,5],[0,1,2,3,4,5,6,7,8,9,10,11])))
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            absorption_df = absorption_fancy.round(4).head(rows).style.set_table_styles([d]).set_caption('Absorption Table: \
            Positive refers to absorptions on the postive side of the surface, negative refers to the negative side of the \
            surface').format("{:.0f}")
            
            # delete the place holder row
            return absorption_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            absorption_raw=pd.DataFrame(absorption_list,columns=['surface','scatter number',\
        'positive_SA','negative_SA','positive_BA','negative_BA','positive_EG',\
        'negative_EG','positive_BT','negative_BT','positive_SF','negative_SF'])
            return absorption_raw
        

    def raytable(self,rows=40,fancy=True):
        '''Creates a data frame listing the ray index, direction, intensity, polarization, and scatter number of the incoming 
            rays.
        
        rows = the number of rows to display in the stylized 'fancy' data frame.
        fancy = If True, will show a stylized data frame with limited data calling functions (useful for understanding the data)
                If False, will show the raw data frame (should be used when calling the data for use in other functions)
        
       %  Array of size [K+1, N, 13] giving the details of each ray's path through the geometry, following positive ray indices 
       %  only.  First index is scatter index, starting with initial condition.  Second index is ray index.  For third index, 
       %  columns 1:3 give XYZ position, 4:13 give ray details (direction, intensity, polarization, as in rays input).
                                  
        
        Scatter numbers:
            0 = scatter from starting ray
            1 = scatter from a ray with 1 prior scatter
            2 = scatter from a ray with 2 prior scatters
            3 = etc. '''
        
        # creates the first row of the array, so appending won't throw us an error
        raytable_list = self.output[2][0]
        raytable_list = np.append(raytable_list,np.transpose([np.zeros(np.shape(self.output[2])[1])]),\
                                      axis=1) 
        
        # loops through all scatter iterations to populate one list with 2 dimensions, but adding a scatter number column
        for i in range(1,np.shape(self.output[2])[0]):
            data = self.output[2][i]
            # The line below adds the scatter number column
            data = np.append(data,np.transpose([np.ones(np.shape(self.output[2][i])[0])])*[i],axis=1)
            raytable_list = np.append(raytable_list,data,axis=0)
        
        # Creates the data frame
        # creates a nice looking table, but not very useful for calling data from
        if fancy == True:
            raytable_fancy=pd.DataFrame(raytable_list,columns=pd.MultiIndex((['Position','Direction Vector','Polarization Axis', \
                    'Stokes Parameters',''],['x','y','z','x_d','y_d','z_d','x_p','y_p','z_p','s0','s1','s2','s3','scatter']),\
                    ([0,0,0,1,1,1,2,2,2,3,3,3,3,4],[0,1,2,3,4,5,6,7,8,9,10,11,12,13]),names=['','ray index']))
        # This is all just styling to make it look a little nicer
            d=dict(selector="th",props=[('text-align', 'center')])
            raytable_df = raytable_fancy.head(rows).style.set_table_styles([d]).format({('','scatter') : "{:.0f}",\
                                                            ('','ray index') : "{:.0f}"}).set_caption('Ray Table')
            return raytable_df
        
        #creates the much more useful (for external functions) dataframe
        if fancy == False:
            raytable_raw=pd.DataFrame(raytable_list,columns=['x','y','z','x_d','y_d','z_d','x_p','y_p','z_p','s0','s1','s2'\
                                                             ,'s3','scatter'])
            return raytable_raw
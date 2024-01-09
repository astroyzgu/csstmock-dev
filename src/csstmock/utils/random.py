import numpy  as np 
import healpy as hp 
def mvpoints_vector(vec, from_vec, to_vec):   
    '''
    This function moves points following the rotation from one point (from_vec) to another (to_vec).
    
    Parameters:
    vec: The vector to be moved.
    from_vec: The initial point of rotation.
    to_vec: The final point of rotation.
    
    Returns:
    A transposed vertical stack of the new x, y, and z coordinates.
    '''
    
    # Import the norm function from numpy.linalg
    from numpy.linalg import norm
    
    # Ensure that from_vec, to_vec, and vec are at least 2D arrays
    from_vec = np.atleast_2d(from_vec) 
    to_vec   = np.atleast_2d(to_vec) 
    vec      = np.atleast_2d(vec) 
    
    # Normalize from_vec and to_vec
    from_vec  = from_vec/norm(from_vec, axis = 1)[:, np.newaxis]
    to_vec    = to_vec/norm(to_vec, axis = 1)[:, np.newaxis]
    
    # Calculate the cosine of the angle between from_vec and to_vec
    cs        = np.sum(from_vec*to_vec, axis = 1)
    
    # Calculate the cross product of from_vec and to_vec
    vec_cross = np.cross(from_vec,to_vec, axis = 1)
    
    # Calculate the sine of the angle between from_vec and to_vec
    sn        = norm(vec_cross, axis = 1)
    
    # Normalize vec_cross
    vec_cross  = vec_cross/norm(vec_cross, axis = 1)[:, np.newaxis]
    
    # Decompose vec_cross into x, y, and z components
    x = vec_cross[:,0]; x0 = vec[:,0]
    y = vec_cross[:,1]; y0 = vec[:,1]
    z = vec_cross[:,2]; z0 = vec[:,2]
    
    # Calculate the new x, y, and z coordinates
    x_ = ( cs+(1-cs)*x**2)*x0 + ( (1-cs)*x*y-sn*z )*y0 + ((1-cs)*x*z+sn*y)*z0
    y_ = ((1-cs)*x*y+sn*z)*x0 + (  cs+(1-cs)*y**2 )*y0 + ((1-cs)*y*z-sn*x)*z0
    z_ = ((1-cs)*x*z-sn*y)*x0 + ( (1-cs)*y*z+sn*x )*y0 + ( cs+(1-cs)*z**2)*z0
    
    # Return the new coordinates
    return np.vstack([x_, y_, z_]).T

import numpy  as np 
import healpy as hp  

class sphrand():
    '''
    This class generates random points on a sphere.
    '''
    def generand_lonlat(n, lonra = [-180, 180], latra = [-90, 90], seed = None):
        """
        This function generates random points on a sphere using the Lambert azimuthal equal-area projection.
        
        Parameters:
        n (int): The number of random points to generate.
        lonra (list): The range of longitude in the format [min_longitude, max_longitude].
        latra (list): The range of latitude in the format [min_latitude, max_latitude].
        seed (int): The seed for the random number generator. Default is None.

        Returns:
        tuple: A tuple containing two numpy arrays representing the longitude and latitude of the generated points.

        Maths: 
        Lambert azimuthal equal-area projection -- generate randoms within r [rad]. 
        α = arctan(y/x); r = sqrt(x^2+y^2)
        δ = 2arcsin(r);  
        """
        # If a seed is provided, set the random seed
        if seed is not None: np.random.seed(seed)
        # Convert the range from degree to radians
        lonrange     = np.atleast_1d(lonra) 
        lonrange_rad = lonrange*np.pi/180
        # Calculate the range of r using the sine of half the sum of the latitude and 90 degrees
        latrange     = np.atleast_1d(latra) 
        r_range   = np.sin( 0.5*(90 + latrange)*np.pi/180 )**2  # (0,1)

        # Generate n random points for alpha and r
        alpha     = np.random.uniform(*lonrange_rad, n)
        r         = np.sqrt( np.random.uniform(*r_range, n) )
        # Calculate delta using the arcsine of r
        delta     = 2*np.arcsin(r) - np.pi/2# arcsin==>(0,90)

        # Return the longitude and latitude of the generated points
        return alpha/np.pi*180, delta/np.pi*180

    def generand_cap(n, lon, lat, radius, seed = None): 
        '''
        This function generates random points within a spherical cap.
        
        Parameters:
        n: The number of random points to generate.
        lon: The longitude of the center of the cap.
        lat: The latitude of the center of the cap.
        radius: The radius of the cap.
        seed: The seed for the random number generator.
        
        Returns:
        The longitude and latitude of the random points.
        '''
        
        # Generate random longitudes and latitudes
        alpha, delta = sphrand.generand_lonlat(n, lonra = [-180, 180], latra = [90-radius, 90], seed = seed)
        
        # Convert the longitudes and latitudes to vectors
        vec      = hp.ang2vec(alpha ,delta, lonlat = True) 
        
        # Define the initial and final points of rotation
        from_vec = hp.ang2vec(0, 90, lonlat = True) 
        to_vec   = hp.ang2vec(lon, lat, lonlat = True) 
        
        # Move the points following the rotation from the initial to the final point
        vec      = mvpoints_vector(vec, from_vec, to_vec)
        
        # Convert the vectors back to longitudes and latitudes and return them
        return hp.vec2ang(vec, lonlat = True) 

    def generand_range(n, range = [-180, 180, -90, 90], seed = None):
        """
        This function generates random points on a sphere using the Lambert azimuthal equal-area projection.
    
        Parameters:
        n (int): The number of random points to generate.
        range (list): The range of the sphere in the format [min_longitude, max_longitude, min_latitude, max_latitude].
        seed (int): The seed for the random number generator. Default is None.

        Returns:
        tuple: A tuple containing two numpy arrays representing the longitude and latitude of the generated points.

        Maths: 
        Lambert azimuthal equal-area projection -- generate randoms within r [rad]. 
        α = arctan(y/x); r = sqrt(x^2+y^2)
        δ = 2arcsin(r);  
        """
        return sphrand.generand_lonlat(n, lonra = range[:2], latra = range[2:], seed = seed)
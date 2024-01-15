import numpy as np
from scipy.spatial import KDTree
import healpy as hp 

def friends_of_friends(data, link_length): 
    """
    This function implements the Friends-of-Friends (FoF) algorithm, which is commonly used in astronomy to identify clusters of galaxies.

    Parameters:
    data (numpy.ndarray): A 2D array where each row represents a galaxy and the columns are the coordinates.
    link_length (float): The maximum distance between two galaxies to consider them as linked.

    Returns:
    list: A list of numpy arrays, where each array contains the indices of galaxies in the same group.

    The function first creates a KDTree from the data for efficient spatial queries. 
    Then it finds all galaxies within the `link_length` of each galaxy. 
    The function identifies isolated galaxies (galaxies with no other galaxies within `link_length`) first to speed up the process.
    Then it identifies linked galaxies to build the group catalog. For each galaxy that is not isolated, it finds all galaxies that are linked to it (directly or indirectly), and assigns them the same group ID.
    Finally, the function returns a list of groups, where each group is represented as a list of indices of galaxies in that group.
    """
    kd     = KDTree(data); nparticle = np.shape(data)[0]
    groups  = kd.query_ball_point(data, r = link_length, return_sorted = True)
    ngroups = kd.query_ball_point(data, r = link_length, return_length = True) 
    ngal    = len(groups); 
    igal    = np.arange(ngal) 
    #
    # identify isolated galaxy first to speed up. 
    #
    is_isolated = ngroups == 1
    num_isolated= np.sum(is_isolated) 
    #--- galaxy catalog, assign 0-num_isolated as groupid 
    galcat  = -np.ones(nparticle).astype(int);    
    galcat[is_isolated] = np.arange(num_isolated) 
    #--- group catalog1
    grpcat1 = igal[is_isolated]; 
    grpcat1 = list( np.atleast_2d(grpcat1).T )
    #
    # identify linked galaxies to build group catalog2.  
    #
    igrp = num_isolated; grpcat2 = []
    for i in igal[~is_isolated]: 
        # have be assigned into FOF groups
        if     galcat[i] != -1: continue 
        imember     = np.unique( np.hstack([i, groups[i]]) ) 
        new_member  = 1
        # print(igrp, imember) 
        while new_member != 0: 
            new_imember = np.unique( np.hstack(groups[imember]) )
            #-- print('#--', igrp, new_imember) 
            new_member  = len(new_imember) - len(imember)
            imember     = new_imember 
        grpcat2.append(imember) 
        galcat[imember] = igrp 
        igrp = igrp + 1 
    #------------------------------------------------------------
    grpcat = grpcat1 + grpcat2
    isort  = np.argsort([ -1.0*len(grpcat_) for grpcat_ in grpcat]) 
    grpcat = [ grpcat[i] for i in isort ]

    galid  = np.hstack(grpcat) 
    grpid  = np.arange(len(grpcat)) 
    num_repeat = np.array([len(grpcat_) for grpcat_ in grpcat])
    grpid  = np.repeat( grpid, num_repeat) 
    grpnum = np.repeat( num_repeat, num_repeat) 

    isort  = np.argsort(galid)
    grpid  = grpid[isort]
    grpnum = grpnum[isort]

    return grpcat, grpid, grpnum

def angular_fof_group(ra, dec, angle): 
    """
    This function identifies groups of close pairs of galaxies based on their angular positions.

    Parameters:
    ra (numpy.ndarray): 1D array of right ascension coordinates of galaxies.
    dec (numpy.ndarray): 1D array of declination coordinates of galaxies.
    angle (float): The maximum angular distance [in degree] between two galaxies to consider them as a close pair.

    Returns:
    list: A list of numpy arrays, where each array contains the indices of galaxies in the same group.

    The function first converts the angular separation to a Euclidean distance, then converts the spherical coordinates (ra, dec) to Cartesian coordinates (x, y, z).
    It then calls the `friends_of_friends` function to identify the groups of galaxies.
    """
    r      = 2.0*np.sin(0.5*angle*np.pi/180) 
    xyz    = hp.ang2vec(ra, dec, lonlat = True) 
    grpcat, grpid, grpnum = friends_of_friends(xyz, r) 
    return grpcat, grpid, grpnum

def close_pair_weight(grpcat, isspec): 
    """
    This function calculates the weight for each galaxy in a group based on whether it has a spectroscopic redshift.

    Parameters:
    grpcat (list): A list of numpy arrays, where each array contains the indices of galaxies in the same group.
    is_spec (numpy.ndarray): 1D boolean array indicating whether each galaxy has a spectroscopic redshift.

    Returns:
    numpy.ndarray: 1D array of weights for each galaxy.

    The function first initializes the weight array with -1.0. Then for each group, it calculates the weight as the fraction of galaxies in the group that have a spectroscopic redshift.
    The weight for galaxies without a spectroscopic redshift is set to 0.0 if they are in a group with at least one galaxy with a spectroscopic redshift, and -1.0 otherwise.
    """
    w = np.hstack(grpcat)*0.0 - 1.0
    for grpcat_ in grpcat: 
        isspec_ = isspec[grpcat_] 
        w_       = 1.0*np.sum(isspec_)/len(isspec_)
        w[grpcat_]    = w_
    indx1 = (w == 0.0)&(~isspec)
    indx2 = (w != 0.0)&(~isspec)
    w[indx1] =  -1  # bad group, shoud remove this region  
    w[indx2] = 0.0  # weighting as 0.0 
    return w


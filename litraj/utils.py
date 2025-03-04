import numpy as np



def get_profile(images, relative = False):
    """
    Calculate energy profile of the trajectory
    
    Parameters
    ----------
    images: list of ase's Atoms
        trajectory
    
    relative: boolean, False by default
        return profile - min(profile)
    
    Returns
    -------
    profile: np.array
        energy profile

    """
    profile = []
    for im in images:
        profile.append(im.get_potential_energy())
    if relative:
        return np.array(profile) - min(profile)
    else:
        return np.array(profile)



def get_barrier(images):
    """
    Calculate migration barrier of the Li-ion migration
    
    Parameters
    ----------
    images: list of ase's Atoms
        trajectory

    Returns
    -------
    barrier: float
        migration barrier

    """
    profile = get_profile(images)
    return max(profile) - min(profile)


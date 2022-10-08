import numpy as np

def reduce_noise (data, cutoff = 5):
    """
    Given data, uses a fourier bandpass filter to reduce noise and estimate uncertainty of each counts per second measurement. 
    
    Arguments:
        * `data` (nx2 numpy array): 2 columns. First is wavelength in angstroms, second is measured counts per second. 
    Keyword arguments: 
        * `cutoff` (float): Value associated with the sensitivity/cutoff value of bandpass filter. Default value: TBD

    Returns:
        * As a 2-tuple:
            - `new_data` (nx2 numpy array): numpy array with the same format as `data`, hopefully with reduced noise
            - `weights` (numpy array length n): weights associated with each new data point.  
    """

    raise NotImplementedError
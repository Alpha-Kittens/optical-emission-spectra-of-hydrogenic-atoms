

def get_max(data, weights):
    """
    Given a list of spectrometer data, returns the max value of the peak and an 1-sigma uncertainty value based on the width of the peak. 

    Arguments:
        * `data` (nx2 numpy array): 2 columns. First is wavelength in angstroms, second is measured counts per second. 
        * `weights` (numpy array length n): weights associated with each data point. May or may not how to use these idk

    Returns: 
        * As a 2-tuple: 
            - `mu`: wavelength (in Angstroms) of the max value
            - `sigma`: 1-sigma uncertainty value associated with width of peak. 
    """

    raise NotImplementedError
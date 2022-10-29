# Imports
from data.data_loader import read_data
import matplotlib.pyplot as plt
import numpy as np


def process_data(data, damping_constant = None, plot_noise_reduction = False, title=None):
    """
    processes data after being intially read from the data file
        1. reduces noise in data and obtains weights
        2. identifies relevant regions in the data
        3. selects the region with the correct peak
        4. slices the data to give only the region we are interested in fitting

    Arguments: 
        * `data` (nx2 numpy array): numpy array (matrix) with 2 columns:
            - First column is x-axis (input to Model)
            - Second column is y-axis (output of Model)

    Returns: 
        * `new_data` (nx2 numpy array): numpy array (matrix) with 2 columns, a slice of the original data
            - First column is x-axis (input to Model)
            - Second column is y-axis (output of Model)
        * 'new_weights' weights for the slice of data of interest
    """

    from noise_reduction import reduce_noise

    # Reducing Noise
    if damping_constant is None and plot_noise_reduction == False:
        new_data, weights = reduce_noise(data)
    elif damping_constant is None:
        if title is None:
            new_data, weights = reduce_noise(data, plot=True)
        else:
            new_data, weights = reduce_noise(data, plot=True, title=title)
    elif plot_noise_reduction == False:
        new_data, weights = reduce_noise(data)
    else:
        if title is None:
            new_data, weights = reduce_noise(data, plot=True)
        else:
            new_data, weights = reduce_noise(data, plot=True, title=title)

    backgrounds, signals = regions(data[:,1])

    maxval = 0
    max_signal = (0,1)
    for signal in signals:
        smax = max(data[:,1][signal[0]:signal[1]+1])
        if smax > maxval:
            maxval = smax
            max_signal = signal
    
    i = np.array(range(len(data[:,1])))

    #return np.logical_and(i >= max_signal[0], i <= max_signal[1]) old return statement
    new_data = data[max_signal[0]:max_signal[1],:]
    new_weights = weights[max_signal[0]:max_signal[1]]
    return new_data, new_weights


def regions(cps, reduce = True):
    """
    Given an array of positive values, determines a cutoff to distinguish signal from noise, then determines regions of indices
    which correspond to background and signal data. 
    
    Arguments:
        * `cps` (np array of positive values): Array for which cutoff is to be chosen
    Keyword Arguments:
        * `reduce` (Boolean): True if regions containing only one point should be removed; False otherwise. 
            Generally should be set to True for wide signals.  

    Returns:
        * As a 2-array
            - `backgrounds` (Array of tuples): Ranges of indices of `cps` classified as background
            - `signals` (Array of tuples):  Ranges of indices of `cps` classified as signal
    """

    cutoff = get_cutoff(cps)
    return get_regions(cps, cutoff, reduce = reduce)




def get_regions(cps, cutoff, reduce = False):
    """
    Given an array of positive values and a cutoff, determines signal and background regions of the data. 
    
    Arguments:
        * `cps` (np array of positive values): Array for which cutoff is to be chosen
        * `cutoff` (float): Classifying value. values of `cps` above `cutoff` are signals; otherwise background.
    Keyword Arguments:
        * `reduce` (Boolean): True if regions containing only one point should be removed; False otherwise. 
            Generally should be set to True for wide signals.  

    Returns:
        * As a 2-array 
            - `backgrounds` (Array of tuples): Ranges of indices of `cps` classified as background
            - `signals` (Array of tuples):  Ranges of indices of `cps` classified as signal
    """

    background_regions = []
    signal_regions = []

    data_below = cps <= cutoff
    
    tracking_background = True
    start = 0

    for i in range(len(cps)):

        if data_below[i] != tracking_background:
            if tracking_background:
                if (not reduce or (i-1) - start > 1) and i-1 >= 0:
                    background_regions.append((start, i-1))
            else:
                if not reduce or (i-1) - start > 1:
                    signal_regions.append((start, i-1))
            tracking_background = not tracking_background
            start = i

    if tracking_background:
        background_regions.append((start, len(cps) - 1))
    else:
        signal_regions.append((start, len(cps) - 1)) 


    return background_regions, signal_regions




# max model, HWHM cutoff
def get_cutoff(cps):
    """
    Given an array of positive values, determines a cutoff do distinguish signal from noise. This algorithm is not particularly refined.
    
    Arguments:
        * `cps` (np array of positive values): Array for which cutoff is to be chosen
    Returns:
        * `cutoff` (float): signal/noise cutoff for `cps` array. 
    """

    hist, binedges = np.histogram(cps, int(max(cps) // min(cps)))

    cutoff = 1 + int(max(cps) // min(cps) / 50)

    #print ("The cutoff is ", binedges[cutoff], ". ", np.sum(hist[cutoff:]), "/", np.sum(hist), " of the data is above the cutoff.")

    return binedges[cutoff]



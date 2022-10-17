import numpy as np
from data.data_loader import read_data
import matplotlib.pyplot as plt

def get_max(data, standard_deviations = 2):
    """
    Given a list of spectrometer data, returns the max value of the peak and an 1-sigma uncertainty value based on the width of the peak. 

    Arguments:
        * `data` (nx2 numpy array): 2 columns. First is wavelength in angstroms, second is measured counts per second. 
        * `standard_deviations` (float): Estimation of number of standard deviations of the mean will be classified as signal. 
            - For wider peaks, this is probably around 2.
            - For narrower peaks, this is probably around 3. 

    Returns: 
        * As a 2-tuple: 
            - `maxes`: Array of wavelengths (in Angstroms) of max values
            - `widths`: Array of 1-sigma uncertainty values associated with width of each peak. 
    """

    #raise NotImplementedError
    
    background, signals = regions(data[:,1])
    print(len(signals))
    print (signals)

    standard_deviations = 2 # For more narrow peaks, may want to change this to 3? 

    maxes = [max(data[signal[0]:signal[1]+1, 1]) for signal in signals]
    widths = [(signal[1] - signal[0])/2/standard_deviations for signal in signals]

    return maxes, widths
    

     

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

    print ("The cutoff is ", binedges[cutoff], ". ", np.sum(hist[cutoff:]), "/", np.sum(hist), " of the data is above the cutoff.")

    return binedges[cutoff]


if __name__ == "__main__":

    data = read_data("data/10.14.22/hydrogen/beta-superfine-1")
    print(get_max(data, 0))
    
    x, y = data[:,0], data[:,1]

    ysize = y.size
    xsize = x.size
    step = (max(x) - min(x) + 1)/xsize
    #print(ysize)
    #print(xsize)
    #print(step)
    from scipy.fft import fft, ifft, fftfreq
    frequencies, amplitudes = fftfreq(ysize, step), fft(y)

    #cps = np.abs(amplitudes)
    #cps = data[:,1]
    #print (max(cps))
    cutoff = get_cutoff(y)
    print (get_regions(y, cutoff, reduce = True))
    plt.plot(x, y, label = "data")
    plt.plot(x, np.where(np.abs(y) > cutoff, np.abs(y), -1), label = "cut data")
    plt.legend()

    amplitudes_new = np.where(np.abs(amplitudes) > cutoff, amplitudes, 0)
    #plt.plot(data[:,0], data[:,1])
    #plt.plot(data[:,0], np.where(cps > cutoff, cps, -1))
    plt.show()
    #hist, binedges = get_cutoffs(data, 0)
    #bincenters = (binedges[1:] + binedges[:-1]) / 2
    #print (bincenters)
    #plt.plot(x, y)
    #plt.plot(x, ifft(amplitudes_new))
    #plt.show()
    #print (hist)
    #print (bincenters)
    
    #plt.bar(bincenters, hist, width = (binedges[-1] - binedges[0]) / len(hist))
    #plt.show()

    #plt.plot()
    #plt.show()
    
import numpy as np
from data.data_loader import read_data
import matplotlib.pyplot as plt

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

def regions(cps, weights, reduce = True):

    cutoff = get_cutoff(cps, weights)
    return get_regions(cps, cutoff, reduce = reduce)

def get_regions(cps, cutoff, reduce = False):

    background_regions = []
    signal_regions = []

    data_below = cps <= cutoff

    #print (data_below)
    
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



def get_cutoff(cps, weights):
    """
    Given an array of positive values and weights, determines a cutoff do distinguish signal from noise. This algorithm is not particularly refined.
    
    Arguments:
        * `cps` (np array of positive values): Array for which cutoff is to be chosen
        * `weights` (np array of same dimension): Weights corresponding to `cps` values. Unused by method. 
    Returns:
        * `cutoff` (float): signal/noise cutoff for `cps` array. 
    """

    #print(np.histogram(cps, max(cps) // min(cps)))

    #print (int(max(cps) // min(cps)))

    print (max(cps))
    print (min(cps))

    print (int(max(cps) // min(cps)))

    hist, binedges = np.histogram(cps, int(max(cps) // min(cps)))

    print (hist, "\n", binedges)

    #print (hist)
    #print (binedges)

    cutoff = 1 + int(max(cps) // min(cps) / 100)

    print ("The cutoff is ", binedges[cutoff], ". ", np.sum(hist[cutoff:]), "/", np.sum(hist), " of the data is above the cutoff.")

    bincenters = (binedges[1:] + binedges[:-1]) / 2
    #plt.bar(bincenters, hist, width = (binedges[-1] - binedges[0]) / len(hist))
    #plt.show()

    return binedges[cutoff]

    #return hist, binedges

if __name__ == "__main__":

    data = read_data("data/DeuteriumScans/SuperFineScans/alpha4H")
    
    x, y = data[:,0], data[:,1]

    ysize = y.size
    xsize = x.size
    step = (max(x) - min(x) + 1)/xsize
    #print(ysize)
    #print(xsize)
    #print(step)
    from scipy.fft import fft, ifft, fftfreq
    frequencies, amplitudes = fftfreq(ysize, step), fft(y)

    cps = np.abs(amplitudes)
    #cps = data[:,1]
    print (max(cps))
    cutoff = get_cutoff(cps, 0)
    print (get_regions(cps, cutoff, reduce = True))
    plt.plot(frequencies, np.abs(amplitudes))
    plt.plot(frequencies, np.where(np.abs(amplitudes) > cutoff, np.abs(amplitudes), -1))

    amplitudes_new = np.where(np.abs(amplitudes) > cutoff, amplitudes, 0)
    #plt.plot(data[:,0], data[:,1])
    #plt.plot(data[:,0], np.where(cps > cutoff, cps, -1))
    plt.show()
    #hist, binedges = get_cutoffs(data, 0)
    #bincenters = (binedges[1:] + binedges[:-1]) / 2
    #print (bincenters)
    plt.plot(x, y)
    plt.plot(x, ifft(amplitudes_new))
    plt.show()
    #print (hist)
    #print (bincenters)
    
    #plt.bar(bincenters, hist, width = (binedges[-1] - binedges[0]) / len(hist))
    #plt.show()

    #plt.plot()
    #plt.show()
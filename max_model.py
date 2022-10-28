import numpy as np
from data.data_loader import read_data
import matplotlib.pyplot as plt
from data_processing import regions
from scipy.fft import fft, ifft, fftfreq


d = lambda data, i : (data[i+1,1] - data[i,1]) / (data[i+1,0] - data[i,0])

def smear(derivatives, i, width = 5):
    min_i = max(i - width, 0)
    max_i = min(i + width, len(derivatives) - 1) + 1 # endpoints don't matter at all
    dists = np.array(range(min_i - i, max_i - i))
    return np.dot(derivatives[min_i : max_i], (1 - (dists / width)**2)) / np.sum(1 - (dists / width)**2)


def hwhm_max(data, weights, noise_reduced = None, plot = True, true_wavelength = None, threshold = 1/2):

    data_half = (max(data[:,1]) + min(data[:,1])) * threshold

    data_above = data[:,1] > data_half
    data_below = data[:,1] <= data_half
    cut_wavelengths = data[data_above,0]
    cut_weights = weights[data_above]
    if noise_reduced is None:
        averaging_weights = data[data_above,1]
    else:
        noise_reduced = np.column_stack((noise_reduced[:,0], [smear(noise_reduced[:,1], i) for i in range(len(noise_reduced[:,1]))]))
        bumpy_derivatives = [d(noise_reduced, i) for i in range(len(noise_reduced[:,0]) - 1)]
        derivatives = np.array([smear(bumpy_derivatives, i) for i in range(len(bumpy_derivatives))])
        d_pass = np.column_stack((noise_reduced[:-1,0], derivatives))
        bumpy_second_derivatives = [d(d_pass, i) for i in range(len(d_pass) - 1)]
        second_derivatives = np.array([smear(bumpy_second_derivatives, i) for i in range(len(bumpy_second_derivatives))])
        
        d_weights = 1 / (2*max(np.abs(derivatives[data_above[:-1]]))/3 + np.abs(derivatives[data_above[:-1]]))
        sd_weights = (1 + -second_derivatives / max(np.abs(second_derivatives)))/2

        averaging_weights = np.sqrt(data[data_above, 1]) * sd_weights[data_above[:-2]]

        #plt.plot(noise_reduced[:,0], noise_reduced[:,1], c = 'r', label = "noise reduced data")
        #plt.plot(noise_reduced[:-1,0], max(noise_reduced[:,1]) / max(np.abs(derivatives)) * np.array(derivatives), c = 'b', label = "scaled first derivatives")
        #plt.plot(noise_reduced[:-2,0], max(noise_reduced[:,1]) / max(np.abs(second_derivatives)) * second_derivatives, c = 'purple', label = "scaled second derivatives")
        #plt.plot(noise_reduced[data_above,0], max(noise_reduced[:,1]) * sd_weights[data_above[:-2]], c = 'lime', label = "scaled second derivative weight factor")
        
        #averaging_weights = data[data_above,1] * np.sqrt(d_weights)
        #plt.plot(noise_reduced[:,0], noise_reduced[:,1] * 10, c = 'r', label = "10x noise reduced data")
        #plt.plot(noise_reduced[:-1,0], bumpy_derivatives, c = 'cyan', label = "derivative of data")
        #plt.plot(noise_reduced[:-1,0], derivatives, c = 'b', label = "smeared derivative of data")
        #plt.plot(noise_reduced[:-2,0], second_derivatives / 100, c = 'purple', label = "1/100 x smeared second derivative")
        #plt.plot(noise_reduced[data_above,0], 1e6 * max(derivatives)/2 * d_weights, c = 'magenta', label = "scaled inverse derivative of data")
        #plt.plot(noise_reduced[data_above,0], averaging_weights * 1e4, c = 'lime', label = "1e4x assigned weights")
        #plt.axhline(y = 0, c = 'g', ls = '--')
        #plt.legend()
        #plt.show()

        
        
    
    
    
    
    w_max = np.dot(cut_wavelengths, averaging_weights) / np.sum(averaging_weights)
    w_err = np.sqrt(np.dot((cut_wavelengths - w_max)**2, averaging_weights) / np.sum(averaging_weights))

    if plot:
        plt.title("hwhm max model on scan for "+str(true_wavelength))
        plt.xlabel("Wavelength reading (A)")
        plt.ylabel("CPS")
        plt.plot(cut_wavelengths, averaging_weights * max(data[:,1]) / max(averaging_weights), label = "scaled averaging weights", c = 'lime')
        plt.errorbar(data[data_below,0], data[data_below,1], yerr = weights[data_below], label = "excluded data", marker = '.', c = 'b', ls='none')
        plt.errorbar(cut_wavelengths, data[data_above,1], yerr = cut_weights, label = "included data", marker = '.', c = 'orange', ls='none')
        plt.axvline(x = w_max, label = "wavelength estimate", c = 'r', linestyle = '--')
        plt.axvline(x = w_max - w_err, label = "error bounds", c = 'magenta', linestyle = '--')
        plt.axvline(x = w_max + w_err, c = 'magenta', linestyle = '--')
        plt.legend()
        plt.show()

    return w_max, w_err

def hwhm_max_for_losers(data, weights, plot = True, true_wavelength = None):

    data_half = (max(data[:,1]) + min(data[:,1])) / 2

    data_above = data[:,1] > data_half
    data_below = data[:,1] <= data_half
    cut_wavelengths = data[data_above,0]
    cut_weights = weights[data_above]

    averaging_weights = cut_weights * data[data_above,1]
    
    w_max = np.dot(cut_wavelengths, averaging_weights) / np.sum(averaging_weights)
    w_err = np.sqrt(np.dot((cut_wavelengths - w_max)**2, averaging_weights) / np.sum(averaging_weights))

    if plot:
        plt.title("hwhm max model on scan for "+str(true_wavelength))
        plt.xlabel("Wavelength reading (A)")
        plt.ylabel("CPS")
        plt.errorbar(data[data_below,0], data[data_below,1], yerr = weights[data_below], label = "excluded data", marker = '.', c = 'b', ls='none')
        plt.errorbar(cut_wavelengths, data[data_above,1], yerr = cut_weights, label = "included data", marker = '.', c = 'orange', ls='none')
        plt.axvline(x = w_max, label = "wavelength estimate", c = 'r', linestyle = '--')
        plt.axvline(x = w_max - w_err, label = "error bounds", c = 'magenta', linestyle = '--')
        plt.axvline(x = w_max + w_err, c = 'magenta', linestyle = '--')
        plt.legend()
        plt.show()

    return w_max, w_err

def do_the_thing():
    from data.data_loader import read_data
    #file = 'data/DeuteriumScans/SuperFineScans/beta2H'
    #from main_calibration_2 import check_3125_668 as d
    #file = d['fp']
    folder = 'data/final_data/sodium/'
    file = folder + "e" + "H"
    data = read_data(file)
    from noise_reduction import reduce_noise
    noise_reduced, weights = reduce_noise(data, damping_constant = 1/3)

    print (hwhm_max(data, weights, plot = True, true_wavelength = "sodium E H", threshold = 1/3))#, true_wavelength = d['true value']))
    print (hwhm_max(data, weights, noise_reduced, plot = True, true_wavelength = "sodium E H", threshold = 1/3))#, true_wavelength = d['true value']))
    #print (hwhm_max_for_losers(data, weights, plot = True))

if __name__ == "__main__":
    do_the_thing()


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
    
    
    background, signals = regions(data[:,1])
    print(len(signals))
    print (signals)

    standard_deviations = 2 # For more narrow peaks, may want to change this to 3? 

    maxes = [max(data[signal[0]:signal[1]+1, 1]) for signal in signals]
    widths = [(signal[1] - signal[0])/2/standard_deviations for signal in signals]

    return maxes, widths
    

#THESE THINGS ARE NOW IN data_processing
#this is called elsewhere so might break stuff
#please check where it breaks and change it to call these methods from data_processing instead of max_model :)

'''
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
'''



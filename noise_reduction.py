import numpy as np
from scipy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt 
from data.data_loader import read_data

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

    x, y = data[:,0], data[:,1]

    ysize = y.size
    xsize = x.size
    step = (max(x) - min(x) + 1)/xsize
    #print(ysize)
    #print(xsize)
    #print(step)

    frequencies, amplitudes = fftfreq(ysize, step), fft(y)

    #print (frequencies)

    #plot(x, y, "Data", "Wavelength (A)", "CPS")

    #print(amplitudes)
    #print(np.abs(amplitudes))

    #plot(frequencies, np.abs(amplitudes), "FFT of data", "frequency", "FFT")

    new_amplitudes, removed_amplitudes = bandpass(frequencies, amplitudes)

    #plot(x, ifft(removed_amplitudes), "Removed stuff", "Wavelength", "CPS", show = False)
    #plot(x, np.abs(ifft(removed_amplitudes)), "Removed stuff", "Wavelength", "CPS")

    new_data = ifft(new_amplitudes)
    removed = ifft(removed_amplitudes)

    #plot(x, y, "Data", "Wavelength (A)", "CPS", show = False)
    #plot(x, new_data, "IFFT data", "Wavelength (A)", "CPS", show = True)

    #raise NotImplementedError

    neighborhood = int(len(x) // 40)

    errors = [np.sqrt(np.sum(removed[max(i - neighborhood, 0):min(i + neighborhood + 1, len(removed))] ** 2) / (min(i + neighborhood + 1, len(removed)) - max(i - neighborhood, 0))) for i in range(len(removed))]
    print (errors)
    #plt.plot(x, errors)
    plt.errorbar(x, new_data, errors)
    plt.plot(x, y)
    plt.show()

    return new_data, 1/errors

def bandpass(frequencies, amplitudes):
    """
    Given frequencies and amplitudes, uses the cutoff method from max_model to retain the dominant frequencies and drop the higher frequencies, which mostly contribute to noise.
    Will consider adding an exponential damping in the future. 

    Arguments:
        * `frequencies` (numpy array): Frequencies of fourier transform of data.
        * `amplitudes` (numpy array): Amplitudes of fourier transform of data. Indices correspond with `frequencies`. 

    Returns:
        * As a 2-tuple:
            `new_amplitudes` (numpy array of same dimension as `amplitudes`): Array of new amplitudes. 
            `removed_amplitudes` (numpy array of same dimension as `amplitudes`): Array of removed amplitudes. These two should sum to the original. 
    """

    #consider: exponential damping of unused amplitudes. 

    from max_model import get_cutoff

    cutoff = get_cutoff(np.abs(amplitudes), 0)

    return np.where(np.abs(amplitudes) > cutoff, amplitudes, 0), np.where(np.abs(amplitudes) <= cutoff, amplitudes, 0)


def bandpass_old(frequencies, amplitudes):

    #use max_model to find width?

    from max_model import regions

    backgrounds, signals = regions(np.abs(amplitudes), 0, reduce = False)

    print(regions(np.abs(amplitudes), 0, reduce = False))

    # rms of cosine is 1/sqrt(2)

    cutoff = 0.1

    fmin = min(frequencies)
    span = max(frequencies) - min(frequencies)
    lower = fmin + span * (1-cutoff)/2
    upper = fmin + span * (1+cutoff)/2

    #indices = [frequency > lower and frequency < upper for frequency in frequencies]
    signal_indices = []
    for signal in signals:
        signal_indices.append(range(signal[0], signal[1]+1))
    indices = [i in signal_indices for i in range(len(amplitudes))]

    return frequencies, np.where(indices, amplitudes, 0)


def plot(x, y, title, xlabel, ylabel, show = True):

    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if show:
        plt.show()

if __name__ == "__main__":
    file = "data/DeuteriumScans/SuperFineScans/alpha4H"
    data = read_data(file)
    reduce_noise(data)
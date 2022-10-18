#Imports
import numpy as np
import csv
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm, chi2 # normal distribution, chi squared distribution
from numpy import exp, pi, sqrt
import lmfit
from max_model import regions
from noise_reduction import reduce_noise
from models import voigt_models

def execute_peak_fit(data, shift = 0, plot = False):

    new_data, weights = reduce_noise(data)
    
    i = primary_signal_region(data)
    
    if shift == 0:
        choice = 'voigt_with_shift'
        params = voigt_models[choice][1](data[i,:], shift)
    else:
        choice = 'two_voigt'
        params = voigt_models[choice][1](data[i,:])
    
    model = lmfit.Model(voigt_models[choice][0])

    result = fit(model, data[i,:], params, weights[i])

    amp, mu, alpha, gamma = result.params['amp'].value, result.params['mu'].value, result.params['alpha'].value, result.params['gamma'].value
    mu_err = result.params['mu'].stderr
    #different handling for double peaks with similar amplitude?

    if 'a' in result.params.keys():
        a = result.params['a'].value
    else:
        a = 0

    if 'mu2' in result.params.keys():
        amp2, mu2, alpha2, gamma2 = result.params['amp2'].value, result.params['mu2'].value, result.params['alpha2'].value, result.params['gamma2'].value
                
    if plot:
        plt.scatter(data[i,0], data[i,1], marker = '.', label = "data")
        plt.plot(data[i,0], voigt_models['voigt_with_shift'][0](data[i,0], amp, mu, alpha, gamma, a), label = "primary peak", color = "orange")
        if choice == 'two_voigt':
            amp2, mu2, alpha2, gamma2 = result.params['amp2'].value, result.params['mu2'].value, result.params['alpha2'].value, result.params['gamma2'].value
            plt.plot(data[i,0], voigt_models[choice][0](data[i,0], amp, mu, alpha, gamma, amp2, mu2, alpha2, gamma2, a), label = "full fit", color = "lime")
        plt.axvline(x = mu, label = "wavelength estimate", linestyle = '--', color = 'r')
        plt.legend()
        plt.show()

    return mu, mu_err


def primary_signal_region (data):

    backgrounds, signals = regions(data[:,1])

    maxval = 0
    max_signal = (0,1)
    for signal in signals:
        smax = max(data[:,1][signal[0]:signal[1]+1])
        if smax > maxval:
            maxval = smax
            max_signal = signal
    
    i = np.array(range(len(data[:,1])))
    return np.logical_and(i >= max_signal[0], i <= max_signal[1])


def fit (model, data, params=None, weights =None):
    """
    implements lmfit to fit given data to a specified model

    Arguments: 
        * `model` (Model): lmfit Model to fit to
        * `data` (nx2 numpy array): numpy array (matrix) with 2 columns:
            - First column is x-axis (input to Model)
            - Second column is y-axis (output of Model)
        * `params` (Parameters): initial parameters values for the model, defaults to None
        * `weights` (n dimensional array): weights for each data point (x,y) in data, defaults to None

    Returns: 
        * `result` (ModelResult) result of fit (returned from model.fit)
    """
    x_axis = data[:, 0]
    y_axis = data[:, 1]
    
    if(params is None and weights is None):
        result = model.fit(y_axis, x=x_axis)
    elif (params is None):
        result = model.fit(y_axis, x=x_axis, weights=weights)
    elif (weights is None):
        result = model.fit(y_axis, params=params, x=x_axis)
    else:
        result = model.fit(y_axis, params=params, x=x_axis, weights=weights)
    return result
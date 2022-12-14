#Imports
from cmath import e
import numpy as np
import csv
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm, chi2 # normal distribution, chi squared distribution
from numpy import exp, pi, sqrt
import lmfit
from data.data_loader import read_data
from max_model import regions
from models_2 import voigt_models
from data_processing import process_data


def execute_peak_fit(data, shift = 0, plot = False):

    from noise_reduction import reduce_noise

    new_data, weights = reduce_noise(data)
    
    i = primary_signal_region(data)
    
    if shift == 0:
        choice = 'voigt_with_shift'
        params = voigt_models[choice][1](data[i,:])
    else:
        choice = 'two_voigt'
        params = voigt_models[choice][1](data[i,:], shift)
    
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
                
    """
    if plot:
        plt.scatter(data[i,0], data[i,1], marker = '.', label = "data")
        plt.plot(data[i,0], voigt_models['voigt_with_shift'][0](data[i,0], amp, mu, alpha, gamma, a), label = "primary peak", color = "orange")
        if choice == 'two_voigt':
            amp2, mu2, alpha2, gamma2 = result.params['amp2'].value, result.params['mu2'].value, result.params['alpha2'].value, result.params['gamma2'].value
            plt.plot(data[i,0], voigt_models[choice][0](data[i,0], amp, mu, alpha, gamma, amp2, mu2, alpha2, gamma2, a), label = "full fit", color = "lime")
        plt.axvline(x = mu, label = "wavelength estimate", linestyle = '--', color = 'r')
        plt.legend()
        plt.show()
    """
    if plot:
        plt.scatter(data[:,0], data[:,1], marker = '.', label = "data")
        plt.plot(data[:,0], voigt_models['voigt_with_shift'][0](data[:,0], amp, mu, alpha, gamma, a), label = "primary peak", color = "orange")
        if choice == 'two_voigt':
            amp2, mu2, alpha2, gamma2 = result.params['amp2'].value, result.params['mu2'].value, result.params['alpha2'].value, result.params['gamma2'].value
            plt.plot(data[:,0], voigt_models[choice][0](data[:,0], amp, mu, alpha, gamma, amp2, mu2, alpha2, gamma2, a), label = "full fit", color = "lime")
        plt.axvline(x = mu, label = "wavelength estimate", linestyle = '--', color = 'r')
        #plt.axvline(x = true, label = "true wavelength", linestyle = '--', color = 'magenta')
        plt.legend()
        plt.show()

    return mu, mu_err

# REDUNDANT WITH Execute Peak fit
# IMPLEMENTS IN DATA PROCESSING NOW
# Also returns the entire result cuz its useful
#

def fit_to_voigt(processed_data, weights, shift = 0, damping_constant = 1/10, plot = False, title = None):
    
    if shift == 0:
        choice = 'voigt_with_shift'
        params = voigt_models[choice][1](processed_data)
    else:
        choice = 'two_voigt'
        params = voigt_models[choice][1](processed_data, shift, plot = plot)
    
    model = lmfit.Model(voigt_models[choice][0])

    result = fit(model, processed_data, params, weights)

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
    #if True:
        #plt.title('Voigt Fit has reduced chi2 = ' + str(round(result.redchi, 2)))
        if title is not None:
            plt.title(title)
        plt.errorbar(processed_data[:,0], processed_data[:,1], yerr = 1/weights, marker = '.', label = "data", ls = 'none')
        plt.plot(processed_data[:,0], voigt_models['voigt_with_shift'][0](processed_data[:,0], amp, mu, alpha, gamma, a), label = "primary peak", color = "orange")
        if choice == 'two_voigt':
            amp2, mu2, alpha2, gamma2 = result.params['amp2'].value, result.params['mu2'].value, result.params['alpha2'].value, result.params['gamma2'].value
            plt.plot(processed_data[:,0], voigt_models[choice][0](processed_data[:,0], amp, mu, alpha, gamma, amp2, mu2, alpha2, gamma2, a), label = "full fit", color = "lime")
            plt.plot(processed_data[:,0], voigt_models['voigt_with_shift'][0](processed_data[:,0], amp2, mu2, alpha2, gamma2, a), label = "secondary peak", color = "green")
        plt.axvline(x = mu, label = "wavelength estimate", linestyle = '--', color = 'r')
        plt.ylabel('Counts per second')
        plt.xlabel('Monochromator Step')
        plt.legend()
        plt.show()
    print ("chi square: "+str(result.chisqr))
    print ("chi square: "+str(result.redchi))

    print(lmfit.fit_report(result))

    return result.params

# essentially, uses max model on noise-reduced data
# used for wide, small, noisy peaks, which are too hard for our fitting to detect
# add dictionary to maxmodel_peaks in main_calibration to use this model
def naive_fit(data, damping_constant):

    from noise_reduction import reduce_noise
    from max_model import get_max
    new_data, weights = reduce_noise(data, damping_constant = damping_constant, plot = True)
    result = get_max(new_data)
    print (result)
    if result[0] == []:
        return data[np.argmax(data[:,1]),0], (max(data[:,0]) - min(data[:,0])) / 4
    return result[0][0], result[1][0]

#
# IMPLEMENTED IN DATA PROCESSING NOW
# Can be deleted now
#
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
        result = model.fit(y_axis, x=x_axis, scale_covar = False)
    elif (params is None):
        result = model.fit(y_axis, x=x_axis, weights=weights, scale_covar = False)
    elif (weights is None):
        result = model.fit(y_axis, params=params, x=x_axis, scale_covar = False)
    else:
        result = model.fit(y_axis, params=params, x=x_axis, weights=weights, scale_covar = False)
    return result






def fit_to_lorentzian(processed_data, weights, plot=True):
    """
    Given data from a scan, fits a lorentzian to it. Outdated--better to use voigt. 

    Arguments:
        * `processed_data` (nx2 numpy array): 2 columns. First is wavelength in angstroms, second is measured counts per second. 
        * `weights`

    Returns:
        * `result` (ModelResult): returned by model.fit
    """
    #type of data is numpy array w/ 2 columns
    def Lorentzian(x, amp, cen, scale, shift):
        return amp * (1/(pi*scale))*(1/(1+(((x-cen)/scale)**2))) + shift
    
    
    model = lmfit.Model(Lorentzian)
    
    x_axis = processed_data[:, 0]
    y_axis = processed_data[:, 1]
    
    start_amp = max(y_axis)
    start_cen = x_axis[np.argmax(y_axis)]
    start_shift = min(y_axis)
    
    result = model.fit(y_axis, x=x_axis, weights=weights, amp=start_amp, cen=start_cen, scale=0.5, shift=start_shift)
    
    amp = result.best_values['amp']
    cen = result.best_values['cen']
    scale = result.best_values['scale']
    shift = result.best_values['shift']


    if plot:            
    #if True:
        plt.title('Lorentzian Fit has reduced chi2 = ' + str(round(result.redchi, 2)))
        plt.errorbar(processed_data[:,0], processed_data[:,1], yerr = 1/weights, marker = '.', label = "data")
        plt.plot(processed_data[:,0], Lorentzian(processed_data[:,0], amp, cen, scale, shift), label = "primary peak", color = "orange")
        plt.axvline(x = cen, label = "wavelength estimate", linestyle = '--', color = 'r')
        plt.ylabel('counts per second')
        plt.xlabel('uncalibrated angstroms')
        plt.legend()
        plt.show()
    print ("chi square: "+str(result.chisqr))
    print ("chi square: "+str(result.redchi))

    print(lmfit.fit_report(result))
    
    return result


def fit_to_Gaussian(processed_data, weights, plot=True):
    """
    Given data from a scan, fits a lorentzian to it. Outdated--better to use voigt. 

    Arguments:
        * `processed_data` (nx2 numpy array): 2 columns. First is wavelength in angstroms, second is measured counts per second. 
        * `weights`

    Returns:
        * `result` (ModelResult): returned by model.fit
    """
    #type of data is numpy array w/ 2 columns
    def Gaussian(x, amp, cen, sigma, shift):
        return amp * exp(-(1/2) * ((x-cen)/sigma)**2) + shift
    
    
    model = lmfit.Model(Gaussian)
    
    x_axis = processed_data[:, 0]
    y_axis = processed_data[:, 1]
    
    start_amp = max(y_axis)
    start_cen = x_axis[np.argmax(y_axis)]
    start_shift = min(y_axis)
    
    result = model.fit(y_axis, x=x_axis, weights=weights, amp=start_amp, cen=start_cen, sigma=0.5, shift=start_shift)
    
    amp = result.best_values['amp']
    cen = result.best_values['cen']
    sigma = result.best_values['sigma']
    shift = result.best_values['shift']


    if plot:            
    #if True:
        plt.title('Gaussian Fit has reduced chi2 = ' + str(round(result.redchi, 2)))
        plt.errorbar(processed_data[:,0], processed_data[:,1], yerr = 1/weights, marker = '.', label = "data")
        plt.plot(processed_data[:,0], Gaussian(processed_data[:,0], amp, cen, sigma, shift), label = "primary peak", color = "orange")
        plt.axvline(x = cen, label = "wavelength estimate", linestyle = '--', color = 'r')
        plt.ylabel('counts per second')
        plt.xlabel('uncalibrated angstroms')
        plt.legend()
        plt.show()
    print ("chi square: "+str(result.chisqr))
    print ("chi square: "+str(result.redchi))

    print(lmfit.fit_report(result))
    
    return result
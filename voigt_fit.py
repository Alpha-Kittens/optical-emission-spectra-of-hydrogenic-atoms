#Imports
import numpy as np
import csv
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm, chi2 # normal distribution, chi squared distribution
from numpy import exp, pi, sqrt
import lmfit
from data import data_loader

def fit_to_voigt(data):
    #type of data is numpy array w/ 2 columns
    '''
    def Lorentzian(x, amp, cen, scale):
        return amp * (1/(pi*scale))*(1/(1+(((x-cen)/scale)**2)))
    '''
    
    model = lmfit.models.VoigtModel()
    
    x_axis = data[:, 0]
    y_axis = data[:, 1]
    
    start_amp = max(y_axis)
    start_cen = x_axis[np.argmax(y_axis)]
    
    result = model.fit(y_axis, x=x_axis, amplitude=start_amp, center=start_cen, sigma = 0.5)
    
    return result


data = data_loader.read_data('data/day_1/calibration_scan_9_4019_4023')
result = fit_to_voigt(data)
#show_results(result)

lmfit.report_fit(result)
plt.figure()
plt.figure(figsize=(100, 100))
result.plot()
plt.show()
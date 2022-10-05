#Imports
import numpy as np
import csv
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm, chi2 # normal distribution, chi squared distribution
from numpy import exp, pi, sqrt
import lmfit

import data.data_loader as loader

loader.read_data()

def fit_to_lorentzian(data):
    """
    Given 
    """
    #type of data is numpy array w/ 2 columns
    def Lorentzian(x, amp, cen, scale):
        return amp * (1/(pi*scale))*(1/(1+(((x-cen)/scale)**2)))
    
    
    model = lmfit.Model(Lorentzian)
    
    x_axis = data[:, 0]
    y_axis = data[:, 1]
    
    start_amp = max(y_axis)
    start_cen = x_axis[np.argmax(y_axis)]
    
    result = model.fit(y_axis, x=x_axis, amp=start_amp, cen=start_cen, scale=0.5)
    
    
    return result


def show_results(result):
    lmfit.report_fit(result)
    plt.figure()
    plt.figure(figsize=(100, 100))
    result.plot()
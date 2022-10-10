#Imports
import numpy as np
import csv
import math
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm, chi2 # normal distribution, chi squared distribution
from numpy import exp, pi, sqrt
import lmfit

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
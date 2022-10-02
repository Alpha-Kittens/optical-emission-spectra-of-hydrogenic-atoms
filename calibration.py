import import_ipynb
import lorentzian_fit

from data.data_loader import read_data
from lmfit import Model, Parameter, fit

reference = {
    'mercury' : {
        'reference' : [],
        'check' : [],
    },

}

def get_band(file):
    """
    Given file, runs Athira's lorenzian_fit in order to extract a mean and "scale" (~width) 
    """
    
    raise NotImplementedError
    return mean, scale


def get_calibration_data(files):
    """
    Given list of files (consdier: element name, and it just knows the files), extracts the band details from each of them. 
    Returs centers and weights
    """

    raise NotImplementedError

def fit_calibration(files, element):

    raise NotImplementedError

def check_calibration(c_function, references):
    """
    Given a list of smaller-line calibration references, checks validity of fit. 
    """

    raise NotImplementedError


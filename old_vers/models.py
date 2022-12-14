import numpy as np
import lmfit
#import random
from max_model import regions

def linear(x, m, b):
    """
    linear model y=mx+b

    Arguments: 
        * `x` : parameter value to evaluate the function at
        * `m` : slope of the line
        * `b` : y-intercept of the line

    Returns: 
        * `y` result of y=mx+b
    """
    return (m*x)+b

def quadratic(x, a, b, c):
    """
    quadratic model y=ax^2 + bx + c

    Arguments: 
        * `x` : parameter value to evaluate the function at
        * `a` : leading coefficient
        * `b` : 2nd term coefficient
        * `c` : y-interecept

    Returns: 
        * `y` result of y=ax^2 + bx + c
    """
    return (a*(x**2)) + (b*x) + c

def eval_voigt_with_shift(x, amplitude, center, sigma, a):

    voigt_model = lmfit.models.VoigtModel()
    params = voigt_model.make_params(amplitude=amplitude, center=center, sigma = sigma)
    return voigt_model.eval(params, x = x) + a

def voigt_with_shift():

    return lmfit.Model(eval_voigt_with_shift)

def shifted_voigt_params(model, data):

    """
    returns initial parameters to input for a voigt_model

    Arguments: 
        * `model` (VoigtModel): the model to make intial values for
        * `data`  (nx2 numpy array): numpy array (matrix) with 2 columns:
            - First column is x-axis (input to Model)
            - Second column is y-axis (output of Model)

    Returns: 
        * `params` (Parameters) with intial values based on data
    """
    x_axis = data[:, 0]
    y_axis = data[:, 1]

    start_amp = max(y_axis)
    start_cen = x_axis[np.argmax(y_axis)]
    start_a = min(y_axis)

    params = lmfit.Parameters()
    params.add('amplitude', value = start_amp)
    params.add('center', value = start_cen)
    params.add('sigma', value = 0.5)
    params.add('a', value = start_a)


    return params  

def eval_two_voigts(x, amplitude1, center1, sigma1, amplitude2, center2, sigma2, a):
    voigt_model = lmfit.models.VoigtModel()
    params1 = voigt_model.make_params(amplitude=amplitude1, center=center1, sigma = sigma1)
    params2 = voigt_model.make_params(amplitude=amplitude2, center=center2, sigma = sigma2)
    return voigt_model.eval(params1, x = x) + voigt_model.eval(params2, x = x) + a

def two_voigts():

    return lmfit.Model(eval_two_voigts)

def two_voigt_params(model, data):

    """
    returns initial parameters to input for a voigt_model

    Arguments: 
        * `model` (VoigtModel): the model to make intial values for
        * `data`  (nx2 numpy array): numpy array (matrix) with 2 columns:
            - First column is x-axis (input to Model)
            - Second column is y-axis (output of Model)

    Returns: 
        * `params` (Parameters) with intial values based on data
    """
    x_axis = data[:, 0]
    y_axis = data[:, 1]

    start_amp = max(y_axis)
    signals = regions(y_axis)[1]
    amplitudes = [max(y_axis[signal[0]:signal[1]+1]) for signal in signals]

    i1 = np.argmax(amplitudes)
    max1 = amplitudes[i1]
    i2 = np.argmin(np.where(amplitudes == max1, 0, amplitudes))
    max2 = amplitudes[i2]
    maxes = [max1, max2]
    indices = [i1, i2]

    start_cen1, start_cen2 = (x_axis[np.argmax(y_axis[signals[indices[i]][0]:signals[indices[i]][1]])] for i in range(2))
    start_sigma1, start_sigma2 = (1/2 * (x_axis[signals[indices[i]][1]] - x_axis[signals[indices[i]][0]]) for i in range(2))
    start_amplitude1, start_amplitude2 = maxes
    start_a = min(y_axis)

    # consider: using regions to establish two independetn peaks
    params = lmfit.Parameters()
    params.add('amplitude1', value = start_amplitude1)
    params.add('center1', value = start_cen1)
    params.add('sigma1', value = start_sigma1)
    params.add('amplitude2', value = start_amplitude2)
    params.add('center2', value = start_cen2)
    params.add('sigma2', value = start_sigma2)
    params.add('a', value = start_a)

    return params

def quadratic_err(x, x_err, a, a_err, b, b_err, c, c_err):
    """
    Error for result of quadratic model

    Arguments:
        * `x` : parameter value to evaluate function
        * `x_err` : error associated with `x`
        * `a` : leading coefficient
        * `a_err` : error associated with `a`
        * `b` : linear coefficient
        * `b_err` : error associated with `b`
        * `c` : y-intercept
        * `c_err` : error associated with `c`
    """

    x_rel = x_err / x
    a_rel = a_err / a
    b_rel = b_err / b
    c_rel = c_err / c

    return ((a * x**2)**2 * (a_rel**2 + 2 * x_rel**2) + (b * x)**2 * (b_rel ** 2 + x_rel**2) + c_err ** 2) ** (1/2)


def voigt_params(model, data):
    """
    returns initial parameters to input for a voigt_model

    Arguments: 
        * `model` (VoigtModel): the model to make intial values for
        * `data`  (nx2 numpy array): numpy array (matrix) with 2 columns:
            - First column is x-axis (input to Model)
            - Second column is y-axis (output of Model)

    Returns: 
        * `params` (Parameters) with intial values based on data
    """
    x_axis = data[:, 0]
    y_axis = data[:, 1]

    start_amp = max(y_axis)
    start_cen = x_axis[np.argmax(y_axis)]

    params = model.make_params(amplitude=start_amp, center=start_cen, sigma = 0.5)

    return params

def extract_voigt(result):
    """
    returns amplitude, center, sigma
    """

    return result.params['amplitude'], result.params['center'], result.params['sigma']
import numpy as np
import lmfit

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

    params = lmfit.Parameters()
    params.add('amplitude', value = start_amp)
    params.add('center', value = start_cen)
    params.add('sigma', value = 0.5)
    params.add('a', value = 400)


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
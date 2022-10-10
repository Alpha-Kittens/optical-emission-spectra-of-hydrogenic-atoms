import numpy as np

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

def voigt_params(model, data):
    """
    returns initial parameters for a voigt_model

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

    
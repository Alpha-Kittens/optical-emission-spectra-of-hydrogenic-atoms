from calibration_models import *
from numpy import sqrt, abs
params = {'a': 4332.627208961621, 'b': 1.0075092615313312, 'c': 9.50624119947782e-07, 'd': 1.397895884872255e-10}
errors = {'a': 0.07866639000040926, 'b': 0.00018992527907119605, 'c': 5.73583395782815e-08, 'd': 1.1676348069758903e-10}
calibrate = lambda x : model_poly(4305.053269035736)(x, **params)
calibrate_error = lambda x, x_err : model_poly_err(4305.053269035736, params, errors)(x, x_err)
calibrate_splitting = lambda x1, x2 : calibrate(x2) - calibrate(x1)
calibrate_splitting_error = lambda x1, x2, x1_err, x2_err : (sqrt(calibrate_error(x1, x1_err)[0]**2 + calibrate_error(x2, x2_err)[0]**2), abs(calibrate_error(x1, x1_err)[1] - calibrate_error(x2, x2_err)[1]))

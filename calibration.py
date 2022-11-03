from calibration_models import *
from numpy import sqrt, abs
params = {'a': 4466.22481894433, 'b': 1.007913925413525, 'c': 1.0527937957154203e-06}
errors = {'a': 0.036516332418422656, 'b': 2.4322042323210507e-05, 'c': 3.2745852719795714e-08}
calibrate = lambda x : model_poly(4437.625305688662)(x, **params)
calibrate_error = lambda x, x_err : model_poly_err(4437.625305688662, params, errors)(x, x_err)
calibrate_splitting = lambda x1, x2 : calibrate(x2) - calibrate(x1)
calibrate_splitting_error = lambda x1, x2, x1_err, x2_err : (sqrt(calibrate_error(x1, x1_err)[0]**2 + calibrate_error(x2, x2_err)[0]**2), abs(calibrate_error(x1, x1_err)[1] - calibrate_error(x2, x2_err)[1]))

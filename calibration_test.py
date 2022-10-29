from calibration_models import *
from numpy import sqrt, abs
params = {'a': 4332.627211720854, 'b': 1.007509272547461, 'c': 9.50626368347023e-07, 'd': 1.3978215205237124e-10}
errors = {'a': 0.07866886122900443, 'b': 0.0001899350068267682, 'c': 5.735578178301609e-08, 'd': 1.1676835780530544e-10}
calibrate = lambda x : model_poly(4305.053269035736)(x, **params)
calibrate_error = lambda x, x_err : model_poly_err(4305.053269035736, params, errors)(x, x_err)
calibrate_splitting = lambda x1, x2 : calibrate(x2) - calibrate(x1)
calibrate_splitting_error = lambda x1, x2, x1_err, x2_err : (sqrt(calibrate_error(x1, x1_err)[0]**2 + calibrate_error(x2, x2_err)[0]**2), abs(calibrate_error(x1, x1_err)[1] - calibrate_error(x2, x2_err)[1]))

#ur momfrom calibration_models import *
from numpy import sqrt, abs
params = {'a': 4466.186982141817, 'b': 1.0080232194275698, 'c': 1.2989189649007526e-06, 'd': -2.2697309873004072e-10, 'e': -1.3643439229457626e-13, 'f': 9.597878626182218e-17}
errors = {'a': 0.05908288140193235, 'b': 0.00034376991417852846, 'c': 9.520697263047814e-08, 'd': 6.238295507735593e-10, 'e': 5.0948724842358545e-14, 'f': 2.466377804657649e-16}
calibrate = lambda x : model_poly(4437.625305688662)(x, **params)
calibrate_error = lambda x, x_err : model_poly_err(4437.625305688662, params, errors)(x, x_err)
calibrate_splitting = lambda x1, x2 : calibrate(x2) - calibrate(x1)
calibrate_splitting_error = lambda x1, x2, x1_err, x2_err : (sqrt(calibrate_error(x1, x1_err)[0]**2 + calibrate_error(x2, x2_err)[0]**2), abs(calibrate_error(x1, x1_err)[1] - calibrate_error(x2, x2_err)[1]))

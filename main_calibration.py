# remaining TODO:
# -handle error extraction of parameters from quadratic fit
# -more nuanced handling of weights post-fit using dy/dx method
# chi-square of calibration test


from data.data_loader import read_data
from fitters import fit
from models import *
from noise_reduction import reduce_noise
import matplotlib.pyplot as plt
import lmfit
import numpy as np

# File Path to main peak data files
folder = "data/mercury/SuperFineScans"

# Define Main Peak Expected Lines
wavelengths = {
    'reference': [],
    'check': [],
}

wavelengths_errors = {
    'reference': [],
    'check' : [],
}

files = {
    'refrence': [],
    'check' : [],
}

for key in files:
    for string in files[key]:
        string = folder + string

# Read data
# Noise reduction and find Uncertainties

data = {
    'reference': [(reduce_noise(read_data(files['reference'][i]))) for i in range(len(files['reference']))],
    'check': [(reduce_noise(read_data(files['check'][i]))) for i in range(len(files['check']))],
}

results_reference = []
error_reference = []

# Voigt Model & Fit
for key in data:
    for entry in data[key]:

        new_data, weights = entry
        
        model = lmfit.models.VoigtModel()
        params = voigt_params(model, new_data)

        result = fit(model, new_data, params, weights)

        cps, wavelength, error = extract_voigt(result)

        results_reference.append(wavelength)
        error_reference.append(error)

# Maximum Counts
#??

# Find Residuals
# no need

# Quadratic Model and Fit

calibration_model = quadratic
calibration_error = quadratic_err

inputs = np.zeros((len(results_reference), 2))
inputs[:,0] = results_reference
inputs[:,1] = wavelengths['reference']

result = fit(model, inputs, weights = error_reference)

a, b, c = result.params['a'], result.params['b'], result.params['c']

calibration = lambda x : calibration_model(x, a, b, c)
a_err = -1 #result.params['a'].stderr
b_err = -1 #result.params['b'].stderr
c_err = -1 #result.params['c'].stderr

calibration_error = lambda x, x_err : quadratic_err(x, x_err, a, a_err, b, b_err, c, c_err)

# File path to small peak data files
# Read data
# Voigt Model & Fit
#completed earlier

# Test quadratic model
# better way of doing this with arrays probably exists

data_check = np.array([data['check'][i][0] for i in range(len(data['check']))])
data_check_err = np.array([data['check'][i][1] for i in range(len(data['check']))])
calibration_test = calibration(data_check)
calibration_test_err = calibration_error(data_check, data_check_err)

residuals = calibration_test - np.array(wavelengths['check'])

plt.plot(data_check, calibration_test)
plt.plot(data_check, wavelengths['check'])
plt.show()


# remaining TODO:
# -handle error extraction of parameters from quadratic fit
# -more nuanced handling of weights post-fit using dy/dx method
# chi-square of calibration test


from data.data_loader import read_data
from fitters import fit
from models_2 import *
from noise_reduction import reduce_noise
import matplotlib.pyplot as plt
import lmfit
import numpy as np

# File Path to main peak data files
folder = "data/10.14.22/mercury/"

# Define Main Peak Expected Lines
wavelengths = {
    # https://physics.nist.gov/PhysRefData/Handbook/Tables/mercurytable2.htm
    # possibly include parameter for hyperfine splitting
    'reference': [6149.475, 5790.663, 5769.598, 5460.735, 4046.563, 3650.153],
    'check': [],
}

# errors not provided for these measurements. 
"""
wavelengths_errors = {
    'reference': [],
    'check' : [],
}
"""

files = {
    'reference': ["6149.50-coarse", "5790.66-superfine", "5769.6-superfine", "5460.74-superfine", "4046.56-superfine_4", "3650.15-superfine"],
    'check' : [],
}

for key in files:
    for i in range(len(files[key])):
        files[key][i] = folder + files[key][i]

files['reference'][4] = "data/mercury/calibration_scan_9_4019_4023"

# Read data
# Noise reduction and find Uncertainties

# noise-reduced version: 
"""
data = {
    'reference': [(reduce_noise(read_data(files['reference'][i]))) for i in range(len(files['reference']))],
    'check': [(reduce_noise(read_data(files['check'][i]))) for i in range(len(files['check']))],
}
"""

# non noise-reduced version:
data = {
    'reference' : [read_data(file) for file in files['reference']],
    'check' : [read_data(file) for file in files['check']]
}

superfine_stepsize = 0.0025

# plot, user request splitting estimate

measured_reference = []
error_reference = []

# Voigt Model & Fit
for key in data:
    for entry in data[key]:

        new_data, weights = reduce_noise(entry)

        plt.scatter(entry[:,0], entry[:,1], marker = '.', label = "data")
        plt.plot(entry[:,0], new_data[:,1], label = "noise-reduced")
        plt.legend()
        plt.show()

        shift = float(input("Estimate the splitting in that plot (Angstroms): "))

        if shift != 0:
            choice = 'two_voigt'
        else:
            choice = 'voigt_with_shift'
        
        model = lmfit.Model(voigt_models[choice][0])
        if choice == 'two_voigt':
            params = voigt_models[choice][1](entry, shift)
        else:
            params = voigt_models[choice][1](entry)

        result = fit(model, entry, params, weights)

        extract = voigt_models[choice][2](result)
        
        amp, mu, alpha, gamma = extract[0][0:4]

        if len(extract[0]) > 4:
            a = extract[0][4]
        else:
            a = 0

        plt.scatter(entry[:,0], entry[:,1], marker = '.', label = "data")
        plt.plot(entry[:,0], voigt_models['voigt_with_shift'][0](entry[:,0], amp, mu, alpha, gamma, a), label = "primary peak")
        if choice == 'two_voigt':
            amp2, mu2, alpha2, gamma2 = extract[0][5:9]
            plt.plot(entry[:,0], voigt_models[choice][0](entry[:,0], amp, mu, alpha, gamma, amp2, mu2, alpha2, gamma2, a), label = "full fit")
        plt.axvline(x = mu, label = "wavelength estimate")
        plt.legend()
        plt.show()

        wavelength = mu
        error = extract[1][1]

        measured_reference.append(wavelength)
        error_reference.append(error)

# Maximum Counts
#??

# Find Residuals
# no need

# Quadratic Model and Fit

calibration_model = lmfit.Model(quadratic)
calibration_error = quadratic_err
"""

inputs = np.zeros((len(results_reference), 2))
inputs[:,0] = measured_reference
inputs[:,1] = wavelengths['reference']

result = fit(model, inputs, weights = error_reference)

a, b, c = result.params['a'], result.params['b'], result.params['c']

calibration = lambda x : quadratic(x, a, b, c)
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
if __name__ == "__main__":
    
    data_check = np.array([data['check'][i][0] for i in range(len(data['check']))])
    data_check_err = np.array([data['check'][i][1] for i in range(len(data['check']))])
    calibration_test = calibration(data_check)
    calibration_test_err = calibration_error(data_check, data_check_err)

    residuals = calibration_test - np.array(wavelengths['check'])


    plt.plot(data_check, calibration_test)
    plt.plot(data_check, wavelengths['check'])
    plt.show()
"""

#score the model/quantify performance


# def apply_model(wavelength_measured, wavelength_error):
""" 
    Applies the found calibration model. Can be imported by other functions. 

    Arguments: 
        * `wavelength_measured` (float): wavelength measured as reported by the spectrometer
        * `error` (float): error measurement of measured wavelength

    Returns:
        * As a 2-tuple:
            - `wavelength` (float): corrected wavelength
            - `wavelength_error` (float): error on corrected wavelength
"""
    
#    return calibration(wavelength_measured), calibration_error(wavelength_measured, wavelength_error)
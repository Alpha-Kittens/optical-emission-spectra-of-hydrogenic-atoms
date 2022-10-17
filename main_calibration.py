from re import M
from data.data_loader import read_data
from fitters import fit
from models_2 import *
from noise_reduction import reduce_noise
import matplotlib.pyplot as plt
import lmfit
import numpy as np

# File Path to main peak data files
folder = "data/10.14.22/mercury/"


# Small lines
small_lines = [4358.328, 6149.475, 3131.55, 3131.84, 5677.105, 4347.494]

# Define Main Peak Expected Lines
wavelengths = {
    # https://physics.nist.gov/PhysRefData/Handbook/Tables/mercurytable2.htm
    # possibly include parameter for hyperfine splitting
    'reference': [5790.663, 5769.598, 5460.735, 4046.563, 3650.153],
    'check': [6562, 4861, 4340.47, 4101.74], #hydrogen
}



# errors not provided for these measurements. 
"""
wavelengths_errors = {
    'reference': [],
    'check' : [],
}
"""

files = {
    'reference': ["5790.66-superfine", "5769.6-superfine", "5460.74-superfine", "4046.56-superfine_4", "3650.15-superfine"],
    'check' : ["alpha-superfine-1", "beta-superfine-1", "gamma-superfine-1", "delta-superfine-1"],
}

for key in ['reference']:
    for i in range(len(files[key])):
        files[key][i] = folder + files[key][i]


for key in ['check']:
    for i in range(len(files[key])):
        files[key][i] = "data/10.14.22/hydrogen/" + files[key][i]

#files['reference'][4] = "data/mercury/calibration_scan_9_4019_4023"

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
measured_check = []
error_check = []


plot = True

# Voigt Model & Fit
for key in data:
    for entry in data[key]:

        new_data, weights = reduce_noise(entry)

        if plot:
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

        amp, mu, alpha, gamma = result.params['amp'].value, result.params['mu'].value, result.params['alpha'].value, result.params['gamma'].value
        mu_err = result.params['mu'].stderr
        #different handling for double peaks with similar amplitude?

        if 'a' in result.params.keys():
            a = result.params['a'].value
        else:
            a = 0
            

        if plot:
            plt.scatter(entry[:,0], entry[:,1], marker = '.', label = "data")
            plt.plot(entry[:,0], voigt_models['voigt_with_shift'][0](entry[:,0], amp, mu, alpha, gamma, a), label = "primary peak", color = "orange")
            if choice == 'two_voigt':
                amp2, mu2, alpha2, gamma2 = result.params['amp2'].value, result.params['mu2'].value, result.params['alpha2'].value, result.params['gamma2'].value
                plt.plot(entry[:,0], voigt_models[choice][0](entry[:,0], amp, mu, alpha, gamma, amp2, mu2, alpha2, gamma2, a), label = "full fit", color = "lime")
            plt.axvline(x = mu, label = "wavelength estimate", linestyle = '--', color = 'r')
            plt.legend()
            plt.show()

        print ("wavelength: "+str(mu))
        print ("true: " + str(wavelengths['reference']))
        #if choice == 'two_voigt':
        #    print(mu2)
        #    print(amp2)
        wavelength = mu
        error = mu_err
        if key == 'reference':
            measured_reference.append(wavelength)
            error_reference.append(error)
        if key == 'check':
            measured_check.append(wavelength)
            error_check.append(error)

# Maximum Counts
#??

# Find Residuals
# no need

# Quadratic Model and Fit

plt.errorbar(measured_reference, wavelengths['reference'], yerr = error_reference, label = "calibration data", color = 'b', fmt = '.')
plt.xlabel("measured wavelength")
plt.ylabel("actual wavelength")
plt.legend()
plt.show()

calibration_model = lmfit.Model(quadratic)
params = lmfit.Parameters()
params.add('a', value = 0)
params.add('b', value = 1)
params.add('c', value = 0)
calibration_error = quadratic_err


inputs = np.zeros((len(measured_reference), 2))
inputs[:,0] = measured_reference
inputs[:,1] = wavelengths['reference']

result = fit(calibration_model, inputs, params, weights = error_reference) # a slight approximation of the true weights but whatever (the slope is roughly 1).

a, b, c = result.params['a'].value, result.params['b'].value, result.params['c'].value

calibration = lambda x : quadratic(x, a, b, c)
a_err = result.params['a'].stderr
b_err = result.params['b'].stderr
c_err = result.params['c'].stderr

chisqr = result.chisqr
redchi = result.redchi

calibration_error = lambda x, x_err : quadratic_err(x, x_err, a, a_err, b, b_err, c, c_err)

calibration_file = "calibration.py"

with open(calibration_file, 'w') as f:
    f.write("from models_2 import quadratic, quadratic_err\n")
    f.write("calibration = lambda x : quadratic(x, "+str(a)+","+str(b)+","+str(c)+")\n")
    f.write("calibration_error = lambda x, x_err : quadratic_err(x,x_err,"+str(a)+","+str(a_err)+","+str(b)+","+str(b_err)+","+str(c)+","+str(c_err)+")\n")

plt.errorbar(measured_reference, wavelengths['reference'], error_reference, label = "calibration data", color = 'b', fmt = '.')
x = np.linspace(min(measured_reference) - 100, max(measured_reference) + 100, 2010)
plt.plot(x, calibration(x), label = "calibration curve", color = 'r')
plt.text(5000, 4000, "Chi square: " + str(chisqr))
plt.text(5000, 3800, "Reduced: " + str(redchi))
plt.xlabel("measured wavelength")
plt.ylabel("actual wavelength")
plt.legend()
plt.show()
# File path to small peak data files
# Read data
# Voigt Model & Fit
#completed earlier

# Test quadratic model
# better way of doing this with arrays probably exists
plt.errorbar(measured_check, wavelengths['check'], error_check, label = "check data", color = 'b', fmt = '.')
x = np.linspace(min(measured_check) - 100, max(measured_check) + 100, 2010)
plt.plot(x, calibration(x), label = "calibration curve", color = 'r')
#plt.text(5000, 4000, "Chi square: " + str(chisqr))
#plt.text(5000, 3800, "Reduced: " + str(redchi))
plt.xlabel("measured wavelength")
plt.ylabel("actual wavelength")
plt.legend()
plt.show()
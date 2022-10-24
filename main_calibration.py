from re import M
from data.data_loader import read_data
from data_processing import process_data
from fitters import fit, execute_peak_fit, naive_fit
from models_2 import *
from noise_reduction import reduce_noise
import matplotlib.pyplot as plt
import lmfit
import numpy as np
from fitters import fit_to_voigt
import os, glob


#this might be buggy. I'm sorry, I totally ran out of energy


'''
# File Path to main peak data files
folder = "data/10.14.22/mercury/"

# Small lines
small_lines = [4358.328, 5677.105, 4347.494, 6149.475, 3131.55, 3131.84]


# Define Main Peak Expected Lines
wavelengths = {
    # https://physics.nist.gov/PhysRefData/Handbook/Tables/mercurytable2.htm
    # possibly include parameter for hyperfine splitting
    'reference': [5790.663, 5769.598, 5460.735, 4046.563, 3650.153],
    'check': []
    #check': [6562, 4861, 4340.47, 4101.74], #hydrogen
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
    #'check' : ["alpha-superfine-1", "beta-superfine-1", "gamma-superfine-1", "delta-superfine-1"],
    'check' : [],
}

for key in ['reference']:
    for i in range(len(files[key])):
        files[key][i] = folder + files[key][i]


for key in ['check']:
    for i in range(len(files[key])):
        files[key][i] = "data/10.14.22/hydrogen/" + files[key][i]

#files['reference'][4] = "data/mercury/calibration_scan_9_4019_4023"
'''



# noise-reduced version: 
"""
data = {
    'reference': [(reduce_noise(read_data(files['reference'][i]))) for i in range(len(files['reference']))],
    'check': [(reduce_noise(read_data(files['check'][i]))) for i in range(len(files['check']))],
}
"""


# File Path to main peak data files
folder_main = "data/final_data/mercury/main/"
folder_check = "data/final_data/mercury/check/"


#main peaks
main_3650_15 = {
    "fp" : folder_main + "3650_15",
    "true value" : 3650.15,
    "splitting" : 0.1,
    "damping" : 0.1,
}
main_4046_56 = {
    "fp" : folder_main + "4046_56",
    "true value" : 4046.56,
    "splitting" : 0.08,
    "damping" : 0.1,
}
main_5460_74 = {
    "fp" : folder_main + "5460_74",
    "true value" : 5460.74,
    "splitting" : 0.1,
    "damping" : 0.1,
}
main_5769_6 = {
    "fp" : folder_main + "5769_6",
    "true value" : 5769.6,
    "splitting" : 0.1,
    "damping" : 0.1,
}
main_5790_66 = {
    "fp" : folder_main + "5790_66",
    "true value" : 5790.66,
    "splitting" : 0.125,
    "damping" : 0.1,    
}

main_peaks = [main_3650_15, main_4046_56, main_5460_74, main_5769_6, main_5790_66]
#main_peaks = []

#check
check_3125_668 = {
    "fp" : folder_check + "3125_668",
    "true value" : 3125.668,
    "splitting" : 0.1,
    "damping" : 0.1,    
}
#THIS NEEDS TO BE SPLIT INTO SEVERAL
check_3131_548 = {
    "fp" : folder_check + "3131_548",
    "true value" : 3131.548,
    "splitting" : 0.3,
    "damping" : 0.1,    
}
check_3654_836 = {
    "fp" : folder_check + "3654_836",
    "true value" : 3654.836,
    "splitting" : 0.1,
    "damping" : 0.1,    
}
# also sus
check_4311_65 = {
    "fp" : folder_check + "4311_65",
    "true value" : 4311.65,
    "splitting" : 0.075,
    "damping" : 0.1,    
}
# also, also a little bit sus
check_4347_494 = {
    "fp" : folder_check + "4347_494",
    "true value" : 4347.494,
    "splitting" : 0.1,
    "damping" : 0.1,  
}
check_4358_328 = {
    "fp" : folder_check + "4358_328",
    "true value" : 4358.328,
    "splitting" : 0.15,
    "damping" : 0.1,  
}

#temporarily not including this one cuz it fails

check_5425_253 = {
    "fp" : folder_check + "5425_253",
    "true value" : 5425.253,
    "splitting" : 0,
    "damping" : 0.3,
}

check_5677_105 = {
    "fp" : folder_check + "5677_105",
    "true value" : 5677.105,
    "splitting" : 0,
    "damping" : 0.1
}
#temporarily not including this one cuz it fails
check_6149_475 = {
    "fp" : folder_check + "6149_475",
    "true value" : 6149.475,
    "splitting" : 0,
    "damping" : 0.3,
}


#check_peaks = [check_3125_668, check_3131_548, check_3654_836, check_4311_65, check_4347_494, check_4358_328, check_5425_253, check_5677_105, check_6149_475]
check_peaks = [check_3125_668, check_3131_548, check_3654_836, check_4311_65, check_4347_494, check_4358_328, check_5677_105]
sus_peaks = [check_3131_548, check_4311_65, check_4358_328]
#sus_peaks = [check_3131_548, check_4347_494]
maxmodel_peaks = [check_5425_253, check_6149_475]
'''
# non noise-reduced version:
data = {
    'reference' : [read_data(file) for file in files['reference']],
    'check' : [read_data(file) for file in files['check']]
}
'''

superfine_stepsize = 0.0025

# plot, user request splitting estimate

measured = {
    'reference': [],
    'check' : [],
}
error = {
    'reference': [],
    'check' : [],
}
true_wavelengths = {
    'reference': [],
    'check' : [],
}


folder_H = "data/final_data/hydrogen/"
H_alpha = {
    'fp' : folder_H + "alpha",
    'true value' : 6562.8518,
    "splitting" : 0.05,
    "damping" : 0.1
}
H_beta = {
    'fp' : folder_H + "beta",
    'true value' : 4861.2786,
    "splitting" : 0.05,
    "damping" : 0.1
}
H_gamma = {
    'fp' : folder_H + "gamma",
    'true value' : 4340.462,
    "splitting" : 0,
    "damping" : 0.1
}
H_delta = {
    'fp' : folder_H + "delta",
    'true value' : 4101.74,
    "splitting" : 0.25,
    "damping" : 0.1
}
H = [H_alpha, H_beta, H_gamma, H_delta]

"""
splittings = {
    'reference': [None, None, None, None, None],
    'check': [],
}
damping_constants = {
    'reference': [None, None, None, None, None],
    'check': [],
}
"""

# Voigt Model & Fit

#for key in data:

enter_splittings = False
check_fits = False

for wavelength in main_peaks + check_peaks + H:
#for wavelength in H:
    #for i in range(len(data[key])):
        #entry = data[key][i]

        #print(wavelength["fp"])

        entry = read_data(wavelength["fp"]) 

        #plt.plot(entry[:,0], entry[:,1])
        #plt.show()

        #Made some changes, if you think this is worse, you can revert though this might mean we should modify process_data to make it more appropriate -Athira
        processed_data,weights = process_data(entry, plot_noise_reduction=enter_splittings)  #Added by Athira, plot set to true so can answer shift questions
        #new_data, weights = reduce_noise(entry)        --commented out by Athira, to implement process_data
        #if enter_splittings or splittings[key][i] is None:
        if enter_splittings:
            '''
            plt.scatter(entry[:,0], entry[:,1], marker = '.', label = "data", color = 'b')
            plt.plot(entry[:,0], processed_data[:,1], label = "noise-reduced", color = 'r')
            plt.legend()
            plt.show()
            '''

            shift = float(input("Estimate the splitting in that plot (Angstroms): "))
            damping_constant_input = input("Give a damping constant (if adequate, just hit enter): ")
            damping_constant = 0.1 if damping_constant_input == "" else float(damping_constant_input)

            '''
            splittings[key][i] = shift
            damping_constants[key][i] = damping_constant
            '''
        else:
            shift = wavelength['splitting']
            damping_constant = wavelength['damping']
        #result = fit_to_voigt(processed_data, weights, shift = splittings[key][i], damping_constant = damping_constants[key][i], plot = enter_splittings)
        
        if wavelength in maxmodel_peaks:
            mu, mu_err = naive_fit(entry, damping_constant)
        else:
            result_params = fit_to_voigt(processed_data, weights, shift, damping_constant, plot = check_fits)
            mu, mu_err = result_params["mu"].value, result_params["mu"].stderr

        '''
        measured[key].append(mu)
        error[key].append(mu_err)
        '''

        print ('--')
        print (wavelength['true value'])
        print (mu)
        print (mu_err)
        
        if wavelength in H:
            measured['check'].append(mu)
            error['check'].append(mu_err)
            true_wavelengths['check'].append(wavelength["true value"])
        else:
            if wavelength not in sus_peaks:
                measured['reference'].append(mu)
                error['reference'].append(mu_err)
                true_wavelengths['reference'].append(wavelength["true value"])

                
                wavelength["measured value"] = mu
                wavelength["error"] = mu_err

print (measured["reference"])
print(len(main_peaks))
print(len(check_peaks))

# Maximum Counts
#??

# Find Residuals
# no need

# Quadratic Model and Fit
#plt.errorbar(measured['reference'], wavelengths['reference'], yerr = error['reference'], label = "calibration data", color = 'b', fmt = '.')
"""
plt.errorbar(measured['reference'], true_wavelengths['reference'], yerr = error['reference'], label = "calibration data", color = 'b', fmt = '.')
plt.xlabel("measured wavelength")
plt.ylabel("actual wavelength")
plt.legend()
plt.show()
"""


quad = True
sin = False
cube = False

if quad:
    calibration_model = lmfit.Model(quadratic)
    params = lmfit.Parameters()
    params.add('a', value = 0)
    params.add('b', value = 1)
    params.add('c', value = 0)
    calibration_error = quadratic_err

if sin:
    calibration_model = lmfit.Model(linear_plus_osc)
    params = lmfit.Parameters()
    params.add('a', value = 1)
    params.add('b', value = 0)
    params.add('n', value = 0.1)
    params.add('omega', value = 1/1000)
    params.add('phi', value = 0, min = -np.pi, max = np.pi)

if cube:
    calibration_model = lmfit.Model(cubic)
    params = lmfit.Parameters()
    params.add('a', value = 0)
    params.add('b', value = 0)
    params.add('c', value = 1)
    params.add('d', value = 0)
    



inputs = np.zeros((len(measured['reference']), 2))
inputs[:,0] = measured['reference']
#inputs[:,1] = wavelengths['reference']
inputs[:,1] = true_wavelengths['reference']

result = fit(calibration_model, inputs, params, weights = 1/np.array(error['reference'])) # a slight approximation of the true weights but whatever (the slope is roughly 1).

calibration_file = "calibration.py"

write = input("Write file? y/n ") == 'y'

if quad:
    a, b, c = result.params['a'].value, result.params['b'].value, result.params['c'].value

    calibration = lambda x : quadratic(x, a, b, c)
    a_err = result.params['a'].stderr
    b_err = result.params['b'].stderr
    c_err = result.params['c'].stderr
    calibration_error = lambda x, x_err : quadratic_err(x, x_err, a, a_err, b, b_err, c, c_err)
    if write:
        with open(calibration_file, 'w') as f:
            f.write("from models_2 import quadratic, quadratic_err, inverse_quadratic\n")
            f.write("calibration = lambda x : quadratic(x, "+str(a)+","+str(b)+","+str(c)+")\n")
            f.write("calibration_error = lambda x, x_err : quadratic_err(x,x_err,"+str(a)+","+str(a_err)+","+str(b)+","+str(b_err)+","+str(c)+","+str(c_err)+")\n")
            f.write("uncalibration = lambda true : inverse_quadratic(true,"+str(a)+","+str(b)+","+str(c)+")\n")

if sin:
    a, b, n, omega, phi = result.params['a'].value, result.params['b'].value, result.params['n'].value, result.params['omega'].value, result.params['phi'].value

    calibration = lambda x : linear_plus_osc(x, a, b, n, omega, phi)
    a_err = result.params['a'].stderr
    b_err = result.params['b'].stderr
    n_err = result.params['n'].stderr
    omega_err = result.params['omega'].stderr
    phi_err = result.params['phi'].stderr

    if write:
        with open(calibration_file, 'w') as f:
            f.write("from models_2 import linear_plus_osc\n")
            f.write("calibration = lambda x : linear_plus_osc(x, "+str(a)+", "+str(b)+", "+str(n)+", "+str(omega)+", "+str(phi)+")\n")

if cube:
    a, b, c, d = result.params['a'].value, result.params['b'].value, result.params['c'].value, result.params['d'].value

    calibration = lambda x : cubic(x, a, b, c, d)
    a_err = result.params['a'].stderr
    b_err = result.params['b'].stderr
    c_err = result.params['c'].stderr   
    d_err = result.params['d'].stderr

    if write:
        with open(calibration_file, 'w') as f:
            f.write("from models_2 import cubic, cubic_err\n")
            f.write("calibration = lambda x : cubic(x, "+str(a)+","+str(b)+","+str(c)+","+str(d)+")\n")
            f.write("calibration_error = lambda x, x_err : cubic(x,x_err,"+str(a)+","+str(a_err)+","+str(b)+","+str(b_err)+","+str(c)+","+str(c_err)+","+str(d)+","+str(d_err)+")\n")



chisqr = result.chisqr
redchi = result.redchi

#probably want to ignore the printed chi square values. I really don't understand them. 

#plt.errorbar(measured['reference'], wavelengths['reference'], error['reference'], label = "calibration data", color = 'b', fmt = '.')
plt.errorbar(measured['reference'], true_wavelengths['reference'], error['reference'], label = "calibration data", color = 'b', fmt = '.')
plt.errorbar(measured['check'], true_wavelengths['check'], error['check'], label = "Hydrogen data", color = 'magenta', fmt = '.')
x = np.linspace(min(min(measured['reference']), min(measured['check'])) - 100, max(max(measured['reference']), max(measured['check'])) + 100, 2010)
plt.plot(x, calibration(x), label = "calibration curve", color = 'r')
plt.text(5000, 4000, "Chi square: %.4f" % chisqr) # truncate decimal
plt.text(5000, 3800, "Reduced: %.4f" % redchi) # truncate decimal
plt.xlabel("measured wavelength")
plt.ylabel("actual wavelength")
plt.legend()
plt.show()
# File path to small peak data files
# Read data
# Voigt Model & Fit
# completed earlier
def chisquare(x, e):
    return np.dot(x-e, (x-e)/e)

print (chisquare(calibration(np.array(measured['reference'])), true_wavelengths['reference']))

print (chisquare(calibration(np.array(measured['check'])), true_wavelengths['check']))

#for i in range(len(measured['check'])):
#    print ("--")
#    print (measured['check'][i])
#    print (error['check'][i])
#    print (calibration_error(measured['check'][i], error['check'][i]))

#get chi square of performance on check data. 
"""
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
"""


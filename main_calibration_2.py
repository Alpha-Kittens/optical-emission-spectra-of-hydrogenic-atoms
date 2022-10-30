from re import M
from data.data_loader import read_data
from data_processing import process_data
from models_2 import *
import matplotlib.pyplot as plt
import lmfit
import numpy as np
from max_model import hwhm_max
from calibration_models import *


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
    "threshold" : 1/2,
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

def do_calibration(measured, true, weights, n_poly, include_exp, plot = True, check = []):
    params = lmfit.Parameters()
    poly_params(params, n_poly)
    if include_exp:
        c_model = lambda x, **params : model_poly(mean)(x, **params) + model_exp(mean)(x, **params)
        exp_params(params)
    else:
        c_model = model_poly(mean)

    model = lmfit.Model(c_model)
    result = model.fit(true, params=params, x=measured, weights=weights)
    #print (lmfit.fit_report(result))
    c_params, c_errs = extract(result.params)

    title = "Calibration model: degree-"+str(n_poly)+" polynomial"
    if include_exp:
        title += " plus exponential"
    print (title)
    print ("Chi square: "+str(result.chisqr))
    print ("Reduced: "+str(result.redchi))
    if plot:
        plt.title(title)
        plt.xlabel("Monochrometer reading (Ã)")
        plt.ylabel("True wavelength (A)")
        plt.errorbar(x = measured, y = true, xerr = 1/weights, marker = '.', ls = 'none', label = "Calibration data", c = 'b')
        xmin = min(measured)
        xmax = max(measured)
        if check != []:
            xmin = min(xmin, min(check[0]))
            xmax = max(xmax, max(check[0]))
            plt.errorbar(x = check[0], y = check[1], xerr = check[2], marker = '.', ls = 'none', label = "H data", c = 'magenta')
        x = np.linspace(xmin - 100, xmax + 100, 1000)
        plt.plot(x, c_model(x, **c_params), label = "Calibration curve", c = 'r')
        plt.text(xmin + 100, 5500, "Chi square: "+str(result.chisqr))
        plt.text(xmin + 100, 5300, "Reduced: "+str(result.redchi))
        plt.legend()
        plt.show()
        plt.title(title)
        plt.xlabel("Monochrometer reading (Ã)")
        plt.ylabel("Residual (A)")
        plt.errorbar(x = measured, y = c_model(np.array(measured), **c_params) - np.array(true), yerr = 1/weights, marker = '.', ls = 'none', label = "Calibration data", c = 'b')
        xmin = min(measured)
        xmax = max(measured)
        if check != []:
            xmin = min(xmin, min(check[0]))
            xmax = max(xmax, max(check[0]))
            plt.errorbar(x = check[0], y = c_model(np.array(check[0]), **c_params) - np.array(check[1]), yerr = check[2], marker = '.', ls = 'none', label = "H data", c = 'magenta')
        x = np.linspace(xmin - 100, xmax + 100, 1000)
        #plt.plot(x, c_model(x, **c_params), label = "Calibration curve", c = 'r')
        plt.axhline(y=0, c = 'r', label = "Calibration curve", linestyle  = '--')
        plt.text(xmin + 100, 5500, "Chi square: "+str(result.chisqr))
        plt.text(xmin + 100, 5300, "Reduced: "+str(result.redchi))
        plt.legend()
        plt.show()
    return c_params, c_errs, result.chisqr, result.redchi

def cwrite(file, key, mean, params, errs):
    # in retrospect, I could have just made the calibration models detect which model is used based on parameter names. But oh well. 
    with open(file, 'w') as f:
        f.write("from calibration_models import *\n")
        f.write("from numpy import sqrt, abs\n")
        f.write("params = " + str(params) + "\n")
        f.write("errors = " + str(errs) + "\n")
        if key[1] == 1:
            f.write("calibrate = lambda x : model_both("+str(mean)+")(x, **params)\n")
            f.write("calibrate_error = lambda x, x_err : model_both_err("+str(mean)+", params, errors)(x, x_err)\n")
        else:
            f.write("calibrate = lambda x : model_poly("+str(mean)+")(x, **params)\n")
            f.write("calibrate_error = lambda x, x_err : model_poly_err("+str(mean)+", params, errors)(x, x_err)\n")
        f.write("calibrate_splitting = lambda x1, x2 : calibrate(x2) - calibrate(x1)\n")
        f.write("calibrate_splitting_error = lambda x1, x2, x1_err, x2_err : (sqrt(calibrate_error(x1, x1_err)[0]**2 + calibrate_error(x2, x2_err)[0]**2), abs(calibrate_error(x1, x1_err)[1] - calibrate_error(x2, x2_err)[1]))\n")

            

if __name__ == '__main__':

    check_fits = False
    sus_peaks = [check_4311_65, check_5677_105]

    for wavelength in main_peaks + check_peaks + H:

            entry = read_data(wavelength["fp"]) 

            #plt.plot(entry[:,0], entry[:,1])
            #plt.show()
            #Made some changes, if you think this is worse, you can revert though this might mean we should modify process_data to make it more appropriate -Athira
            processed_data,weights, noise_reduced = process_data(entry, plot_noise_reduction=False, noise_reduced = True)  #Added by Athira, plot set to true so can answer shift questions
            #new_data, weights = reduce_noise(entry)        --commented out by Athira, to implement process_data
            #if enter_splittings or splittings[key][i] is None:
            
            # execute hwhm
            if wavelength not in H:
                if wavelength not in sus_peaks:
                    if 'threshold' in wavelength.keys():
                        threshold = wavelength['threshold']
                    else:
                        threshold = 1/3
                    mu, mu_err = hwhm_max(processed_data, weights, plot = check_fits, noise_reduced = noise_reduced, true_wavelength = wavelength['true value'], threshold = threshold)
                    measured['reference'].append(mu)
                    error['reference'].append(mu_err)
            else:
                if wavelength not in sus_peaks:
                    mu, mu_err = hwhm_max(processed_data, weights, plot = check_fits, noise_reduced = noise_reduced, true_wavelength = wavelength['true value'], threshold = 1/3)
                    measured['check'].append(mu)
                    error['check'].append(mu_err)

    print("Summary:")
    print ("Reference (Hg) peaks:")
    true_wavelengths = []
    offset = 0
    for i, wavelength in enumerate(main_peaks + check_peaks):
        if wavelength not in sus_peaks:
            true_wavelengths.append(wavelength['true value'])
            print ("True value: " + str(wavelength['true value']))
            print ("Estimated peak: " + str(measured['reference'][i - offset]))
            print ("Uncertainty: " + str(error['reference'][i - offset]))
            print ("---")
        else:
            offset += 1
    print ("Check (H) peaks:")
    true_H_wavelengths = []
    offset = 0
    for i, wavelength in enumerate(H):
        if wavelength not in sus_peaks:
            true_H_wavelengths.append(wavelength['true value'])
            print ("True value: " + str(wavelength['true value']))
            print ("Estimated peak: " + str(measured['check'][i - offset]))
            print ("Uncertainty: " + str(error['check'][i - offset]))
            print ("---")    
        else:
            offset += 1

    mean = np.mean(measured['reference'])
    weights = 1/np.array(error['reference'])

    check = (measured["check"], true_H_wavelengths, error['check'])

    x = np.linspace(min(min(measured['reference']),min(measured['check'])), max(max(measured['reference']),max(measured['check'])), 1000)

    results = {}
    for include_exp in (False, True):
        ekey = 1 if include_exp else 0
        for degree in range(1, 7):
            results[(ekey, degree)] = do_calibration(measured['reference'], true_wavelengths, weights, degree, include_exp, check = check, plot = True)
            params, errs, chisqr, redchi = results[(ekey, degree)] 
            #plt.plot(x, model_poly(mean)(x, results[(ekey, degree)]), label = "calibration curve", c = 'r')
            #plt.errorbar(measured['reference'], true_wavelengths, xerr = )

    minchi = 10
    mindetails = None
    for details, result in results.items():
        if result[3] < minchi:
            minchi = result[3]
            mindetails = details
    print (minchi)
    print (mindetails)

    winner = (0,3)
    slope = results[winner][0]['b']
    print (slope)
    params, errs, chisqr, redchi = do_calibration(measured['reference'], true_wavelengths, weights * slope, winner[1], winner[0] == 1, check = check, plot = True)
    cwrite("calibration_test.py", winner, mean, params, errs)
    #cwrite("calibration.py", winner, mean, results[winner][0], results[winner][1])



        
    #winner selection

    #file writing

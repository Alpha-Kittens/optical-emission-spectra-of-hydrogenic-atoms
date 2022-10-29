# Imports
from data import data_loader
import lmfit
from data_processing import get_cutoff
from fitters import fit, fit_to_voigt
import matplotlib.pyplot as plt
import models_2 as models
import noise_reduction
import math
from calibration import calibrate
from calibration import calibrate_error
from data_processing import process_data
from  max_model import hwhm_max
import numpy as np


# File Path to peak data files
folder = 'data/final_data/hydrogen/'
alpha_H_fp  = folder + 'alpha'
beta_H_fp = folder + 'beta'
gamma_H_fp = folder + 'gamma'
delta_H_fp = folder + 'delta'


#Constants
m_p = 1.67262192E-27 #kilograms
m_e = 9.1093837E-31 #kilograms
c = 2.998E8 #m/s^2
k_B = 1.38E-23 #boltzmann constant
nf = 2
R_inf = 10973731.568160 #m^-1

# Wavelengths (note that Hydrogen is the higher peak in isotope splitting)
alpha = {
    "fp" : alpha_H_fp,
    "ni"   : 3,
    "reference" : 6562.79
}
beta = {
    "fp" : beta_H_fp,
    "ni"   : 4,
    "reference" : 4861.35

}
gamma = {
    "fp" : gamma_H_fp,
    "ni"   : 5,
    "reference" : 4340.47

}
delta = {
    "fp" : delta_H_fp,
    "ni"   : 6,
    "reference" : 4101.74

}


information = {
    'alpha' : alpha,
    'beta' : beta,
    'gamma' : gamma,
    'delta' : delta,
}

damping_constants = {
    'alpha' : 1/10,
    'beta' : 1/3,
    'gamma' : 1/10,
    'delta' : 1/10
}

# Run fits for each wavelength in information
#for key in ['alpha']:
for key in information:
    print ("=======")
    print (key)

    wavelength = information[key]

    # Read data
    data = data_loader.read_data(wavelength["fp"])

    processed_data,weights, noise_reduced = process_data(data, plot_noise_reduction=True, noise_reduced=True)

    uncalibrated_wavelength, uncalibrated_uncertainty = hwhm_max(processed_data, weights, noise_reduced=noise_reduced, plot=True, threshold = 1/3)

    print ("true: "+str(wavelength['reference']))
    # Calibration to convert the wavelengths
    true_wavelength = calibrate(uncalibrated_wavelength)
    true_uncertainty_stat, true_uncertainty_sys = calibrate_error(uncalibrated_wavelength,uncalibrated_uncertainty)

    wavelength["wavelength"] = true_wavelength
    wavelength["wavelength_unc_sys"] = true_uncertainty_sys
    wavelength["wavelength_unc_stat"] = true_uncertainty_stat

    print("Systematic uncertainty: " + str(wavelength["wavelength_unc_sys"]))
    print("Statistical uncertainty: " + str(wavelength["wavelength_unc_stat"]))



    # Bohr Formula to find Rydhberg Const.
    wavelength["R_H"] = 1/(((wavelength["wavelength"])*(10**(-10)))*((1/(nf**2)) - (1/(wavelength["ni"]**2))))
    wavelength["R_H_unc_sys"] = (wavelength["R_H"]/(wavelength["wavelength"]))*(wavelength["wavelength_unc_sys"])
    wavelength["R_H_unc_stat"] = (wavelength["R_H"]/(wavelength["wavelength"]))*(wavelength["wavelength_unc_stat"])
    wavelength["R_H_unc_tot"] = math.sqrt(wavelength["R_H_unc_sys"]**2 + wavelength["R_H_unc_stat"]**2)


    # Rydberg Calculated from Reference values
    wavelength["ref_R_H"] = 1/(((wavelength["reference"])*(10**(-10)))*((1/(nf**2)) - (1/(wavelength["ni"]**2))))



    '''
    OLD CODE
    # Noise reduction and find Uncertainties
    new_data, weights = noise_reduction.reduce_noise(data, damping_constants[key]) # damping constant default is 1/10


    # Voigt Model & Fit
    x_axis = data[:, 0]
    y_axis = data[:, 1]

    #model = lmfit.models.VoigtModel()
    model = lmfit.Model(models.voigt_with_shift)
    #model = models.two_voigts() #not working

    #params = models.voigt_params(model, data)
    params = models.voigt_shift_params(data)
    #params = models.two_voigt_params(model, data) #not working

    result = fit(model, data, params, weights)
    #result = fit(model, data, params)

    #lmfit.report_fit(result)
        
    print(lmfit.fit_report(result))


    result.plot(datafmt = '.')
    plt.plot(x_axis, new_data[:,1], '-', color = 'r', label = 'noise-reduced')
    plt.legend()
    plt.show()
    '''



# Weighted Average
numerator = 0
denominator = 0
for key in information:
    wavelength = information[key]
    numerator += wavelength["R_H"] / ((wavelength["R_H_unc_tot"])**2)
    denominator += 1/((wavelength["R_H_unc_tot"])**2)

mean = numerator/denominator
print(mean)

uncertainty = 1/math.sqrt(denominator)

#Standard Deviation
sum = 0
for key in information:
    wavelength = information[key]
    sum += (wavelength["R_H"] - mean)**2
std_dev = math.sqrt(sum/(len(information)))

std_error = std_dev/math.sqrt(len(information))

results = {
    "mean" : mean,
    "uncertainty" : uncertainty,
    "standard deviation": std_dev,
    "standard error" : std_error
}

print(results)


#Make a plot
n = []
wavelengths = []
errors = []
references = []
residuals = []

for key in information:
    wavelength = information[key]
    n.append(wavelength["ni"])
    wavelengths.append(wavelength["wavelength"])
    references.append(wavelength["reference"])
    residuals.append(wavelength["wavelength"] - wavelength["reference"])
    errors.append(math.sqrt(wavelength["wavelength_unc_stat"]**2 + wavelength["wavelength_unc_sys"]**2))

#plt.errorbar(n, wavelengths, errors, fmt='o')
#plt.errorbar(n, references, fmt='o')
plt.errorbar(n, residuals,errors, fmt='o')
plt.axhline(y = 0, c = 'r', ls = '--')
plt.xlabel('initial n')
plt.ylabel('(measured wavelength - reference value) (A)')
plt.show()


measured_rydbergs = []
true_rydbergs = []
rydberg_error = []
sum = 0
for key in information:
    wavelength = information[key]
    sum += wavelength["ref_R_H"]
    true_rydbergs.append(wavelength["ref_R_H"])
    measured_rydbergs.append(wavelength["R_H"])

# Average ref rydberg
ref_rydberg = sum/len(information)

print('Reference rydberg: ' + str(ref_rydberg))


# Rydberg Model
rydberg_model = lmfit.Model(models.rydberg_model)

weights = []
for i in errors:
    weights.append(1/(i * (1e-10)))

balmer_data = np.zeros(shape=(len(n), 2))
for i in range(0, len(n)):
    balmer_data[i, 0] = n[i]
    balmer_data[i, 1] = wavelengths[i] * (1e-10)


result = fit(model=rydberg_model, data= balmer_data, weights=weights)

# Errors rescaled
errors_m = []
for i in errors:
    errors_m.append(i * (1e-8)) # errors times 100

plt.errorbar(balmer_data[:, 0], balmer_data[:, 1], errors_m, fmt='o')
x=np.linspace(min(balmer_data[:,0]), max(balmer_data[:,0]), 1000)
y=[]
for i in x:
    y.append(models.rydberg_model(i, result.best_values['a']))
plt.plot(x, y)
plt.show()


print(lmfit.fit_report(result))

balmer_ref_data = np.zeros(shape=(len(n), 2))
for i in range(0, len(n)):
    balmer_ref_data[i, 0] = n[i]
    balmer_ref_data[i, 1] = references[i] * (1e-10)


result = fit(model=rydberg_model, data= balmer_ref_data)
print(lmfit.fit_report(result))




# Write to a file (to be read by isotope shift)
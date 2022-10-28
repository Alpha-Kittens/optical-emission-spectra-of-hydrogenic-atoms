# Imports
from data import data_loader
import lmfit
from data_processing import get_cutoff
from fitters import fit, fit_to_voigt
import matplotlib.pyplot as plt
import models_2 as models
import noise_reduction
import math
from calibration import calibration
from calibration import calibration_error
from data_processing import process_data
from  max_model import hwhm_max

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
}
beta = {
    "fp" : beta_H_fp,
    "ni"   : 4,
}
gamma = {
    "fp" : gamma_H_fp,
    "ni"   : 5,
}
delta = {
    "fp" : delta_H_fp,
    "ni"   : 6,
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
for key in information:

    wavelength = information[key]

    # Read data
    data = data_loader.read_data(wavelength["fp"])

    processed_data,weights = process_data(data, plot_noise_reduction=True)

    uncalibrated_wavelength, uncalibrated_uncertainty = hwhm_max(processed_data, weights, plot=True)


    # Calibration to convert the wavelengths
    true_wavelength = calibration(uncalibrated_wavelength)
    true_uncertainty = calibration_error(uncalibrated_wavelength,uncalibrated_uncertainty)

    wavelength["wavelength"] = true_wavelength
    wavelength["wavelength_unc"] = true_uncertainty

    


    # Bohr Formula to find Rydberg Const.
    wavelength["R_H"] = 1/(((wavelength["wavelength"])*(10**(-10)))*((1/(nf**2)) - (1/(wavelength["ni"]**2))))
    wavelength["R_H  unc"] = (wavelength["R_H"]/(wavelength["wavelength"]))*(wavelength["wavelength_unc"])


    # Determining the Temperature
    alpha = calibration(wavelength["params"]['alpha'].value + wavelength["params"]['mu'].value) - true_wavelength

    constant_term = ((m_p)*(c**2)/(k_B))
    ratio_term = alpha/(true_wavelength)

    wavelength["temperature"] = constant_term*(ratio_term**2)/(2*math.log(2))
    print('temperature: ' + str(wavelength["temperature"]))



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
    numerator += wavelength["R_H"] / ((wavelength["R_H  unc"])**2)
    denominator += 1/((wavelength["R_H  unc"])**2)

mean = numerator/denominator
print(mean)

uncertainty = 1/math.sqrt(denominator)

results = {
    "mean" : mean,
    "uncertainty" : uncertainty,
}

print(results)


# Write to a file (to be read by isotope shift)
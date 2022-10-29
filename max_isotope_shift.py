# Imports
import wave
from data import data_loader
import lmfit
from fitters import fit
import matplotlib.pyplot as plt
import models_2 as models
import noise_reduction
import math
from fitters import fit_to_voigt
from calibration import calibrate, calibrate_error, calibrate_splitting, calibrate_splitting_error
from data_processing import process_data
from  max_model import hwhm_max




#Constants
m_p = 1.67262192E-27 #kilograms
m_e = 9.1093837E-31 #kilograms
nf = 2
R_inf = 10973731.568160 #m^-1


# File Path to peak data files
folder = 'data/final_data/deuterium/'


# Properties for each line
alpha = {
    "reference_hydrogen" : 6562.79,
    "ni"   : 3,
    "damping constant": 1/10
}
beta = {
    "reference_hydrogen" : 4861.35,
    "ni"   : 4,
    "damping constant": 1/10
}
gamma = {
    "reference_hydrogen" : 4340.47,
    "ni"   : 5,
    "damping constant": 1/10
}
delta = {
    "reference_hydrogen" : 4101.74,
    "ni"   : 6,
    "damping constant": 1/10
}

# A dictionary with all the lines
information = {
    'alpha' : alpha,
    'beta' : beta,
    'gamma' : gamma,
#    'delta' : delta,
}


for key in information:
    wavelength = information[key]

    wavelength['fp_D'] = folder + key + 'D'
    wavelength['fp_H'] = folder + key + 'H'


    if(key == "delta"):
        wavelength['fp'] = folder + key

    # Read data
    data_D = data_loader.read_data(wavelength["fp_D"])
    data_H = data_loader.read_data(wavelength["fp_H"])

    #Deuterium
    processed_data_D, weights_D, noise_reduced_D = process_data(data_D, plot_noise_reduction=True, noise_reduced=True)
    uncalibrated_wavelength_D, uncalibrated_uncertainty_D = hwhm_max(processed_data_D, weights_D, noise_reduced= noise_reduced_D, plot=True)

    #Hydrogen
    processed_data_H, weights_H, noise_reduced_H = process_data(data_H, plot_noise_reduction=True, noise_reduced=True)
    uncalibrated_wavelength_H, uncalibrated_uncertainty_H = hwhm_max(processed_data_H, weights_H, noise_reduced= noise_reduced_H, plot=True)


    # Calibration
    wavelength["deuterium"] = calibrate(uncalibrated_wavelength_D)
    wavelength['deuterium_error_stat'] = calibrate_error(uncalibrated_wavelength_D, uncalibrated_wavelength_D)[0]
    wavelength['deuterium_error_sys'] = calibrate_error(uncalibrated_wavelength_D, uncalibrated_wavelength_D)[1]
    print('deuterium: ' + str(wavelength["deuterium"]))

    wavelength['hydrogen'] = calibrate(uncalibrated_wavelength_H)
    wavelength['hydrogen_error_stat'] = calibrate_error(uncalibrated_wavelength_H, uncalibrated_wavelength_H)[0]
    wavelength['hydrogen_error_sys'] = calibrate_error(uncalibrated_wavelength_H, uncalibrated_wavelength_H)[1]
    print('hydrogen: ' + str(wavelength["hydrogen"]))


    #Shift
    wavelength["shift"] = calibrate_splitting(uncalibrated_wavelength_D, uncalibrated_wavelength_H)
    wavelength["shift_error_stat"] = calibrate_splitting_error(uncalibrated_uncertainty_D, uncalibrated_uncertainty_H, uncalibrated_uncertainty_D, uncalibrated_uncertainty_H)[0]
    wavelength["shift_error_sys"] = calibrate_splitting_error(uncalibrated_uncertainty_D, uncalibrated_uncertainty_H, uncalibrated_uncertainty_D, uncalibrated_uncertainty_H)[1]

    print('shift: ' + str(wavelength["shift"]))

    # Balmer Formula to find ratio of deuterium to proton mass
    wavelength["ratio"] = 1/(1 - R_inf*(m_p/m_e)*wavelength["shift"]*(10**(-10))*((1/(nf**2)) - (1/(wavelength["ni"]**2))))
    wavelength["ratio_error_stat"] = ((wavelength["ratio"])**2)*R_inf*(m_p/m_e)*((1/(nf**2)) - (1/(wavelength["ni"]**2))) * wavelength["shift_error_stat"]*(10**(-10))
    wavelength["ratio_error_sys"] = ((wavelength["ratio"])**2)*R_inf*(m_p/m_e)*((1/(nf**2)) - (1/(wavelength["ni"]**2))) * wavelength["shift_error_sys"]*(10**(-10))
    wavelength['ratio_error_tot'] = math.sqrt(wavelength["ratio_error_stat"]**2 + wavelength["ratio_error_sys"]**2)

    print("ratio: " + str(wavelength['ratio']))


# Mean (Weighted Average)
numerator = 0
denominator = 0
for key in information:
    wavelength = information[key]
    numerator += wavelength["ratio"] / ((wavelength["ratio_error_tot"])**2)
    denominator += 1/((wavelength["ratio_error_tot"])**2)
    print(wavelength["ratio_error_tot"])
    print(denominator)

mean = numerator/denominator
print(mean)

uncertainty = 1/math.sqrt(denominator)

results = {
    "mean" : mean,
    "uncertainty" : uncertainty,
}

print(results)

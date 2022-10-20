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
from calibration import calibration, calibration_error
from data_processing import process_data



#Constants
m_p = 1.67262192E-27 #kilograms
m_e = 9.1093837E-31 #kilograms
nf = 2
R_inf = 10973731.568160 #m^-1


# File Path to peak data files
folder = 'data/final_data/deuterium/'


# Properties for each line
alpha = {
    "hydrogen" : 6562.79,
    "ni"   : 3,
    "damping constant": 1/10
}
beta = {
    "hydrogen" : 4861.35,
    "ni"   : 4,
    "damping constant": 1/10
}
gamma = {
    "hydrogen" : 4340.47,
    "ni"   : 5,
    "damping constant": 1/10
}
delta = {
    "hydrogen" : 4101.74,
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

    wavelength['fp'] = folder + key + 'D'

    if(key == "delta"):
        wavelength['fp'] = folder + key

    # Read data
    data = data_loader.read_data(wavelength["fp"])

    #Process data
    processed_data, weights = process_data(data, plot_noise_reduction=True)


    #Fitting
    result_params  = fit_to_voigt(processed_data, weights, plot=True)

    # Calculating the Shift
    wavelength["deuterium"] = calibration(result_params["mu"].value)
    wavelength["shift"] = wavelength["hydrogen"] - wavelength["deuterium"]

    wavelength["shift_error"] = calibration_error(result_params["mu"].value, result_params["mu"].stderr)

    # Balmer Formula to find ratio of deuterium to proton mass
    wavelength["ratio"] = 1/(1 - R_inf*(m_p/m_e)*wavelength["shift"]*(10**(-10))*((1/(nf**2)) - (1/(wavelength["ni"]**2))))

    wavelength["ratio_error"] = ((wavelength["ratio"])**2)*R_inf*(m_p/m_e)*((1/(nf**2)) - (1/(wavelength["ni"]**2))) * wavelength["shift_error"]*(10**(-10))

# Mean (Weighted Average)
numerator = 0
denominator = 0
for key in information:
    wavelength = information[key]
    numerator += wavelength["ratio"] / ((wavelength["ratio_error"])**2)
    denominator += 1/((wavelength["ratio_error"])**2)

mean = numerator/denominator
print(mean)

# Uncertainty (currently this is wrong, pls ignore)
sum = 0
for key in information:
    wavelength = information[key]
    sum += (wavelength["ratio"] - mean)**2

std_dev = math.sqrt(sum/(len(information) - 1))

stat_unc = std_dev/math.sqrt(len(information))


results = {
    "mean" : mean,
    "standard deviation" : std_dev,
    "stat uncertainty" : stat_unc,
}

print(results)
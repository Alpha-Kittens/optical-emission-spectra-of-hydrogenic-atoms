from cmath import sqrt
from operator import truediv
from calibration import *
import numpy as np
from data import data_loader
from data_processing import process_data
from max_model import hwhm_max
import matplotlib.pyplot as plt

"""
Expected peaks:
h:  4494.18-4497.66 => 4961.200458672298-4964.64922227521
Measured range of interest:
4464.5-4466 | 4468-4470
g:  4664.811-4668.560 => 4634.672076705228-4638.390033986312
Measured range of interest:
4633.5-4635 | 4637.5-4639
f:  4748-4751.9 => 4717.165261188159-4721.032281414788
Measured range of interest:
4716-4717.75 | 4720-4721.5
e:  4978.541-4982.813 => 4945.701518371859-4949.935310465465
Measured range of interest:
4945-4946.1 | 4949-4950.5
d:  5148.838-5153.402 => 5114.445582941578-5118.9671221146655
Measured range of interest:
5113.4-5114.6 | 5118-5119.5
++5115.9-5117.3
c:  5682.633-5688.193-5688.205 => 5642.97692960492-5648.4908234835775
Measured range of interest:
5642-5643.5 | 5647.5 5649
b:  5889.950-5895.924 => 5848.088048059503-5853.997146099862
Measured range of interest:
5846-5850 | 5852-5855.5
a:  6154.225-6160.747 => 6109.419360517577-6115.867875353023
Measured range of interest:
6108-6109 | 6114.8-6116.3
"""

# constants
h_val = 4.135667696e-15
c_val = 299792458
en = lambda x : h_val*c_val/(x * 1e-10)
w = en
en_err = lambda x, x_err : (h_val*c_val/(x * 1e-10)) * x_err/x
w_err = en_err

# File Path to peak data files
folder = 'data/final_data/sodium/'

a = {
    'transition' : "7s --> 4p",
    "n" : 7
}
b = {
    'transition' : "4p --> 3s",
    "n" : 3
}
c = {
    'transition' : "5d --> 4p",
    'n': 5
}
d = {
    'transition' : "6s --> 4p",
    'n' : 6
}
e = {
    'transition' : "6d --> 4p",
    'n' : 6
}
f = {
    'transition' : "7s --> 4p",
    'n' : 7
}
g = {
    'transition' : "7d --> 4p",
    'n' : 7
}
h = {
    'transition' : "",
    'n' : 0
}

#Collect all the splittings
information = {
    'a' : a,
    'b' : b,
    "c" : c,
    "d" : d,
    "e" : e,
    "f" : f,
    "g" : g,
    "h" : h
}

slicent_peaks = [d]


for key in information:
    print ("=============")

    wavelength = information[key]
    wavelength["fp_1"] = folder + key + "L"
    wavelength["fp_2"] = folder + key + "H"


    energies = [[0,0], [0,0]]
    werr = [0,0]
    for i in range(1, 3):

    #for i in [2]:
        # Read data
        data = data_loader.read_data(wavelength["fp_" + str(i)])

        #Process data
        processed_data, weights, noise_reduced = process_data(data, plot_noise_reduction=False, noise_reduced=True, title=wavelength["transition"], slice = wavelength not in slicent_peaks)

        # Obtain wavelength and error from max model
        uncalibrated_wavelength, uncalibrated_uncertainty = hwhm_max(processed_data, weights, plot=False, noise_reduced=noise_reduced)
        print ("uncalibrated: "+str(uncalibrated_wavelength))

        # Calibration
        wavelength["unc_wavelength_" + str(i)] = uncalibrated_wavelength
        wavelength["unc_uncertainty_" + str(i)] = uncalibrated_uncertainty
        werr[i-1] = uncalibrated_uncertainty


        unc = calibrate_error(uncalibrated_wavelength, uncalibrated_uncertainty)
        wavelength["wavelength_" + str(i)] = calibrate(uncalibrated_wavelength)
        wavelength["stat_uncertainty_" + str(i)], wavelength["sys_uncertainty_" + str(i)] = unc


        
    

        energies[0][i-1] = en(calibrate(uncalibrated_wavelength))
        print ("energy: "+str(energies[0][i-1]))
        energies[1][i-1] = [en_err(calibrate(uncalibrated_wavelength), calibrate_error(uncalibrated_wavelength, uncalibrated_uncertainty)[j]) for j in range(2)]
    
    print ("===")
    
    print (energies)

    #print (unc[0])
    print (calibrate_splitting_error(wavelength["wavelength_2"], wavelength["wavelength_1"], werr[1], werr[0]))
    #print (energies[1])
    e_split = energies[0][0] - energies[0][1]
    e_err = np.sqrt(energies[1][1][0]**2 + energies[1][0][0]**2), abs(energies[1][1][1] - energies[1][1][0])
    e_err_tot = np.sqrt(e_err[0]**2+e_err[1]**2)

    wavelength["e"] = e_split
    wavelength["e_err"] = e_err
    wavelength["e_err_tot"] = e_err_tot
    wavelength["e_err_sys"] = e_err[1]

    print ("results:")
    print (e_split)
    print (e_err)
    print (e_err_tot)
    print (e_err_tot/e_split)
    #wavelength['splitting'] = calibrate_splitting(wavelength["unc_wavelength_1"], wavelength["unc_wavelength_2"])
    #wavelength["stat_splitting_uncertainty"] = calibrate_splitting_error(wavelength["unc_wavelength_1"], wavelength["unc_wavelength_2"], wavelength["unc_uncertainty_1"], wavelength["unc_uncertainty_2"])[0]
    #wavelength["sys_splitting_uncertainty"] = calibrate_splitting_error(wavelength["unc_wavelength_1"], wavelength["unc_wavelength_2"], wavelength["unc_uncertainty_1"], wavelength["unc_uncertainty_2"])[1]
    #wavelength["splitting_uncertainty_tot"] = sqrt(wavelength["stat_splitting_uncertainty"]**2 + wavelength["sys_splitting_uncertainty"]**2)

    #print("Splitting: " + str(wavelength['splitting']))
    #print("Spliting Error: " + str(wavelength['splitting_uncertainty_tot']))

#Weigted Average
numerator = 0
denominator = 0
for key in information:
    wavelength = information[key]
    numerator += wavelength['e'] /(wavelength['e_err_tot'] **2)
    denominator += 1/(wavelength['e_err_tot'] **2)

weighted_average = numerator/denominator
weighted_uncertainty = 1/sqrt(denominator)

print("Average splitting: " + str(weighted_average))
print("Uncertainty: " + str(weighted_uncertainty))


# Plotting
wavelengths = []
splittings = []
errors = []
n = []
syserrors = []

for key in information:
    wavelength = information[key]
    wavelengths.append(wavelength["wavelength_1"])
    splittings.append(wavelength["e"])
    errors.append(wavelength["e_err_tot"])
    syserrors.append(wavelength["e_err_sys"])
    n.append(wavelength['n'])

true = en(5889.950) - en(5895.924)

plt.title("Measurements of 3P fine structure splitting of sodium")
plt.xlabel("~Monochromator reading of doublet (Ã)")
plt.ylabel("Energy splitting (eV)")
plt.axhline(y = true, c = 'r', label = "True value: %6f eV" % true, ls = '--')
plt.errorbar(wavelengths, splittings, errors, fmt='o', color = 'blue', label = "Total uncertainties")
plt.errorbar(wavelengths, splittings, syserrors, fmt = 'o', label = "Systematic uncertainties", color = 'orange')
plt.text(4500, 0.001, "Average value: (%.6f ± " % weighted_average + "%6f) eV" % weighted_uncertainty, fontsize = 12)
plt.legend(loc = "upper left")
plt.show()

print (true)
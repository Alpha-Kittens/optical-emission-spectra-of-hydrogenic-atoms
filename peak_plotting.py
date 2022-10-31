from cProfile import label
from cmath import pi, sqrt
from calibration import *
import numpy as np
from data import data_loader
from data_processing import process_data
from fitters import fit_to_voigt, fit_to_lorentzian, fit_to_Gaussian
from max_model import hwhm_max
from max_model import check_against_voigt_pretty
import matplotlib.pyplot as plt
import os
import lmfit

'''
This file is created for making plots of things you want to show in presentation/paper
'''




###
# A Hydrogen Line
###

planck = 6.626e-34
c = 2.998e8
to_eV = 6.242e18


# An arbritrary peak
fp_alpha = 'data/final_data/hydrogen/alpha'

#Plot raw data
data = data_loader.read_data(fp_alpha)
#plt.plot(data[:, 0], data[:, 1])
#plt.title("Hydrogen alpha = 6562.79 A")
#plt.xlabel('Uncalibrated Angstroms')
#plt.ylabel('Counts per second')
'''
processed_data, weights, noise_reduced = process_data(data, plot_noise_reduction=True, noise_reduced = True, title="Hydrogen alpha = 6562.79 A")

hwhm_max(processed_data, weights, noise_reduced, plot=True, true_wavelength='Hydrogen alpha = 6562.79 A')
'''


processed_data, weights = process_data(data, plot_noise_reduction=True, title="Hydrogen alpha = 6562.79 A")

result_params = fit_to_voigt(processed_data, weights, plot=True)

print(result_params['mu'].value)

def to_Energy(wavelength):
    energy = (planck * c/(wavelength * (1e-10)))*to_eV

    return energy


line_width = to_Energy(result_params['mu'].value - result_params['gamma'].value) - to_Energy(result_params['mu'].value + result_params['gamma'].value)


print("Line Width: " + str(line_width))

result = fit_to_lorentzian(processed_data, weights)

result_2 = fit_to_Gaussian(processed_data, weights)




fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Hydrogen alpha = 6562.79 A')
ax1.plot(data[:, 0], data[:, 1])
ax1.set(xlabel='Uncalibrated Angstroms', ylabel='Counts per second')

energies = []
for i in range(0, len(data[:, 0])):
    energy = planck * c / (data[i, 0] * (1e-10))
    energies.append(energy)

#ax2.plot(energies, data[:, 1])
#ax2.set(xlabel='Uncalibrated Energy', ylabel='Counts per second')

ax2.plot(data[:,0], energies)
ax2.set(xlabel='Uncalibrated Wavelength', ylabel='Uncalibrated Energy')



plt.show()





###
# PLOTTING DEUTREIUM SCANS
###
#Plots of alpha from deuterium
'''
fp_alpha1 = 'data/final_data/interesting/alpha7L'
fp_alpha2 = 'data/final_data/interesting/alpha6L'
fp_alpha3 = 'data/final_data/interesting/alpha5L'
fp_alpha4 = 'data/final_data/interesting/alpha3L'
fp_alpha5 = 'data/final_data/interesting/alpha2L'
fp_alpha6 = 'data/final_data/interesting/alpha1L'

alpha_scans = [fp_alpha1, fp_alpha2, fp_alpha3, fp_alpha4, fp_alpha5, fp_alpha6]
alpha_scans = [fp_alpha1, fp_alpha2, fp_alpha3, fp_alpha6]


plt.title("Deuterium alpha = 6562.79 A")
plt.xlabel('Uncalibrated Angstroms')
plt.ylabel('Counts per second')
for i in range(len(alpha_scans)):
    data = data_loader.read_data(alpha_scans[i])
    if(i < 2):
        plt.plot(data[:, 0], data[:, 1], label='day 1 - scan ' + str(i + 1))
    else:
        plt.plot(data[:, 0], data[:, 1], label='day 2 - scan ' + str(i - 1))


plt.legend()
plt.show()
'''


###
# PLOTTING MERCURY
###
'''
fp_mercury1 = 'data/final_data/interesting/mercury_4046_56_bulb1'
fp_mercury2 = 'data/final_data/interesting/mercury_4046_56_bulb2'


mercury_scans_4046_56 = [fp_mercury1, fp_mercury2]

plt.title("Mercury scans for 4046.56 A Line")
plt.xlabel('Uncalibrated Angstroms')
plt.ylabel('Counts per second')
for i in range(len(mercury_scans_4046_56)):
    data = data_loader.read_data(mercury_scans_4046_56[i])
    plt.plot(data[:, 0], data[:, 1], label="bulb " + str(i+1))


plt.legend()
plt.show()
'''

'''
fp_mercury1 = 'data/final_data/interesting/mercury_3650_15_1'
fp_mercury2 = 'data/final_data/interesting/mercury_3650_15_2'


mercury_scans_3650_15 = [fp_mercury1, fp_mercury2]

plt.title("Mercury scans for 3650.15 A Line")
plt.xlabel('Uncalibrated Angstroms')
plt.ylabel('Counts per second')
for i in range(len(mercury_scans_3650_15)):
    data = data_loader.read_data(mercury_scans_3650_15[i])
    if i==1:
        plt.plot(data[:, 0], 10*data[:, 1], label="scan " + str(i+1) + "(cps * 10)")
    else:
        plt.plot(data[:, 0], data[:, 1], label="scan " + str(i+1))




plt.legend()
plt.show()
'''

'''
#An interesting glitch in fine structure b scan
fp_sodium_1 = 'data/final_data/interesting/sodium_bfine'
fp_sodium_2 = 'data/final_data/interesting/sodium_bL'
fp_sodium_3 = 'data/final_data/interesting/sodium_bH'

sodium_b_scans = [fp_sodium_1,fp_sodium_2,fp_sodium_3]

plt.title("Sodium b scans")
plt.xlabel('Uncalibrated Angstroms')
plt.ylabel('Counts per second')
for i in range(len(sodium_b_scans)):
    data = data_loader.read_data(sodium_b_scans[i])
    plt.plot(data[:, 0], data[:, 1], label="scan " + str(i+1))


plt.legend()
plt.show()
'''


'''
#Plotting all mercury
directory = 'data/final_data/mercury/'
folders = os.listdir(directory)


for folder in folders:
    folder_fp = directory + folder + '/'
    files = os.listdir(folder_fp)
    for file in files:
        fp = folder_fp + file
        data = data_loader.read_data(fp)
        plt.plot(data[:, 0], data[:, 1], label=folder)
        #plt.title(file)
plt.legend()
plt.show()
'''

# Plotting splits
'''
split_1 = 'data/final_data/interesting/mercury_3131_548'
data = data_loader.read_data(split_1)
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('Uncalibrated Angstroms')
plt.ylabel('Counts per second')
plt.title("Mercury Scan for 3131.48 A Line")
plt.show()

processed_data, weights = process_data(data, plot_noise_reduction=True)
fit_to_voigt(processed_data, weights, shift=0.4, plot=True)
'''

'''
split_1 = 'data/final_data/interesting/mercury_3650_15_1'
data = data_loader.read_data(split_1)
#plt.plot(data[:, 0], data[:, 1])
#plt.xlabel('Uncalibrated Angstroms')
#plt.ylabel('Counts per second')
#plt.title("Mercury Scan for 3650.15 A Line")
#plt.show()

processed_data, weights, noise_reduced = process_data(data, noise_reduced=True, plot_noise_reduction=True)
#fit_to_voigt(processed_data, weights, shift=0, plot=True)
#fit_to_voigt(processed_data, weights, shift=0.2, plot=True)

check_against_voigt_pretty(processed_data, weights, shift=0.2, noise_reduced=noise_reduced, true_wavelength= 'Mercury 3650.15 A line')



# Plotting bumpy peaks
'''
'''
folder = 'data/final_data/interesting/'
bumpy_1 = folder + 'sodium_eH'
bumpy_2 = folder + 'sodium_eL'
bumpy_3 = folder + 'mercury_4311_65'
bumpy_4 = folder + 'sodium_bfine'
bumpy_5 = folder + 'sodium_aH'

bumpy_peaks = [bumpy_1, bumpy_2, bumpy_3, bumpy_4, bumpy_5]

for i in bumpy_peaks:
    data = data_loader.read_data(i)

    processsed_data, weights = process_data(data, plot_noise_reduction=True)
'''
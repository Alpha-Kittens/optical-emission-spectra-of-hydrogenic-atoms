from cmath import sqrt
from calibration import *
import numpy as np
from data import data_loader
from data_processing import process_data
from max_model import hwhm_max
import matplotlib.pyplot as plt

'''
This file is created for making plots of things you want to show in presentation/paper
'''

'''
# An arbritrary peak
fp_alpha = 'data/final_data/hydrogen/alpha'

#Plot raw data
data = data_loader.read_data(fp_alpha)
plt.plot(data[:, 0], data[:, 1])
plt.title("Hydrogen alpha = 6562.79 A")
plt.xlabel('Uncalibrated Angstroms')
plt.ylabel('Counts per second')
plt.show()

process_data(data, plot_noise_reduction=True, title="Hydrogen alpha = 6562.79 A")
'''

#Two plots of alpha from deuterium
fp_alpha1 = 'data/final_data/interesting/alpha7L'
fp_alpha2 = 'data/final_data/interesting/alpha6L'
fp_alpha3 = 'data/final_data/interesting/alpha5L'
fp_alpha4 = 'data/final_data/interesting/alpha3L'
fp_alpha5 = 'data/final_data/interesting/alpha2L'
fp_alpha6 = 'data/final_data/interesting/alpha1L'

alpha_scans = [fp_alpha1, fp_alpha2, fp_alpha3, fp_alpha4, fp_alpha5, fp_alpha6]

plt.title("Deuterium alpha = 6562.79 A")
plt.xlabel('Uncalibrated Angstroms')
plt.ylabel('Counts per second')
for i in alpha_scans:
    data = data_loader.read_data(i)
    plt.plot(data[:, 0], data[:, 1])

plt.show()



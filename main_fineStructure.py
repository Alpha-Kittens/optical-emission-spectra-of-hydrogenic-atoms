from calibration import *
import numpy as np
from data import data_loader
from data_processing import process_data
from fitters import fit_to_voigt

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

# File Path to peak data files
folder = 'data/final_data/sodium/'

a = {
    'transition' : ""
}
b = {
    'transition' : ""
}
c = {
    'transition' : ""
}
d = {
    'transition' : ""
}
e = {
    'transition' : ""
}
f = {
    'transition' : ""
}
g = {
    'transition' : ""
}
h = {
    'transition' : ""
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

for key in information:

    wavelength = information[key]
    wavelength["fp_1"] = folder + key + "L"
    wavelength["fp_2"] = folder + key + "H"

    for i in range(1, 3):
        # Read data
        data = data_loader.read_data(wavelength["fp_" + str(i)])

        #Process data
        processed_data, weights = process_data(data, plot_noise_reduction=True)


        #Fitting
        result_params  = fit_to_voigt(processed_data, weights, plot=True)

        wavelength["wavelength_" + str(i)] = result_params['mu'].value
        wavelength["wavelength_error_" + str(i)] = result_params['mu'].stderr


for key in information:
    print(information[key])
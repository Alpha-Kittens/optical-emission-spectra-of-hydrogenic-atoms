# Imports
from data import data_loader
import lmfit
from fitters import fit
import matplotlib.pyplot as plt
import models
import noise_reduction

# File Path to peak data files
folder = 'data/DeuteriumScans/SuperFineScans'
alpha_L_fp  = folder + '/alpha5L'
alpha_H_fp  = folder + '/alpha4H'

beta_L_fp = folder + '/beta1L'
beta_H_fp = folder + '/beta2H'

gamma_L_fp = folder + '/gamma1L'
gamma_H_fp = folder + '/gamma2H'

delta_L_fp = folder + '/delta1L'
delta_H_fp = folder + '/delta2H'

# Wavelengths
alpha = {
    "1_fp" : alpha_L_fp,
    "2_fp" : alpha_H_fp,
}
beta = {
    "1_fp" : beta_L_fp,
    "2_fp" : beta_H_fp,
}
'''
gamma = {
    "1_fp" : gamma_L_fp,
    "2_fp" : gamma_H_fp,
}
delta = {
    "1_fp" : delta_L_fp,
    "2_fp" : delta_H_fp,
}
'''

#information = [alpha, beta, gamma, delta]
information = [alpha, beta]

# Run fits for each wavelength in information
for wavelength in information:
    for i in range(1, 3):
        # Read data
        data = data_loader.read_data(wavelength[str(i) + "_fp"])

        # Noise reduction and find Uncertainties
        #new_data, weights = noise_reduction.reduce_noise(data)


        # Voigt Model & Fit
        x_axis = data[:, 0]
        y_axis = data[:, 1]

        model = lmfit.models.VoigtModel()

        params = models.voigt_params(model, data)

        #result = fit(model, data, params, weights)
        result = fit(model, data, params)

        #lmfit.report_fit(result)

        result.plot()
        plt.show()

        wavelength[str(i) + "_result"] = result

        # Maximum Counts

    # Isotope shifts (deuterium wavelength - hydrogen wavelength)
    wavelength["shift"] = abs(wavelength["1_result"].best_values['center'] - wavelength["2_result"].best_values['center'])
    print(wavelength["shift"])


# Balmer Formula to find ratio of deuterium to proton mass


# Average and Standard Dev.
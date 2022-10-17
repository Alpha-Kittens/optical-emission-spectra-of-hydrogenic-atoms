# Imports
from data import data_loader
import lmfit
from fitters import fit
import matplotlib.pyplot as plt
import models_2 as models
import noise_reduction
import math

# File Path to peak data files
folder = 'data/DeuteriumScans/SuperFineScans'
alpha_L_fp  = folder + '/alpha5L'
alpha_H_fp  = folder + '/alpha4H'

beta_L_fp = folder + '/beta1L'
beta_H_fp = folder + '/beta2H'

gamma_L_fp = folder + '/gamma2L'
gamma_H_fp = folder + '/gamma1H'

delta_L_fp = folder + '/delta1L'
delta_H_fp = folder + '/delta2H'


#Constants
m_p = 1.67262192E-27 #kilograms
m_e = 9.1093837E-31 #kilograms
nf = 2
R_inf = 10973731.568160 #m^-1

# Wavelengths
alpha = {
    "1_fp" : alpha_L_fp,
    "2_fp" : alpha_H_fp,
    "ni"   : 3,
}
beta = {
    "1_fp" : beta_L_fp,
    "2_fp" : beta_H_fp,
    "ni"   : 4,
}
gamma = {
    "1_fp" : gamma_L_fp,
    "2_fp" : gamma_H_fp,
    "ni"   : 5,
}
'''
delta = {
    "1_fp" : delta_L_fp,
    "2_fp" : delta_H_fp,
    "ni"   : 6,
}
'''

#information = [alpha, beta, gamma, delta]
#information = [alpha, beta]
information = {
    'alpha' : alpha,
    'beta' : beta,
    'gamma' : gamma,
    #'delta' : delta,
}
damping_constants = {
    'alpha' : 1/10,
    'beta' : 1/3,
    'gamma' : 1/10,
    #'delta' : 1/10
}

# Run fits for each wavelength in information
for key in information:

    wavelength = information[key]

    for i in range(1, 3):
        # Read data
        data = data_loader.read_data(wavelength[str(i) + "_fp"])

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


        wavelength[str(i) + "_result"] = result

        # Maximum Counts Method
    
    # Calibration to convert the wavelengths

    # Isotope shifts (deuterium wavelength - hydrogen wavelength)
    wavelength["shift"] = abs(wavelength["1_result"].best_values['mu'] - wavelength["2_result"].best_values['mu'])

    # Uncertainty in Shift
    wavelength["shift_unc"] = math.sqrt((wavelength["1_result"].params['mu'].stderr)**2 + (wavelength["1_result"].params['mu'].stderr)**2)

    # Balmer Formula to find ratio of deuterium to proton mass
    wavelength["ratio"] = 1/(1 - R_inf*(m_p/m_e)*wavelength["shift"]*(10**(-10))*((1/(nf**2)) - (1/(wavelength["ni"]**2))))

    wavelength["ratio_unc"] = ((wavelength["ratio"])**2)*R_inf*(m_p/m_e)*((1/(nf**2)) - (1/(wavelength["ni"]**2))) * wavelength["shift_unc"]*(10**(-10))

# Mean
sum = 0
for key in information:
    wavelength = information[key]
    print("ratio: " + str(wavelength["ratio"]))
    sum += wavelength["ratio"]

mean = sum/len(information)
print(mean)

# Statistical Uncertainty
sum = 0
for key in information:
    wavelength = information[key]
    sum += (wavelength["ratio"] - mean)**2

std_dev = math.sqrt(sum/(len(information) - 1))

stat_unc = std_dev/math.sqrt(len(information))


# Systematic Uncertainty
sum = 0
for key in information:
    wavelength = information[key]
    sum += (wavelength["shift_unc"])**2
total_shift_unc = math.sqrt(sum)/math.sqrt(len(information))
print(total_shift_unc)


sum = 0
for key in information:
    wavelength = information[key]
    sum += (wavelength["ratio_unc"])**2
sys_unc = math.sqrt(sum)/math.sqrt(len(information))
print(sys_unc)



results = {
    "mean" : mean,
    "standard deviation" : std_dev,
    "stat uncertainty" : stat_unc,
    "sys uncertainty" : sys_unc
}

print(results)
# Imports
from data import data_loader
import lmfit
from fitters import fit

# File Path to peak data files
folder = 'data/DeuteriumScans/SuperFineScans'

alpha_L_fp  = folder + '/alpha5L'
alpha_H_fp  = folder + '/alpha4H'

beta_L_fp = folder + '/beta1L'
beta_H_fp = folder + '/beta2H'

# Read data
data = data_loader.read_data(beta_L_fp)


# Noise reduction and find Uncertainties


# Voigt Model & Fit
model = lmfit.models.VoigtModel()


result = fit(model, data)

lmfit.report_fit(result)

# Maximum Counts

# Isotope shifts (deuterium wavelength - hydrogen wavelength)


# Balmer Formula to find ratio of deuterium to proton mass


# Average and Standard Dev.
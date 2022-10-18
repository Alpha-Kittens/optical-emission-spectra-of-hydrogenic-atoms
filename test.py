# Imports
from data import data_loader
import lmfit
from data_processing import get_cutoff
from fitters import fit
import matplotlib.pyplot as plt
import models_2 as models
import noise_reduction
import math
from calibration import calibration
from calibration import calibration_error

'''
This file is created for running tests when needed
no worries about keeping it organized in any way
'''


folder = 'data/10.14.22/hydrogen'
alpha_H_fp  = folder + '/alpha-superfine-1'


#Before Calibration
data = data_loader.read_data(alpha_H_fp)
new_data, weights = noise_reduction.reduce_noise(data, plot=False) # damping constant default is 1/10

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title('uncalibrated')
#ax1.plot(data[:,0], data[:, 1])
#ax1.plot(new_data[:,0], new_data[:, 1])


model = lmfit.Model(models.voigt_with_shift)
params = models.voigt_shift_params(data)
result = fit(model, data, params, weights)        
print(lmfit.fit_report(result))
ax1.plot(data[:,0], data[:, 1], 'o')
ax1.plot(data[:,0], result.init_fit, '--', label='initial fit')
ax1.plot(data[:,0], result.best_fit, '-', label='best fit')


calibrated_data = data
for i in range(len(data[:, 0])):
    calibrated_data[i, 0] = calibration(data[i, 0])
ax2.set_title('calibrated')
new_calibrated_data, calibrated_weights = noise_reduction.reduce_noise(calibrated_data, plot=False)

#ax2.plot(calibrated_data[:,0], calibrated_data[:, 1])
#ax2.plot(new_calibrated_data[:,0], new_calibrated_data[:, 1])


model_calibrated = lmfit.Model(models.voigt_with_shift)
calibrated_params = models.voigt_shift_params(calibrated_data)
calibrated_result = fit(model_calibrated, calibrated_data, calibrated_params, calibrated_weights)        
print(lmfit.fit_report(calibrated_result))
ax2.plot(calibrated_data[:,0], calibrated_data[:, 1], 'o')
ax2.plot(calibrated_data[:,0], calibrated_result.init_fit, '--', label='initial fit')
ax2.plot(calibrated_data[:,0], calibrated_result.best_fit, '-', label='best fit')



plt.show()
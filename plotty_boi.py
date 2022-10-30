import matplotlib.pyplot as plt
from data.data_loader import read_data
from noise_reduction import reduce_noise
import numpy as np

#mercury = "data/day_1/"
#coarse = "calibration_scan_3_5350_5500"
true = 5460.74
#fine = "data/10.14.22/mercury/5460.74-fine"
fine = "data/final_data/mercury/main/5460_74"

data = read_data(fine)
noise_reduced, weights = reduce_noise(data)
plt.xlabel("Monochrometer reading (Ã)")
plt.ylabel("PMT counts per second")
plt.errorbar(data[:,0], data[:,1], yerr = 1/np.array(weights), label = "Monochrometer data", c = 'b', ls = 'none', marker = '.')
plt.xticks([5422, 5422.5, 5423, 5423.5, 5424])
#plt.axvline(x = 5422.5, ls = '--', c = 'orange', label = "Rough peak bounds")
#plt.axvline(x = 5424, ls = '--', c = 'orange')
#plt.axvline(x = true, ls = '--', c = 'r', label = "True wavelength: "+str(true)+" (A)")
plt.legend()
plt.show()
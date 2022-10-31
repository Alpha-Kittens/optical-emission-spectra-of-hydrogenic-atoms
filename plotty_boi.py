import matplotlib.pyplot as plt
from data.data_loader import read_data
from noise_reduction import reduce_noise
import numpy as np

mercury = "data/day_1/"
coarse = mercury + "calibration_scan_3_5350_5500"
true = 5460.74
#fine = "data/10.14.22/mercury/5460.74-fine"
#hfine = "data/final_data/deuterium/alphaH"
#dfine = "data/final_data/deuterium/alphaD"

#plt.title("Splitting of Sodium yellow doublet")
plt.xlabel("Monochromator reading (Ã)")
plt.ylabel("PMT counts per second")
#for fine in [hfine, dfine]:

data = read_data(coarse)
noise_reduced, weights = reduce_noise(data)

plt.errorbar(data[:,0], data[:,1], yerr = np.sqrt(data[:,0]), label = "Monochrometer data", c = 'b', marker = '.')
#plt.axvline(x = 5422.5, ls = '--', c = 'orange', label = "Rough peak bounds")
#plt.axvline(x = 5424, ls = '--', c = 'orange')
plt.axvline(x = true, ls = '--', c = 'r', label = "True wavelength: "+str(true)+" (A)")
plt.legend()
plt.show()
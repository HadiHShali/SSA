# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 11:28:44 2025

@author: GeodesyLab
"""

import numpy as np
rand = np.random.default_rng(42)

#To detect periodic signals in unevenly spaced observations, consider the following data:
t = 100 * rand.random(100)
y = np.sin(2 * np.pi * t) + 0.1 * rand.standard_normal(100)


#These are 100 noisy measurements taken at irregular times, with a frequency of 1 cycle per unit time.

#The Lomb-Scargle periodogram, evaluated at frequencies chosen automatically based on the input data, can be computed as follows using the LombScargle class
from astropy.timeseries import LombScargle
frequency, power = LombScargle(t, y).autopower(minimum_frequency=0.1,
                                                   maximum_frequency=1.9,
                                                   )
#Plotting the result with Matplotlib gives:
import matplotlib.pyplot as plt  
plt.figure()
plt.plot(frequency, power)

plt.figure()
plt.plot(t,y,'*b')
plt.show()

#The periodogram shows a clear spike at a frequency of 1 cycle per unit time, as we would expect from the data we constructed.



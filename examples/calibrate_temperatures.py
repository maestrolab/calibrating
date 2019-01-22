# -*- coding: utf-8 -*-
"""
Created on Jan 21 2019
@author: Pedro Leal
"""

import matplotlib.pyplot as plt

from calibration.temperature import processing_raw, fitting
from calibration.temperature.tangent import Tangent

optimizer = 'differential_evolution'
raw_data = processing_raw("../data/filtered_data_50MPa.txt")

for transformation in ['Austenite', 'Martensite']:
    f = Tangent(transformation, raw_data[transformation])
    f = fitting(f, optimizer)
    f.plotting()

plt.grid()
x, y, z = f.raw_data.T
plt.legend(loc="lower left")
plt.xlabel("Temperature (C)")
plt.ylabel("Strain (m/m)")
plt.show()

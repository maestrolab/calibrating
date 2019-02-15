# -*- coding: utf-8 -*-
"""
Created on Jan 21 2019
@author: Pedro Leal
"""

import matplotlib.pyplot as plt

from calibration.temperature import processing_raw, fitting_temperatures, fitting_transformations
from calibration.temperature.tangent import Tangent
from calibration.temperature.plotting import plot_tangents

optimizer = 'differential_evolution'
driven = 'temperature'
constant_stress = 200  # MPa
standard = True
filename = "../data/NiTiHf_UNT/filtered_data_"+str(constant_stress)+"MPa.txt"
raw_data = processing_raw(filename, driven, constant_stress)
# raw_data = processing_raw("../data/NiTi_flexinol/filtered_data_50MPa.txt")

colors = {'Austenite': ['--r', 'r'], 'Martensite': ['--b', 'b']}
surfaces = {}
for transformation in ['Austenite', 'Martensite']:
    surfaces[transformation] = Tangent(transformation,
                                       raw_data[transformation], standard)
    surfaces[transformation] = fitting_temperatures(surfaces[transformation],
                                                    optimizer)
    surfaces[transformation].plotting(transformation, colors[transformation])
# a, m = fitting_transformations(surfaces['Austenite'],
#                                surfaces['Martensite'], optimizer)

plt.grid()
plt.show()
# x, y, z = f.raw_data.T

# plot_tangents({'Austenite': a, 'Martensite': m}, constant_stress)

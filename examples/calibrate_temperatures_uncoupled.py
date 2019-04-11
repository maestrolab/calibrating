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

constant_stress = 100
standard = True
hot_to_cold = True
filter_rate = 40
filename = "../data/NiTiHf_Karaman/filtered_data_"+str(constant_stress)+"MPa.txt"

raw_data = processing_raw(filename, driven, constant_stress, filter_rate,
                          hot_to_cold)

colors = {'Austenite': ['--r', 'r'], 'Martensite': ['--b', 'b']}
surfaces = {}
for transformation in ['Austenite', 'Martensite']:
    surfaces[transformation] = Tangent(transformation,
                                       raw_data[transformation], standard)
    surfaces[transformation] = fitting_temperatures(surfaces[transformation],
                                                    optimizer)
    surfaces[transformation].plotting(transformation, colors[transformation])

plt.grid()
plt.show()

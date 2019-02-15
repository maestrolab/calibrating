# -*- coding: utf-8 -*-
"""
Created on Jan 21 2019
@author: Pedro Leal
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle

from calibration.temperature import processing_raw, fitting_temperatures, Transformation_Surface, fitting_transformation, fitting_transformations
from calibration.temperature.plotting import plot_tangents
from calibration.temperature.tangent import Tangent


# Inputs
optimizer = 'differential_evolution'
driven = 'temperature'
constant_stresses = [100, 200, 300]
standard = False
location_of_file = "../data/NiTiHf_UNT/filtered_data_"

# Calibrating tangents
calibrated_tangents = {'Austenite': [], 'Martensite': []}
for constant_stress in constant_stresses:
    filename = location_of_file + str(constant_stress) + "MPa.txt"
    raw_data = processing_raw(filename, driven, constant_stress)
    surfaces = {'Austenite': [], 'Martensite': []}
    for transformation in ['Austenite', 'Martensite']:
        surfaces[transformation] = Tangent(transformation,
                                           raw_data[transformation], standard)

    a, m = fitting_transformations(surfaces['Austenite'],
                                   surfaces['Martensite'], optimizer)
    calibrated_tangents['Austenite'].append(a)
    calibrated_tangents['Martensite'].append(m)

plot_tangents(calibrated_tangents, constant_stresses)

pickle.dump(calibrated_tangents,
            open('../data/NiTiHf_UNT/temperatures.p', 'wb'))

# Calibrating surfaces
calibrated_surfaces = {}
for trans in ['Austenite', 'Martensite']:
    f = Transformation_Surface(trans, calibrated_tangents[trans])
    f = fitting_transformation(f, optimizer)
    f.plotting()

plt.show()

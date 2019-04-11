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


# # Inputs
optimizer = 'differential_evolution'
driven = 'temperature'

# constant_stresses = [50, 100, 150, 172, 200]
constant_stresses = [100, 200, 300,400]
# constant_stresses = [50, 100, 200, 300, 400, 500, 600]
standard = True
hot_to_cold = True
filter_rate = 40
# location_of_file = "../data/NiTi_flexinol/filtered_data_"
# location_of_file = "../data/NiTiHf_UNT/filtered_data_"
# location_of_file = "../data/NiTiHf_Karaman/filtered_data_"
# location_of_file = "../data/Alex/MD_NI50_8TI49_2VF0/filtered_data_"
location_of_file = "../data/Alex/MD_NI50_8TI49_2VF0_V2/filtered_data_"

colors = {'Austenite': ['--r', 'r'], 'Martensite': ['--b', 'b']}


# Calibrating tangents
calibrated_tangents = {'Austenite': [], 'Martensite': []}
for constant_stress in constant_stresses:
    filename = location_of_file + str(constant_stress) + "MPa.txt"
    raw_data = processing_raw(filename, driven, constant_stress, filter_rate,
                              hot_to_cold)
    surfaces = {'Austenite': [], 'Martensite': []}
    for transformation in ['Austenite', 'Martensite']:
        surfaces[transformation] = Tangent(transformation,
                                           raw_data[transformation],
                                           constant_stress,
                                           standard)
        surfaces[transformation] = fitting_temperatures(surfaces[transformation],
                                                        optimizer)
        # surfaces[transformation].plotting(transformation,
        #                                   colors[transformation])
        calibrated_tangents[transformation].append(surfaces[transformation])
    plt.show()
# plt.grid()
# plt.show()

pickle.dump(calibrated_tangents,
            open('../data/NiTiHf_Karaman/temperatures.p', 'wb'))
# calibrated_tangents = pickle.load(open('../data/NiTiHf_Karaman/temperatures.p', 'rb'))
plot_tangents(calibrated_tangents, constant_stresses)


# Calibrating surfaces
calibrated_surfaces = {}
for trans in ['Austenite', 'Martensite']:
    f = Transformation_Surface(trans, calibrated_tangents[trans])
    f = fitting_transformation(f, optimizer)
    f.plotting()

plt.show()

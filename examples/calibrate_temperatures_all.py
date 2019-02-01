# -*- coding: utf-8 -*-
"""
Created on Jan 21 2019
@author: Pedro Leal
"""

import matplotlib.pyplot as plt

from calibration.temperature import processing_raw, fitting
from calibration.temperature.tangent import Tangent
from calibration.temperature.combine_plots import combine_plots


# Inputs
optimizer = 'differential_evolution'
driven = 'temperature'
constant_stresses = [100, 200, 300]
location_of_file = "../data/NiTiHf_UNT/filtered_data_" # filename format: "filtered_data_[constant_stress].MPa.txt"

combine_plots(driven,constant_stresses,location_of_file, optimizer)

# -*- coding: utf-8 -*-
"""
Created on Jan 21 2019
@author: Pedro Leal
"""

import matplotlib.pyplot as plt
import numpy as np
import pickle

from calibration.temperature import processing_raw, fitting_temperatures, Transformation_Surface, fitting_transformation
from calibration.temperature.plotting import plot_tangents
from calibration.temperature.tangent import Tangent

optimizer = 'differential_evolution'
calibrated_tangents = pickle.load(open('../data/NiTiHf_UNT/temperatures.p', 'rb'))

# Calibrating surfaces
calibrated_surfaces = {}
for trans in ['Austenite', 'Martensite']:
    f = Transformation_Surface(trans, calibrated_tangents[trans])
    f = fitting_transformation(f, optimizer)
    f.plotting()
plt.show()

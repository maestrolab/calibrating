# -*- coding: utf-8 -*-
"""
Created on Jan 21 2019
@author: Pedro Leal
"""

import matplotlib.pyplot as plt
import numpy as np


class Tangent():
    """Class for tangent lines
    - transformation: Austenite or Martensite
    - raw_data: numpy.array with data for (temperatture, strain, stress)"""

    def __init__(self, transformation, raw_data):
        if transformation == 'Austenite':
            self.T_1, self.T_4 = raw_data[0, 0], raw_data[-1, 0]
        elif transformation == 'Martensite':
            self.T_4, self.T_1 = raw_data[0, 0], raw_data[-1, 0]
        self.raw_data = raw_data.copy()
        self.transformation = transformation

        # Default values for bounds and x0
        self.bounds = [(30, 120), (10, 60)] + \
            4*[(min(self.raw_data[:, 1]), max(self.raw_data[:, 1])), ]
        self.x0 = [(x[0]+x[1])/2. for x in self.bounds]

    def lines(self, T_i):
        """Calculates tangent line function for a value T_i
        - T_i: float to calculate estimate value of strain"""
        [T, strain] = self.props.T
        for j in range(3):
            if T_i - T[j+1] < 1e-5:
                diff = (strain[j+1] - strain[j])/(T[j+1] - T[j])
                return diff*(T_i-T[j]) + strain[j]

    def update(self, x):
        """Update properties based on design vector from optimizer
        - x: [T_1, T_2, strain_1, strain_2, strain_3, strain_4]"""
        T = [self.T_1, x[0], x[0] + x[1], self.T_4]
        strain = x[2:6]
        self.props = np.vstack([T, strain]).T

    def error(self, x):
        """Calculate root mean squared"""
        self.update(x)
        f = np.array([self.lines(T_i) for T_i in self.raw_data[:, 0]])
        strain = self.raw_data[:, 1]
        root_mean = np.sqrt(np.sum((f-strain)**2)/len(strain))
        return root_mean

    def plotting(self):
        """Plot raw and tangent lines"""
        T, epsilon, sigma = self.raw_data.T

        f = []
        for i in range(len(T)):
            f.append(self.lines(T[i]))

        plt.plot(T, f, 'b', label="Tangents")
        plt.plot(T, epsilon, 'g', label="Experimental data")

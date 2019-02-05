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
    - raw_data: numpy.array with data for (temperature, strain, stress)"""

    def __init__(self, transformation, raw_data, standard=False):
        if transformation == 'Austenite':
            self.T_1, self.T_4 = raw_data[0, 0], raw_data[-1, 0]
        elif transformation == 'Martensite':
            self.T_4, self.T_1 = raw_data[0, 0], raw_data[-1, 0]
        self.raw_data = raw_data.copy()
        self.transformation = transformation
        self.standard = standard

        # Default values for bounds and x0
        if not standard:
            n_strains = 4
        else:
            n_strains = 3

        self.bounds = 2*[(min(self.raw_data[:, 0]), max(self.raw_data[:, 0]))] + \
            n_strains*[(min(self.raw_data[:, 1]), max(self.raw_data[:, 1])), ]
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
        - x: [T_2, T_3, strain_1, strain_2, strain_3, strain_4]"""
        if not self.standard:
            T = [self.T_1, x[0], x[0] + x[1], self.T_4]
            strain = x[2:6]
        else:
            T = [self.T_1, x[0], x[0] + x[1], self.T_4]
            T_50 = (T[1] + T[2])/2.
            strain_50 = np.interp(T_50, self.raw_data[:,0], self.raw_data[:,1])
            strain_2 = x[3]
            #a = (strain_2-strain_50)/(T[1] - T_50)
            #b = strain_50 - a*T_50
            a = (strain_50-strain_2)/(T_50-T[1])
            b = strain_50 - a*T_50
            #strain_3 = a+b*T[2]
            strain_3 = a*T[2]+b
            strain = [x[2], strain_2, strain_3, x[-1]]
            plt.figure()
            plt.plot(self.raw_data[:,0], self.raw_data[:,1])
            plt.scatter(T_50, strain_50, color='m')
            plt.scatter(T, strain)
            plt.show()
        self.props = np.vstack([T, strain]).T

    def error(self, x):
        """Calculate root mean squared"""
        self.update(x)
        f = np.array([self.lines(T_i) for T_i in self.raw_data[:, 0]])
        strain = self.raw_data[:, 1]
        root_mean = np.sqrt(np.sum((f-strain)**2)/len(strain))
        return root_mean

    def plotting(self, label=None, color=['r','b']):
        """Plot raw and tangent lines"""
        T, epsilon, sigma = self.raw_data.T

        n = 200
        T_tangent = np.linspace(min(T), max(T), n)
        f = []
        for i in range(n):
            f.append(self.lines(T_tangent[i]))

        plt.plot(T_tangent, f, color[0], label="Tangent (" + label + ")")
        plt.plot(T, epsilon, color[1], label="Experimental data (" + label + ")")

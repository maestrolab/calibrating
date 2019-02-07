"""
Created on Jan 21 2019
@author: Pedro Leal
"""

from scipy.optimize import differential_evolution, minimize
from calibration.filehandling import output_reader
import numpy as np
import matplotlib.pyplot as plt


class Transformation_Surface():
    """Class for tangent lines
    - transformation: Austenite or Martensite
    - raw_data: numpy.array with data for (temperature, strain, stress)"""

    def __init__(self, transformation, f_list):
        self.raw_start = []
        self.raw_finish = []
        self.stress = []
        for f in f_list:
            self.raw_start.append(f.start)
            self.raw_finish.append(f.finish)
            self.stress.append(f.stress)
        self.raw_start = np.vstack([self.raw_start, self.stress]).T
        self.raw_finish = np.vstack([self.raw_finish, self.stress]).T

        self.transformation = transformation
        self.bounds = np.array([(0., 200.), (0., 150.), (0.1, 30.)])
        self.x0 = [(x[0]+x[1])/2. for x in self.bounds]

    def error(self, x):
        self.update(x)

        T_raw, sigma_raw = self.raw_start.T
        T = self.lines(sigma_raw)
        raw = [self.raw_start, self.raw_finish]
        rmse = 0
        for i in range(2):
            T_raw, sigma_raw = raw[i].T
            T_i = T[i].T
            rmse += np.sqrt(np.sum((T_i-T_raw)**2)/len(sigma_raw))
        return rmse

    def update(self, x):
        if self.transformation == 'Austenite':
            self.start = x[0]
            self.finish = x[0] + x[1]
            self.slope = x[2]
        else:
            self.start = x[0] + x[1]
            self.finish = x[0]
            self.slope = x[2]

    def lines(self, sigma):
        T_s = sigma/self.slope + self.start
        T_f = sigma/self.slope + self.finish
        return T_s, T_f

    def plotting(self, surfaces=True):
        if self.transformation == 'Austenite':
            color = 'r'
        else:
            color = 'b'

        T, sigma = self.raw_start.T
        plt.scatter(T, sigma, c=color, marker='o',
                    label=self.transformation + " Start")
        T, sigma = self.raw_finish.T
        plt.scatter(T, sigma, c=color, marker='s',
                    label=self.transformation + " Finish")

        if surfaces:
            sigma = np.array([0]+list(sigma))
            T_s, T_f = self.lines(sigma)
            plt.plot(T_s, sigma, color,
                     label=self.transformation + " Start")
            plt.plot(T_f, sigma, color,
                     label=self.transformation + " Finish")
        plt.legend(loc="lower left", prop={'size': 8})
        plt.xlabel("Temperature (C)")
        plt.ylabel("Stress (N/m^2)")


def fitting_temperatures(f, optimizer='differential_evolution'):
    """Optimize properties for class f to represent raw data.
       - f: any class with attributes .error, .x0, and .bound
       - optimizer: BFGS (gradient) or differential_evolution"""

    print('Fitting ' + f.transformation)

    if optimizer == 'BFGS':
        result = minimize(f.error, f.x0, method='BFGS')
    elif optimizer == 'differential_evolution':
        result = differential_evolution(f.error, f.bounds, popsize=100,
                                        maxiter=100)

    f.update(result.x)
    if f.transformation == 'Austenite':
        print('As=', f.start, 'Af=', f.finish)
    elif f.transformation == 'Martensite':
        print('Ms=', f.start, 'Mf=', f.finish)
    return(f)


def fitting_transformation(f, optimizer='differential_evolution'):
    """Optimize properties for class f to represent raw data.
       - f: any class with attributes .error, .x0, and .bound
       - optimizer: BFGS (gradient) or differential_evolution"""

    print('Fitting ' + f.transformation)

    if optimizer == 'BFGS':
        result = minimize(f.error, f.x0, method='BFGS')
    elif optimizer == 'differential_evolution':
        result = differential_evolution(f.error, f.bounds, popsize=100,
                                        maxiter=100)

    f.update(result.x)
    if f.transformation == 'Austenite':
        print('As=', f.start, 'Af=', f.finish, 'C=', f.slope)
    elif f.transformation == 'Martensite':
        print('Ms=', f.start, 'Mf=', f.finish, 'C=', f.slope)
    return(f)


def processing_raw(filename, driven='temperature', constant_stress=None):
    """Convert .txt file to a numpy array (temperature, strain, sigma) for
       Austenite and Martensite.
       - filename: string for file to process"""

    raw_data = output_reader(filename)
    try:
        strain = raw_data['Strain']
        temperature = raw_data['Temperature']
        try:
            stress = raw_data['Stress']
        except KeyError:
            stress = np.ones(len(strain))*constant_stress
    except KeyError:
        raise "Formating error: should have 'Temperature' and 'Strain' columns"

    if driven == 'temperature':
        i = temperature.index(max(temperature))
    else:
        raise NotImplementedError

    austenite = np.vstack([temperature[:i+1], strain[:i+1], stress[:i+1]]).T
    martensite = np.vstack([temperature[i:], strain[i:], stress[i:]]).T
    return({'Austenite': austenite, 'Martensite': martensite})

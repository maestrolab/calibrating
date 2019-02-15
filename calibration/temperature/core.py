"""
Created on Jan 21 2019
@author: Pedro Leal
"""

from scipy.optimize import differential_evolution, minimize, bisect
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


def fitting_transformations(a, m, optimizer='differential_evolution'):
    """Optimize properties for class f to represent raw data.
       - f: any class with attributes .error, .x0, and .bound
       - optimizer: BFGS (gradient) or differential_evolution"""

    def error(x):
        """
        Originally (16 DOF):
        x = [Ta_2, Ta_3, Ea_1, Ea_2, Ea_3, Ea_4,
             Tm_2, Tm_3, Em_1, Em_2, Em_3, Em_4]
        Now:
        x = [Ta_2, Ta_3, Ea_1, Ea_2, Ea_4, Em_2]"""
        x_a, x_m = update(x, output=True)

        return(a.error(x_a) + m.error(x_m))

    def update(x, output=False):
        """
        Originally (16 DOF):
        x = [Ta_2, Ta_3, Ea_1, Ea_2, Ea_3, Ea_4,
             Tm_2, Tm_3, Em_1, Em_2, Em_3, Em_4]
        Now:
        x = [Ta_2, Ta_3, Ea_1, Ea_2, Ea_4, Em_2]"""
        x_a = x[:-1]

        Ta_2, Ta_3_delta, Ea_1, Ea_2, Ea_4, Em_2 = x

        Ta_3 = Ta_2 + Ta_3_delta
        Tm_1 = m.T_1
        Tm_4 = m.T_4

        Em_1 = Ea_1
        Em_4 = Ea_4
        Ea_3 = a._strain_3(Ta_2, Ta_3, Ea_2)
        Tm_2 = Tm_1 + ((Ta_2-Tm_1)*(Em_2-Em_1))/(Ea_2-Ea_1)

        def aa(Tm_3):
            return m._strain_3(Tm_2, Tm_3, Em_2)

        def bb(Tm_3):
            return Em_4 - ((Tm_4-Tm_3)*(Ea_4-Ea_3))/(Tm_4-Ta_3)

        def cc(Tm_3):
            return aa(Tm_3)-bb(Tm_3)

        # plt.figure()
        # T = np.linspace(0, 300, 200)
        # plt.plot(T, cc(T))
        # plt.show()
        res = minimize(cc, Ta_3)
        Tm_3 = res.x

        x_m = [Tm_2, Tm_3-Tm_2, Em_1, Em_2, Em_4]
        a.update(x_a)
        m.update(x_m)
        # plt.figure()
        # plt.plot([a.T_1, Ta_2, Ta_3, a.T_4], [Ea_1, Ea_2, Ea_3, Ea_4], 'r')
        # x, y, z = a.raw_data.T
        # plt.plot(x, y, '--r')
        # plt.scatter([a.T_1, Ta_2, Ta_3, a.T_4], [Ea_1, Ea_2, Ea_3, Ea_4], c='r')
        # plt.plot([m.T_1, Tm_2, Tm_3, m.T_4], [Em_1, Em_2, Em_3, Em_4], 'b')
        # x, y, z = m.raw_data.T
        # plt.plot(x, y, '--b')
        # plt.scatter([m.T_1, Tm_2, Tm_3, m.T_4], [Em_1, Em_2, Em_3, Em_4], c='b')
        # plt.show()
        if output:
            return x_a, x_m

    bounds = a.bounds + m.bounds[-1:]
    x0 = a.x0 + m.x0[-1:]  # we want this to be 6
    print('bounds', bounds)
    print('x0', x0)
    # [Ta_2, Ta_3, Ea_1, Ea_2, Ea_4, Em_2]
    error([150, 50, 4.2, 4.1, 0, 4.15])

    if optimizer == 'BFGS':
        result = minimize(error, x0, method='BFGS')
    elif optimizer == 'differential_evolution':
        result = differential_evolution(error, bounds, popsize=15,
                                        maxiter=100)

    update(result.x, output=True)
    print('As=', a.start, 'Af=', a.finish)
    print('Ms=', m.start, 'Mf=', m.finish)
    return(a, m)


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

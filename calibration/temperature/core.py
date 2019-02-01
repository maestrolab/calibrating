"""
Created on Jan 21 2019
@author: Pedro Leal
"""

from scipy.optimize import differential_evolution, minimize
from calibration.filehandling import output_reader
import numpy as np


def fitting(f, optimizer='differential_evolution'):
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
        f.start, f.finish = f.props[1, 0], f.props[-2, 0]
        print('As=', f.start, 'Af=', f.finish)
    elif f.transformation == 'Martensite':
        f.start, f.finish = f.props[-2, 0], f.props[1, 0]
        print('Ms=', f.start, 'Mf=', f.finish)
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

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
        print('As=', f.props[0, 0], 'Af=', f.props[-1, 0])
    elif f.transformation == 'Martensite':
        print('Ms=', f.props[-1, 0], 'Mf=', f.props[0, 0])
    return(f)


def processing_raw(filename, driven='temperature'):
    """Convert .txt file to a numpy array (temperature, strain, sigma) for
       Austenite and Martensite.
       - filename: string for file to process"""
    raw_data = output_reader(filename, header=['Temperature', "Strain",
                                               "Stress"],)
    temperature = raw_data['Temperature']
    strain = raw_data['Strain']
    stress = raw_data['Stress']

    if driven == 'temperature':
        i = temperature.index(max(temperature))
    else:
        raise NotImplementedError

    austenite = np.vstack([temperature[:i+1], strain[:i+1], stress[:i+1]]).T
    martensite = np.vstack([temperature[i:], strain[i:], stress[i:]]).T
    return({'Austenite': austenite, 'Martensite': martensite})

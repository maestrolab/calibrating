"""
Created on Jan 21 2019
@author: Pedro Leal
"""

from scipy.optimize import differential_evolution, minimize
from calibration.filehandling import output_reader
import numpy as np
import matplotlib.pyplot as plt

def fitting(f, optimizer='differential_evolution'):
    """Optimize properties for class f to represent raw data.
       - f: any class with attributes .error, .x0, and .bound
       - optimizer: BFGS (gradient) or differential_evolution"""
       #add a flag for true/false standard vs. non standard

    print('Fitting ' + f.transformation)

    if optimizer == 'BFGS':
        result = minimize(f.error, f.x0, method='BFGS')
    elif optimizer == 'differential_evolution':
        result = differential_evolution(f.error, f.bounds, popsize=100,
                                        maxiter=100)

    f.update(result.x)
    if f.transformation == 'Austenite':
        f.start, f.finish = f.props[1, 0], f.props[-2, 0]
        # A50 = f.finish + ((f.start-f.finish)/2)
        #
        # AstartTempmin = raw_data[:,0][raw_data[:,0] > f.start].min()
        # AfinishTempmax = raw_data[:,0][raw_data[:,0] < f.finish].max()
        #
        # AstartTempminloc = raw_data[:,0].tolist().index(AstartTempmin) #temperature
        # AfinishTempmaxloc = raw_data[:,0].tolist().index(AfinishTempmax) #temperature
        #
        # AstartStrain = raw_data[:,1][AstartTempminloc]
        # AfinishStrain = raw_data[:,1][AfinishTempmaxloc]
        #
        # A50min = raw_data[:,0][raw_data[:,0] > A50].min() #temperature
        # A50minloc = raw_data[:,0].tolist().index(A50min) #location of Temp and strain
        # A50minStrain = raw_data[:,1][A50minloc] #strain
        #
        # A50max = raw_data[:,0][raw_data[:,0] < A50].max() #temperature
        # A50maxloc = raw_data[:,0].tolist().index(A50min) #location of Temp and strain
        # A50maxStrain = raw_data[:,1][A50maxloc] #strain
        #
        # x1 = [A50min, A50max] #temperature
        # y1 = [A50minStrain, A50maxStrain] #strain
        # xvals = np.linspace(f.start, f.finish) #temperature
        #
        # A50_strain = np.interp(A50,raw_data[:,0],raw_data[:,1]) #temp, strain
        #
        # #linear interpolate strain
        print('As=', f.start, 'Af=', f.finish)
    elif f.transformation == 'Martensite':
        f.start, f.finish = f.props[-2, 0], f.props[1, 0]
        # M50 = (f.start + ((f.finish-f.start)/2))
        # M50min = raw_data[:,0][raw_data[:,0] > M50].min()
        # M50minloc = raw_data[:,0].tolist().index(M50min) #location of Temp and strain
        # M50minStrain = raw_data[:,1][M50minloc]
        #
        # M50max = raw_data[:,0][raw_data[:,0] < M50].max()
        # M50maxloc = raw_data[:,0].tolist().index(M50max) #location of Temp and strain
        # M50maxStrain = raw_data[:,1][M50maxloc]

        print('Ms=', f.start,  'Mf=', f.finish)
    return(f)

    # f.update(result.x)
    # if f.transformation == 'Austenite':
    #     f.start, f.mid, f.finish = f.props[1, 0], f.props[0,1], f.props[-2, 0]
    #     print('As=', f.start, 'As50=', f.mid, 'Af=', f.finish)
    # elif f.transformation == 'Martensite':
    #     f.start, f.mid, f.finish = f.props[-2, 0], f.props[0,1], f.props[1, 0]
    #     print('Ms=', f.start, 'Ms50=', f.mid, 'Mf=', f.finish)
    # return(f)

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

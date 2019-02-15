import matplotlib.pyplot as plt
import numpy as np


def plot_tangents(calibrated_tangents, constant_stresses):
    try:
        for i in range(len(constant_stresses)):
            plt.subplot(2, 2, i+1)
            plt.title(str(constant_stresses[i]) + ' MPa',
                      fontdict={'fontsize': 8, 'fontweight': 'medium'})
            for transformation in ['Austenite', 'Martensite']:
                calibrated_tangents[transformation][i].plotting()
            plt.legend(loc="lower left", prop={'size': 8})
            plt.xlabel("Temperature (C)")
            plt.ylabel("Strain (m/m)")
    except TypeError:
        plt.figure()
        for transformation in ['Austenite', 'Martensite']:
            calibrated_tangents[transformation].plotting()
        plt.legend(loc="lower left", prop={'size': 8})
        plt.xlabel("Temperature (C)")
        plt.ylabel("Strain (m/m)")
    plt.show()

import matplotlib.pyplot as plt
import numpy as np

from calibration.temperature import processing_raw, fitting
from calibration.temperature.tangent import Tangent
from calibration.temperature.plotting import plotstraintemp, plotstresstemp
from calibration.temperature.regression_fit import regfit

def combine_plots(driven,constant_stresses,location_of_file, optimizer):
  colors = {'Austenite': ['--r','r'],
            'Martensite': ['--b','b']}

  # Calibrating temperatures

  calibrated = {'Austenite': {'start':[], 'finish':[]},
                'Martensite': {'start':[], 'finish':[]}}

  fig = plt.figure()
  plt.grid()
  plotstraintemp(driven,constant_stresses,location_of_file, optimizer, colors, calibrated)

  plt.subplot(2,2,4)
  plotstresstemp(calibrated, constant_stresses)

  plt.tight_layout()
  # regfit(calibrated['Martensite']['start'], calibrated['Martensite']['finish'], calibrated['Austenite']['start'], calibrated['Austenite']['finish'], constant_stresses)

  plt.show()

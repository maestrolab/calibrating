import matplotlib.pyplot as plt

from calibration.temperature import processing_raw, fitting
from calibration.temperature.tangent import Tangent


def combine_plots(driven,constant_stresses,location_of_file, optimizer):
  colors = {'Austenite': ['--r','r'],
            'Martensite': ['--b','b']}

  # Calibrating temperatures

  calibrated = {'Austenite': {'start':[], 'finish':[]},
                'Martensite': {'start':[], 'finish':[]}}
  n = 1
  fig = plt.figure()

  plt.grid()
  for constant_stress in constant_stresses:
      filename = str(location_of_file) + str(constant_stress) + "MPa.txt"
      raw_data = processing_raw(filename, driven, constant_stress)

      plt.subplot(2,2,n)
      for transformation in ['Austenite', 'Martensite']:
          f = Tangent(transformation, raw_data[transformation])
          f = fitting(f, optimizer)

          plt.title('Temperature vs. Strain ' + str(constant_stress), fontdict={'fontsize': 8, 'fontweight': 'medium'})
          f.plotting(str(constant_stress)+' MPa', color=colors[transformation])


          calibrated[transformation]['start'].append(f.start)
          calibrated[transformation]['finish'].append(f.finish)
      plt.legend(loc="lower left", prop={'size': 8})
      plt.xlabel("Temperature (C)")
      plt.ylabel("Strain (m/m)")

      n = n+1

  # curve fit with regression lines
  #make everything a function!
  plt.subplot(2,2,4)
  plt.title('Temperature vs. Stress', fontdict={'fontsize': 8, 'fontweight': 'medium'})
  plt.scatter(calibrated['Austenite']['start'], constant_stresses, c='r',marker='o', label="Austenite Start")
  plt.scatter(calibrated['Austenite']['finish'], constant_stresses, c='r', marker = 's', label="Austenite Finish")
  plt.scatter(calibrated['Martensite']['start'], constant_stresses, c='b', marker = 'o', label="Martensite Start")
  plt.scatter(calibrated['Martensite']['finish'], constant_stresses, c='b', marker = 's', label="Martensite Finish")
  plt.legend(loc="lower left", prop={'size': 8})
  plt.xlabel("Temperature (C)")
  plt.ylabel("Stress (N/m^2)")

  plt.tight_layout()
  plt.show()

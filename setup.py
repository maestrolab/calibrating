"""Calibration: A tool to calibrate SMA properties."""

from setuptools import setup, find_packages

setup(name='calibration',
      version='0.0',
      description='A tool to calibrate SMA properties.',
      url='NA',
      author='leal26',
      author_email='leal26@tamu.edu',
      license='MIT',
      packages=['calibration', 'calibration.full', 'calibration.temperature'],
      zip_safe=False
      )

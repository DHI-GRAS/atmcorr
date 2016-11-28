#!/usr/bin/env python

from distutils.core import setup

setup(name='atmospheric_correction',
      version='0.1',
      description='atmospheric_correction toolkit',
      author='Jonas Solvsteen',
      author_email='josl@dhi-gras.com',
      url='https://www.dhi-gras.com',
      packages=['atmospheric_correction'],
      package_dir={'atmospheric_correction': ''},
      package_data={'atmospheric_correction': ['dependency\sensorResponseCurves\*.txt', 'dependency\sixsV1.1']}
      )

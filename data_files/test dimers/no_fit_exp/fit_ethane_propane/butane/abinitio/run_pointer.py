#!/usr/bin/env python
from pointer import FitFFParameters as Fit

settings = {}
settings['mon1'] = 'butane'
settings['mon2'] = 'butane'
settings['ofile_prefix'] = 'output/'
settings['ofile_suffix'] = '_unconstrained'
settings['fit_bii'] = False

settings_files = ['mastiff']

Fit(settings_files,**settings)


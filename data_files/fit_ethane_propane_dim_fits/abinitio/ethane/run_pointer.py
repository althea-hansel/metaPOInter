#!/usr/bin/env python
from pointer import FitFFParameters as Fit

settings = {}
settings['mon1'] = 'ethane'
settings['mon2'] = 'ethane'
settings['ofile_prefix'] = 'output/'
settings['ofile_suffix'] = '_unconstrained'
settings['fit_bii'] = False

settings_files = ['mastiff']

Fit(settings_files,**settings)


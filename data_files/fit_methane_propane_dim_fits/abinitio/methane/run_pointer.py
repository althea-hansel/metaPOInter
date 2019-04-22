#!/usr/bin/env python
from pointer import FitFFParameters as Fit

settings = {}
settings['mon1'] = 'methane'
settings['mon2'] = 'methane'
settings['ofile_prefix'] = 'output/'
settings['ofile_suffix'] = '_unconstrained'
settings['fit_bii'] = True

settings_files = ['mastiff']

Fit(settings_files,**settings)


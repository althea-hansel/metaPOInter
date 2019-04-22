#!/usr/bin/env python
from pointer import FitFFParameters as Fit

settings = {}
settings['mon1'] = 'propane'
settings['mon2'] = 'propane'
settings['ofile_prefix'] = 'output/fit_exp_'
settings['ofile_suffix'] = '_unconstrained'
settings['fit_bii'] = True

settings_files = ['mastiff']

Fit(settings_files,**settings)


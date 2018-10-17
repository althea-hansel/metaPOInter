#!/usr/bin/env python
from fit_ff_parameters import FitFFParameters as Fit

settings = {}
settings['mon1'] = 'ethane'
settings['mon2'] = 'ethane'
settings['ofile_prefix'] = 'output/fit_exp_'
settings['ofile_suffix'] = '_unconstrained'

settings_files = ['mastiff']

Fit(settings_files,**settings)


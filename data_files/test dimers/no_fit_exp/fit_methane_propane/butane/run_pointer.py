#!/usr/bin/env python
from pointer import FitFFParameters as Fit

settings = {}
settings['mon1'] = 'butane'
settings['mon2'] = 'butane'
settings['constraint_files'] = ['optimized_params.constraints']
settings['constrained_atomtypes'] = ['C', 'H']

settings_files = ['default']


Fit(settings_files,**settings)


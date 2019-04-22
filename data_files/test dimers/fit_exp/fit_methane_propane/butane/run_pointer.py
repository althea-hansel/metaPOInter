#!/usr/bin/env python
from pointer import FitFFParameters as Fit

settings = {}
settings['constraint_files'] = ['optimized_params.constraints']
settings['constrained_atomtypes'] = ['C', 'H']

settings_files = ['mastiff', 'config']


Fit(settings_files,**settings)


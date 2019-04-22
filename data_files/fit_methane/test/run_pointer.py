#!/usr/bin/env python
from pointer import FitFFParameters as Fit

settings = {}
settings['mon1'] = 'methane'
settings['mon2'] = 'methane'
settings['constraint_files'] = ['input/methane_methane_unconstrained.constraints']
settings['constrained_atomtypes'] = ['C', 'H']

settings_files = ['mastiff']


Fit(settings_files,**settings)


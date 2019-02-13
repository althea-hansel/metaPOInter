#######################################
# Known Bugs:

#######################################
"""
This Parameter Generalizer is designed to generate general force field parameters from libraries of ab initio MASTIFF force fields.

This program requires the MASTIFF fitting package POInter stored as a local module named pointer in python path.

This program was designed to be run using the Anaconda python distribution based on python 2.7

Required input directories and files (in same working directory as this script):
################################################################################
dir generalizerInput (or other specified in settings file):
	.sapt files for all ab initio FF
	dir input:
		all input files required by POInter

dir abInitioConstraints (or other specified in settings file):
	.constraints files with parameters for all ab initio FF

settings.py: see example

Output directories and files:
#############################
optimized_params.constraints

global_exponents.exp

dir dim_fits (or other specified in settings file):
	dir (for each dimer):
		coeffs.out
		dhf.dat (all .dat files contain DFT-SAPT energy vs FF energy for each dimer pair)
		dispersion.dat
		electrostatics.dat
		exchange.dat
		induction.dat
		total_energy.dat

Program Flow Description
#####################

The basic purpose of this program is to take ab initio parameters and fit them to multiple dimers. The program starts by reading in the parameters from the ab initio FF and storing them in the objects created below. Currently, the average of all ab initio parameters for a given atomtype is used as the starting point for optimization. The parameters are optimized for the minimal least-squared-error summed across all dimers for a given atomtype. 

Written by Althea Hansel, Harvey Mudd College '19 as part of Undergraduate Senior Thesis

"""
#TODO: where are "optimized" C params in output from?

# Standard packages
import numpy as np
import sys
import os
from scipy.optimize import minimize
import json
import glob
import math
from sympy.utilities.iterables import flatten

# Local modules
from pointer import fit_ff_parameters as Pointer
import settings as generalizer_settings

##################
# Global variables
##################

scriptDir = os.getcwd()

#variable used so store imported JSON dictionaries when working with library of ab initio FF's
inputJsons = {}

#list of different types of parameters we care about
parameterTypes = ['A', 'B', 'C', 'aniso']

#names of energy components
energyComponents = ['Exchange', 'Electrostatics', 'Induction', 'Dhf', 'Dispersion']

#names of disperision (C) components
dispersionTerms = ['C6', 'C8', 'C10', 'C12']

#use to store the parameters we will give to the minimize function as initial guesses when importing directly from a .params file
initialMinimizeValues = {}

#stores intial average parameter values
initialAverageValues = {}

#input list for minimize function
initialValuesList = []

#list of all atomtypes from all imported dimers
allAtomtypes = []

#stores atomtypes held by each dimers
dimer_atomtypes = {}

#list of aniso atomtypes
anisoAtomtypes = []

#contains all imported dimecular parameters
dimers = {}

#store a list of spherical harmonics for different atomtypes
atomtypeSphericalHarmonics = {}

#track what component is being fit
fitting_component = 0
# components are, in order, exchange (0),
            #                           electrostatics (1),
            #                           induction (2),
            #                           dhf (3),
            #                           dispersion (4),
            #                           residual error (5),
            #                           total energy (6)

#track how many parameters each atomtype has
atomtype_nparams = {}

#stores POInter objects for each dimer
dimer_POInter_objects = {}

#stores final parameters by component
optimized_params = {}

#store qm energies for each dimer
qm_energies = {}

#store ff energy for each dimer
ff_energies = {}

#stores initial B parameters for each atomtype
atomtype_B_params = {}

#stores C parameters for each atomtype
atomtype_C_params = {}

#stores drude charge for each atomtype
atomtype_drude_charges = {}

#flag for if drude calculation has already been performed for a given dimer
drude_flags = {}

#collects total lsq_error for each component
component_lsq_error = {component : 0.0 for component in energyComponents}

#flag for if isotropic dispersion should be scaled
scale_iso_disp = generalizer_settings.scale_iso_disp

####################################################################################################
# Class for holding the parameters from a single ab initio FF imported from a JSON .constraints file
####################################################################################################

class dimecularParameters:
	#creates object to hold parameters for an individual dimecular ab initio FF

	def __init__(self, jsonInput, dimer):
		self.jsonInput = jsonInput
		self.dimerName = dimer
		self.anisoFlag = False
		self.atomtypes = {}
		self.atomtypesList = []
		self.populateAtomtypes()

	def __repr__(self):
		return str(self.atomtypes)

	def __getitem__(self, atomtype):
		return self.atomtypes[atomtype]

	def setAnisoFlag(self, atomtype):
		if atomtype["sph_harm"] != []:
			self.anisoFlag = True
			return True
		else:
			return False

	def populateAtomtypes(self):
		"""takes input JSON file and makes dictionary of parameters by atomtype. Dictionary is stored in self.atomtypes instance variable

		Parameters
		----------
		none

		Returns 
		-------
		none

		"""
		self.atomtypesList = list(self.jsonInput.keys())
		#add list of atomtypes to entry in dimer_atomtypes and sort
		dimer_atomtypes[self.dimerName] = [atomtype for atomtype in self.atomtypesList]
		dimer_atomtypes[self.dimerName].sort()
		#populate dimer with atomtype objects
		for atomtype in self.atomtypesList:
			self.atomtypes[str(atomtype)] = self.Atomtype(self.jsonInput[atomtype], self.setAnisoFlag(self.jsonInput[atomtype]))
			#store spherical harmonics and drude charges for later access
			global atomtypeSphericalHarmonics
			atomtypeSphericalHarmonics[str(atomtype)] = self.atomtypes[str(atomtype)].getSphericalHarmonics()
			atomtype_drude_charges[str(atomtype)] = self.atomtypes[str(atomtype)].get_drude_charge()
			#add atomtype to list of all atomtypes if not already there
			if str(atomtype) not in allAtomtypes:
				allAtomtypes.append(str(atomtype))
			#if atomtype is anisotropic, add to list of anisotropic atomtypes
			if self.anisoFlag and atomtype not in anisoAtomtypes:
				anisoAtomtypes.append(str(atomtype))

	class Atomtype:
		#creates object to store parameters for an atomtype within a dimer
		def __init__(self, inputDict, anisoFlag = False):
			self.parameters = {}
			self.parameters["A"] = self.Parameters(inputDict['A'], 'A')
			self.parameters["B"] = inputDict['B'] #since only one B value for each atomtype
			self.parameters["C"] = self.Parameters(inputDict['C'], 'C')
			self.drude_charge = inputDict['drude_charge']
			self.anisoFlag = anisoFlag
			self.sphericalHarmonics = []
			if anisoFlag:
				self.parameters["aniso"] = self.Parameters(inputDict["aniso"], 'aniso', inputDict["sph_harm"])
				self.sphericalHarmonics = inputDict["sph_harm"]
			else:
				self.parameters["aniso"] = {}

		def __repr__(self):
			return str(self.parameters)

		def __getitem__(self, paramType):
			return self.parameters[paramType]

		def isAniso(self):
			"""return true if atomtype is anisotropic"""
			return self.anisoFlag

		def getSphericalHarmonics(self):
			return self.sphericalHarmonics

		def get_drude_charge(self):
			return self.drude_charge

		class Parameters:
			#creates dictionary of parameter values for individual energy terms
			def __init__(self, inputList, paramType, sphericalHarmonics = []):
				self.type = paramType
				self.values = {}

				if paramType == 'A':
					self.values = self.populateAvalues(inputList)

				elif paramType == 'aniso':
					self.values = self.populateAnisoValues(inputList, sphericalHarmonics)

				elif paramType == 'C':
					self.values = self.populateCvalues(inputList)

			def __repr__(self):
				return str(self.values)

			def __getitem__(self, energyTerm):
				if self.type == "aniso":
					return self.values[energyTerm]
				else:
					return self.values[energyTerm]


			def populateAvalues(self, inputList):
				"""takes in list of A parameter values and returns as dictionary

				Parameters
				----------
				inputList: list of A parameter values, acquired from JSON input dictionary of outer atomtype object

				Returns 
				-------
				Dictionary of A parameter values with energy component as key

				"""
				output = {}
				output["Exchange"] = inputList[0]
				output["Electrostatics"] = inputList[1]
				output["Induction"] = inputList[2]
				output["Dhf"] = inputList[3]
				output["Dispersion"] = inputList[4]
				return output

			def populateAnisoValues(self, inputList, sphericalHarmonics):
				"""takes in list of aniso parameter values and returns as dictionary

				Parameters
				----------
				inputList: list of aniso parameter values, acquired from JSON input dictionary of outer atomtype object

				Returns 
				-------
				Dictionary of aniso parameter values with spherical harmonic as key

				"""
				output = {}
				for i in range(len(energyComponents)):
					innerOutput = {}
					for j in range(len(sphericalHarmonics)):
						innerOutput[str(sphericalHarmonics[j])] = inputList[i][j]
					output[energyComponents[i]] = innerOutput
				return output

			def populateCvalues(self, inputList):
				"""takes in list of C parameter values and returns as dictionary

				Parameters
				----------
				inputList: list of C parameter values, acquired from JSON input dictionary of outer atomtype object

				Returns 
				-------
				Dictionary of C parameter values with order as key

				"""
				output = {}
				output["C6"] = inputList[0]
				output["C8"] = inputList[1]
				output["C10"] = inputList[2]
				output["C12"] = inputList[3]
				return output

##########################
# JSON interface functions
##########################

def readJson(filePath):
	"""read in JSON file and return dictionary from file

	Parameters
	----------
	filePath: path to JSON file to read

	Returns 
	-------
	Dictionary format of JSON file

	"""
	with open (filePath, 'r') as json_file:
		json_data =json.load(json_file)
	return json_data

def importParameters():
	"""import all parameters from JSON files. Finds all .constraints files in directory of initial constrait files, converts JSON to dictionary and stores in global dictionary of initial parameters with dimer name as key

	Parameters
	----------
	None

	Returns 
	-------
	None

	"""
	inputDir = generalizer_settings.constraints_input_directory

	try:
		os.chdir(scriptDir + inputDir)
		#get list of all constraint (ab initio) files
		file_list = glob.glob("*.constraints")
		file_list.sort()
		
		global inputJsons
		inputJsons = {}
		prefix = generalizer_settings.constraints_prefix
		suffix = generalizer_settings.constraints_suffix
		for dimer in file_list:
			dimerJson = readJson(dimer)
			#get name of dimer
			name = dimer.replace(prefix, '')
			name = name.replace(suffix, '')
			dimerName = name 
			inputJsons[dimerName] = dimerJson

		os.chdir(scriptDir)

	except OSError:
		print("need to create abInitioConstraints directory!")

##################################################################################################################
# functions for use with imported inital parameters
# to be used when we already have parameters we want to start minimize() with, do not need to average from library
##################################################################################################################

def useInitialParamsFile(dimerName):
	"""check to see if .params file already exists, import initial values if so

	Parameters
	----------
	dimerName: name of dimer to find .params file for

	Returns 
	-------
	None, but stores given parameters in initialMinimizeValues

	"""
	if os.path.isfile(dimerName + ".params"):
		global initialMinimizeValues
		initialMinimizeValues = dimecularParameters(readJson(dimerName + ".params"), dimerName) #need to make this a list to interface with minimize

def minimizeDictionaryToList():
	"""convert dictionary of imported initial values to list

	Parameters
	----------
	None

	Returns 
	-------
	None, but stores given parameters in initialValuesList to be given to minimize. List will be of the format [Atomtype1 A,B,C,Aniso; Atomtype2 A,B,C,Aniso; ...]

	"""
	output = []

	#loop over all atomtypes
	for atomtype in allAtomtypes:
		currentAtomtype = initialMinimizeValues[atomtype]
		#loop over A energy components
		aParams = currentAtomtype["A"]
		for component in energyComponents:
			output.append(aParams[component])
		#B
		output.append(currentAtomtype["B"])
		#loop over C energy components
		cParams = currentAtomtype["C"]
		for component in dispersionTerms:
			output.append(cParams[component])
		#loop over aniso energy components
		anisoParams = currentAtomtype["aniso"]
		sphericalHarmonics = currentAtomtype.getSphericalHarmonics()
		for component in energyComponents:
			#loop over aniso spherical harmonics
			for sh in sphericalHarmonics:
				output.append(anisoParams[component][sh])

	global initialValuesList
	initialValuesList = initialValuesList + output

######################################################
# fuctions to use when given library of ab initio FF
######################################################

def importdimecularParameters():
	"""import all ab initio parameters

	Parameters
	----------
	None

	Returns 
	-------
	None, but stores atomtypes of imported dimers in global list of atomtypes

	"""
	importParameters()
	for dim in list(inputJsons.keys()):
		print(("Importing dimecular parameters: " + dim))
		dimers[dim] = dimecularParameters(inputJsons[dim], dim)
		dimAtomtypes = dimers[dim].atomtypesList
		print((dim + " parameters successfully imported!"))

def averageAvalues():
	"""return averages of all A values in dictionary

	Parameters
	----------
	None

	Returns 
	-------
	Dictionary of average A values for all atomtypes, averaged from all previously imported dimers. Nested dictionary of {atomtype {energy component : value}} 

	"""
	allAvalsRaw = {}
	allAvalsOut = {}
	for dim in list(dimers.keys()):
		currentdimer = dimers[dim]
		for atomtype in allAtomtypes:
			if atomtype in allAvalsRaw: #if we already have this atomtype, add parameters to existing list
				for component in energyComponents:
					try: #try/except allows for multiple atomtypes in different input files
						allAvalsRaw[atomtype][component].append(currentdimer[atomtype]["A"][component])
					except KeyError:
						pass
			else:
				componentsDict = {}
				for component in energyComponents:
					try:
						componentsDict[component] = [currentdimer[atomtype]["A"][component]]
						allAvalsRaw[atomtype] = componentsDict
					except KeyError:
						pass
	for atomtype in allAvalsRaw: #calculate averages from list for each atomtyp for each energy type
		componentAvgs = {}
		for component in energyComponents:
			componentAvgs[component] = sum(allAvalsRaw[atomtype][component])/len(allAvalsRaw[atomtype][component])
		allAvalsOut[atomtype] = componentAvgs
	return allAvalsOut

def averageBvalues():
	"""return average of all B values

	Parameters
	----------
	None

	Returns 
	-------
	Dictionary of average B values for all atomtypes, averaged from all previously imported dimers.

	"""
	allBvalsRaw = {}
	allBvalsOut = {}
	for dim in list(dimers.keys()):
		currentdimer = dimers[dim]
		for atomtype in allAtomtypes:
			if atomtype in allBvalsRaw: #if we already have this atomtype, add parameters to existing list
				try:
					allBvalsRaw[atomtype].append(currentdimer[atomtype]["B"])
				except KeyError:
					pass
			else:
				try:
					allBvalsRaw[atomtype] = [currentdimer[atomtype]["B"]]
				except KeyError:
					pass
	for atomtype in allBvalsRaw:
		allBvalsOut[atomtype] = sum(allBvalsRaw[atomtype])/len(allBvalsRaw[atomtype])
	return allBvalsOut

def averageCvalues():
	"""return averages of all C values in dictionary

	Parameters
	----------
	None

	Returns 
	-------
	Dictionary of average C values for all atomtypes, averaged from all previously imported dimers. Nested dictionary of {atomtype {order : value}} 

	"""
	allCvalsRaw = {}
	allCvalsOut = {}
	for dim in list(dimers.keys()):
		currentdimer = dimers[dim]
		for atomtype in allAtomtypes:
			if atomtype in allCvalsRaw: #if we already have this atomtype, add parameters to existing list
				for term in dispersionTerms:
					try:
						allCvalsRaw[atomtype][term].append(currentdimer[atomtype]["C"][term])
					except KeyError:
						pass
			else: #if we have not seen this atomtype yet
				termsDict = {}
				for term in dispersionTerms:
					try:
						termsDict[term] = [currentdimer[atomtype]["C"][term]]
						allCvalsRaw[atomtype] = termsDict
					except KeyError:
						pass
	for atomtype in allCvalsRaw:
		termAvgs = {}
		for term in dispersionTerms:
			termAvgs[term] = sum(allCvalsRaw[atomtype][term])/len(allCvalsRaw[atomtype][term])
		allCvalsOut[atomtype] = termAvgs
	return allCvalsOut

def averageAnisoValues():
	"""return averages of all aniso values in dictionary

	Parameters
	----------
	None

	Returns 
	-------
	Dictionary of average Aniso values for all atomtypes, averaged from all previously imported dimers. Nested dictionary of {atomtype {energy component {spherical harmonic : value}}} 

	"""
	#make dictionaries to store parameters
	allAnisoValsRaw = {}
	allAnisoValsOut = {}

	#get parameters for each ab initio dimer
	for dim in list(dimers.keys()):
		currentdimer = dimers[dim]
		#for each atomtype, check if already in dictionary
		for atomtype in allAtomtypes:
			try:
				currentAtomtype = currentdimer[atomtype]
				if currentAtomtype.isAniso(): #check if atomtype is anisotropic
					currentAnisoParams = currentAtomtype["aniso"]
					currentSphericalHarmonics = currentAtomtype.getSphericalHarmonics()
					if atomtype in allAnisoValsRaw: #check if we already have a list for this atomtype
						for energycomponent in energyComponents:
							currentValues = currentAnisoParams[energycomponent]
							for sh in currentSphericalHarmonics:
								allAnisoValsRaw[atomtype][energycomponent][sh].append(currentValues[sh])
					else:
						energycomponentValues = {}
						for energycomponent in energyComponents:
							currentValues = currentAnisoParams[energycomponent]
							shValues = {} #store values for spherical harmonics for a given energy component
							for sh in currentSphericalHarmonics:
								shValues[sh] = [currentValues[sh]]
							energycomponentValues[energycomponent] = shValues
						allAnisoValsRaw[atomtype] = energycomponentValues
			except KeyError:
				pass
	#perform averaging
	for atomtype in allAnisoValsRaw:
		componentsAvg = {}
		for component in energyComponents:
			energyComponentsAvgs = {}
			currentValues = allAnisoValsRaw[atomtype][component]
			for sh in currentValues:
				energyComponentsAvgs[str(sh)] = sum(currentValues[sh])/len(currentValues[sh])
			componentsAvg[str(component)] = energyComponentsAvgs
		allAnisoValsOut[atomtype] = componentsAvg
	return allAnisoValsOut

def getInitalAverageParameters():
	"""calculate average parameters and store in initialAverageValues

	Parameters
	----------

	None

	Returns
	-------

	None, but stores average parameters in initialAverageValues dictionary with parameter type (i.e. A, B, C, aniso) as key
	"""
	#calculate average values
	aParams = averageAvalues()
	bParams = averageBvalues()
	cParams = averageCvalues()
	anisoParams = averageAnisoValues()

	#store average values in initialAverageValues dict
	global initialAverageValues
	initialAverageValues["A"] = aParams
	initialAverageValues["B"] = bParams
	initialAverageValues["C"] = cParams
	initialAverageValues["aniso"] = anisoParams

def averageParamsDictToList():
	"""takes initialAverageValues dictionary of averaged inital guesses and makes into list of format [Atomtype1 A,B,C,Aniso; Atomtype2 A,B,C,Aniso; ...]

	Parameters
	----------

	None

	Returns
	-------

	None, but stores list in initialValuesList
	"""

	output = []
	#loop over atomtypes
	allAtomtypes.sort()
	for atomtype in allAtomtypes:
		nparams = 0
		#loop over A components
		aParams = initialAverageValues["A"][atomtype]
		for component in energyComponents:
			output.append(aParams[component])
			nparams += 1
		#B
		output.append(initialAverageValues["B"][atomtype])
		atomtype_B_params[atomtype] = initialAverageValues["B"][atomtype]
		nparams += 1
		#loop over C components
		cParams = initialAverageValues["C"][atomtype]
		for component in dispersionTerms:
			output.append(cParams[component])
			nparams += 1
		#loop over aniso components
		try:
			anisoParams = initialAverageValues["aniso"][atomtype]
			sphericalHarmonics = atomtypeSphericalHarmonics[atomtype]
			for component in energyComponents:
				#loop over spherical harmonics
				for sh in sphericalHarmonics:
					output.append(anisoParams[component][sh])
					nparams += 1
		except KeyError: #catch if atomtype is not anisotropic
			pass
		atomtype_nparams[atomtype] = nparams

	global initialValuesList
	initialValuesList = initialValuesList + output

def write_global_exponents_file():
	""" writes initial B params to .exp file so that all dimers get same initial B param values

	Parameters
	----------
	none

	Returns
	-------
	None
	"""
	generalizer_input_directory = generalizer_settings.generalizer_input_directory

	os.chdir(scriptDir + generalizer_input_directory)
	with open("global_exponents.exp", 'w') as f:
		for atomtype in allAtomtypes:
			current_atomtype_B_param = atomtype_B_params[atomtype]
			f.write(atomtype + '\t' + str(current_atomtype_B_param) + '\n')
	os.chdir(scriptDir)

##########################################################################################################################
# metaPOInter and helper functions
# metaPOInter interfaces with POInter to get lsq_error, dlsq_error, and FF energy for each parameter set for each dimer
##########################################################################################################################

def map_params(inputList, fileName):
	"""writes .constraints file in JSON format for metaPOInter from a given input list of parameters

	Parameters
	----------

	inputList: list of parameters in format [Atomtype1 A,B,C,Aniso; Atomtype2 A,B,C,Aniso; ...]

	Returns:

	None, but writes .constraints file for metaPOInter
	"""

	#make sure we have a python list and not a numpy array
	if isinstance(inputList, np.ndarray):
		inputList = inputList.tolist()

	output = {}
	i = 0 #use to keep track of where we are in the parameter list, updated after adding each parameter or set of parameters

	for atomtype in allAtomtypes:
		atomtypeParams = {}
		atomtypeParams["params_order"] = ["Exchange", "Electrostatics", "Induction", "Dhf", "Dispersion"]
		#get A params
		atomtypeParams["A"] = inputList[i:i+5]
		i = i + 5
		#get B param
		atomtypeParams["B"] = inputList[i]
		i += 1
		#get C params
		atomtypeParams["C"] = inputList[i: i + len(dispersionTerms)]
		i = i + len(dispersionTerms)
		#check if atomtype is anisotropic. If so, fetch spherical harmonics and save in atomtypeParams. Also save aniso A values, with empty lists if not anisotropic
		if atomtype in anisoAtomtypes:
			atomtypeParams["sph_harm"] = atomtypeSphericalHarmonics[atomtype]
			anisoList = []
			for j in range(5):
				anisoList.append(inputList[i : i + len(atomtypeSphericalHarmonics[atomtype])])
				i = i + len(atomtypeSphericalHarmonics[atomtype])
			atomtypeParams["aniso"] = anisoList
		else:
			atomtypeParams["sph_harm"] = []
			atomtypeParams["aniso"] = [[],[],[],[],[]]
		#set drude charge
		global atomtype_drude_charges
		atomtypeParams["drude_charge"] = atomtype_drude_charges[atomtype]
		output[atomtype] = atomtypeParams

	#write .constraints file with parameters from list
	with open(fileName + ".constraints", 'w') as f:
		json.dump(output, f, indent = 4)

def calc_harmonic_constraint_penalty(params, k = 1e-5):
	"""
	calculates harmonic penalty for parameters differing from initial constraints

	Parameters:
	-----------

	params: tuple or list, current parameters

	k: float, Weighitng parameter; higher k values employ more stringent
            constraints on the values the parameters can take. k=1e-5

    also used initialValuesList global variable

	returns:
	--------

	harmonic_error : float, Energy penalty associated with employed harmonic constraints.

	parameter_errors: list, contains floats of errors for each individual parameter
        
    dharmonic_error : list
            Derivatives of the harmonic error with respect to each
            parameter.
    """
	harmonic_error = 0
	parameter_errors = []
	dharmonic_error = []
	dharmonic_error = [ 0 for _ in params]
	for i in range(len(params)):
		b = params[i]
		b0 = initialValuesList[i]
		param_error = k*(b - b0)**2
		harmonic_error += param_error
		parameter_errors.append(param_error)
		dharmonic_error[i] = 2*k*(b - b0)*b0

	return harmonic_error, dharmonic_error, parameter_errors

def make_bounds_list(parameterList, include_B = True):
	"""
	makes list of tuples containing bounds for parameters. Order corresponds with parameter order [Atomtype1 A,B,C,Aniso; Atomtype2 A,B,C,Aniso; ...]

	Parameters: 
	----------

	parameterList: list, input parameter values
	include_B: boolean, dictates whether to include a bound for the B param

	Returns:
	--------

	bounds_list: list of bound tuples

	"""
	i = 0 #keep track of position in bounds_list
	#TODO: remove i

	bounds_list = []

	global fitting_component

	for atomtype in allAtomtypes:
		#A params
		if fitting_component != 4:
			bounds_list.append((0.0,1e3))
			i += 1
		else:
			bounds_list.append((0.7,1.3))
			i += 1
		#B param
		if include_B:
			bounds_list.append((1e-2,1e2))
			i += 1
		#aniso
		if atomtype in anisoAtomtypes:
			sphericalHarmonics = atomtypeSphericalHarmonics[atomtype]
			for sh in sphericalHarmonics:
				bounds_list.append((-1.0,1.0))
				i += 1

	return bounds_list

def get_fit_parameters(all_params, fit_dimer):
	"""
	Fetches the parameters to be fit for a given component.

	Parameters
	----------

	all_params: list, parameters for all atomtypes for the given component, in form [atomtype1 A, B, aniso; atomtype2 A, B, aniso...]
	fit_dimer: string, name of dimer being worked with

	Returns
	-------

	fit_params: list, params (for the current component) to be input to POInter in format [atomtype1 A, aniso, B; atomtype2 A, aniso, B], only including atomtypes for that dimer

	"""

	current_dimer_atomtypes = dimer_atomtypes[fit_dimer]

	fit_params = []

	#keep track of where we are in the input list. Will be positioned at the front of whatever set of parameters we are trying to get
	input_placeholder = 0

	global fitting_component

	#get parameters for each atomtype
	for pot_atomtype in allAtomtypes:
		# determine number of parameters
		n_params = 1 # for A param
		if pot_atomtype in anisoAtomtypes:
			n_spherical_harmonics = len(atomtypeSphericalHarmonics[pot_atomtype])
			n_params += n_spherical_harmonics
		if fitting_component == 0 and generalizer_settings.fit_b:
			n_params += 1 #if we are including the B param (exchange)

		#if our dimer has this atomtype, add these parameters to the output list	
		if pot_atomtype in current_dimer_atomtypes:
			atomtype_params = all_params[input_placeholder : input_placeholder + n_params]
			for param in atomtype_params:
				fit_params.append(param)

		input_placeholder += n_params

	return fit_params

def remove_C_params(parameterList):
	"""
	removes C parameters from parameterList of all parameters (for all atomtypes) in format [atomtype1 A, B, C, aniso; atomtype2 A, B, C, aniso...]

	parameters
	----------
	parameterList: list, all parameters in format [atomtype1 A, B, C, aniso; atomtype2 A, B, C, aniso...]

	returns
	-------
	outputList: list, all A, B, and aniso parameters in format [atomtype1 A, B, aniso; atomtype2 A, B, aniso...]
	"""

	input_index = 0
	outputList = parameterList
	for atomtype in allAtomtypes:
		#move index counter past A and B params
		input_index += 6
		#add C params to global dict and delete from list
		global atomtype_C_params
		atomtype_C_params[atomtype] = outputList[input_index : input_index + 4]
		del outputList[input_index : input_index + 4]
		#move index counter past aniso params
		if atomtype in anisoAtomtypes:
			n_aniso_params = len(atomtypeSphericalHarmonics[atomtype])*5
			input_index += n_aniso_params
	return outputList

def get_component_parameters(all_params):
	""" will take all parameters and return only those for a given component

	parameters
	----------
	all_params: list, all parameters (except C params) for all atomtypes and energy components in format [atomtype1 A, B, aniso; atomtype2 A, B, aniso...]

	returns
	-------
	output_params: list, has parameters for a given energy component for all atomtypes in format [atomtype1 A, aniso, B; atomtype2 A, aniso, B]
	"""
	input_index = 0
	output_params = []
	global fitting_component
	print(all_params)

	for atomtype in allAtomtypes:
		n_params = atomtype_nparams[atomtype] - 4 #subtract out C params
		#get A param
		output_params.append(all_params[input_index + fitting_component])
		input_index += 5
		#get aniso params
		if atomtype in anisoAtomtypes:
			init_input_index = input_index
			n_spherical_harmonics = len(atomtypeSphericalHarmonics[atomtype])
			for i in range(fitting_component):
				input_index += n_spherical_harmonics
			aniso_params = all_params[input_index + 1 : input_index + 1 + n_spherical_harmonics]
			for aaniso in aniso_params:
				output_params.append(aaniso)
			input_index = init_input_index
		#get B if required (exchange)
		if fitting_component == 0 and generalizer_settings.fit_b:
			output_params.append(1.0)

		#move index counter to start of next atomtype
		input_index -= 5
		input_index += n_params

	return output_params

def add_dlsq_error(all_dlsq_error, new_dlsq_error, dimerName):
	""" add dlsq_error returned by POInter to list of all dlsq_error

	parameters
	----------
	all_dlsq_error: list, all dlsq_error for all atomtypes and parameters in format corresponding to [atomtype1 A, aniso, B; atomtype2 A, aniso, B]
	new_dlsq_error: list, dlsq_error returned from POInter to be added to all_dlsq_error (in same format at all_dlsq_error, just without some atomtypes not found in that dimer)
	dimerName: name of dimer

	returns
	-------
	all_dlsq_error updated with new_dlsq_error added
	"""

	global fitting_component
	dim_atomtypes = dimer_atomtypes[dimerName]
	all_dlsq_error_index = 0

	#add new dlsq_error to all_dlsq_error
	new_dlsq_error_index = 0

	for pot_atomtype in allAtomtypes:
		n_params = atomtype_nparams[pot_atomtype] - 4 # excluding C params which are not fit
		if fitting_component != 0:
			n_params -= 1 # subtract B param if not fitting exchange
		n_params = int(n_params / 5) # to find numer of parameters for a single component energy
		if pot_atomtype in dim_atomtypes:
			for i in range(n_params):
				all_dlsq_error[all_dlsq_error_index + i] += new_dlsq_error[new_dlsq_error_index + i]
			new_dlsq_error_index += n_params
		all_dlsq_error_index += n_params

	return all_dlsq_error

def write_output_params_to_list(parameterDict):
	"""
	write parameters from optimization function to a list in format [atomtype1 A, B, C, aniso; atomtype2 A, B, C, aniso...]

	parameters
	----------
	parameterDict: dictionary, parameters to be placed in list

	returns
	-------
	outputList: list, in format [atomtype1 A, B, C, aniso; atomtype2 A, B, C, aniso...]
	"""
	outputList = []
	input_atomtype_i = 0
	shift = 1
	exchange_shift = 1
	for atomtype in allAtomtypes:
		
		if atomtype in anisoAtomtypes:
			n_spherical_harmonics = len(atomtypeSphericalHarmonics[atomtype])
		else:
			n_spherical_harmonics = 0

		shift = 1 + n_spherical_harmonics
		#get A params	
		for component in energyComponents:
			if component == "Exchange" and generalizer_settings.fit_b:
				b_scale = parameterDict[component][input_atomtype_i + n_spherical_harmonics + exchange_shift]
				exchange_shift += 1
			else:
				b_scale = 1.0
			try:
				a_param = parameterDict[component][input_atomtype_i]
				outputList.append(a_param)
			except IndexError:
				if component == "Dispersion": #add dispersion params for anisotropic atomtypes
					a_param = [1.0 for sh in range(n_spherical_harmonics)]
					outputList += a_param
				else:
					raise IndexError
		#get B param
		b_init_param = atomtype_B_params[atomtype]
		new_b_param = b_init_param * b_scale
		outputList.append(new_b_param)
		#get C params
		c_params = atomtype_C_params[atomtype]
		for param in c_params:
			outputList.append(param)
		#get aniso
		if atomtype in anisoAtomtypes:
			for component in energyComponents:
				aniso_params = parameterDict[component][input_atomtype_i + 1 : input_atomtype_i + 1 + n_spherical_harmonics]
				for param in aniso_params:
					outputList.append(param)

		input_atomtype_i += shift

	return outputList

def get_updated_B_from_exchange(exchange_params):
	""" get B parameters from list output from exchange fitting


	"""
	exchange_shift = 1
	input_atomtype_i = 0
	shift = 1
	output_B_params = []

	for atomtype in allAtomtypes:
		if atomtype in anisoAtomtypes:
			n_spherical_harmonics = len(atomtypeSphericalHarmonics[atomtype])
		else:
			n_spherical_harmonics = 0

		b_scale = exchange_params[input_atomtype_i + n_spherical_harmonics + exchange_shift]
		exchange_shift += 1

		b_init_param = atomtype_B_params[atomtype]
		new_b_param = b_init_param * b_scale
		output_B_params.append(new_b_param)

		input_atomtype_i += shift

	return output_B_params

def calc_disp_lsq_error():
	"""
	calculate the lsq_error for the dispersion component

	parameters
	----------
	none

	returns
	-------
	total_lsq_error: float, sum of all lsq_error for each individual dimecular dispersion component
	"""
	#set fitting component for dispersion
	global fitting_component
	fitting_component = 4

	total_lsq_error = 0.0

	for dim in list(dimers.keys()):
		#get POInter object and setting fitting_component
		pointerModel = dimer_POInter_objects[dim]
		pointerModel.component = fitting_component
		#get qm energies
		pointerModel.qm_fit_energy = np.array(pointerModel.subtract_hard_constraint_energy())
		#set POInter instance variables
		pointerModel.fit_bii = False
		pointerModel.n_isotropic_params = 1
		n_aiso = pointerModel.n_isotropic_params if not pointerModel.fit_bii else pointerModel.n_isotropic_params - 1
		n_aaniso = n_aiso
		#calc lsq_error, dlsq_error
		current_dimer_atomtypes = dimer_atomtypes[dim]
		currentParams = tuple(1.0 for atomtype in current_dimer_atomtypes)
		print(("Getting lsq_error for component", fitting_component, "with current parameters", currentParams))
		pointerModel.get_num_eij = pointerModel.generate_num_eij(currentParams)
		pointerOutput = pointerModel.calc_leastsq_ff_fit(currentParams)
		total_lsq_error += pointerOutput[0]

	return total_lsq_error

def make_output_directories():
	"""make directories to store output files from parameter generalizer_settings

	Parameters
	----------
	None

	Returns
	-------
	None
	"""
	global scriptDir
	output_directory = generalizer_settings.output_directory
	if not os.path.isdir(scriptDir + output_directory):
		os.mkdir(scriptDir + output_directory)
	os.chdir(scriptDir + output_directory)
	for dim in list(dimers.keys()):
		if not os.path.isdir(scriptDir + output_directory + dim):
			os.mkdir(scriptDir + output_directory + dim)
	os.chdir(scriptDir)

def metaPOInter(parameterList):
	"""
	Interface function with POInter code. Called by minimize(), will get lsq_error and dlsq_error for each dimer in library, returns total lsq_error and total dlsq_error

	Parameters
	----------

	None

	Returns
	-------

	total_lsq_error: float, sum lsq_error for entire library of dimers

	dlsq_error: list, list of total dlsq_error's for each parameter being optimized

	"""

	total_lsq_error = 0
	dlsq_error = [0.0 for i in range(len(parameterList))]
	global fitting_component

	#get lsq_error and dlsq_error for each dimer from POInter, add to running totals
	generalizer_input_directory = generalizer_settings.generalizer_input_directory
	os.chdir(scriptDir + generalizer_input_directory)
	for dim in list(dimers.keys()):
		print(("Calculating for " + str(dim)))
		if dim not in dimer_POInter_objects:
			saptFile = dim + ".sapt"
			#make POInter object
			pointerModel = Pointer.FitFFParameters(fit = False, energy_file = saptFile)
			pointerModel.inputdir = os.getcwd()

			#call required POInter start-up functions
			exp_file = generalizer_settings.exponent_file
			monomers = dim.split('_')
			kwargs = {"mon1": monomers[0], "mon2": monomers[1], "exponent_files" : exp_file}
			pointerModel.read_settings(['default'], kwargs)
			pointerModel.read_energies()
			pointerModel.read_params()
			pointerModel.initialize_parameters()
			#set required POInter instance variables
			pointerModel.n_isotropic_params = pointerModel.default_n_isotropic_params
			pointerModel.final_energy_call = True
			pointerModel.component = fitting_component

			#make sure the order of the parameters in POInter and metaPOInter match
			pointerModel.fit_isotropic_atomtypes.sort()
			pointerModel.fit_anisotropic_atomtypes.sort()
			#store POInter object for later use by subsequent calls to metaPOInter
			dimer_POInter_objects[dim] = pointerModel
			#set dimer drude flag to false
			drude_flags[dim] = False
		else:
			pointerModel = dimer_POInter_objects[dim]
			pointerModel.component = fitting_component

		#get drude energy if component 1
		if fitting_component == 1 and not drude_flags[dim]:
			pointerModel.get_drude_oscillator_energy()
			drude_flags[dim] = True
		#subtract hard constraints energy from energy to be fit
		dim_qm_dict = qm_energies[dim]
		if fitting_component not in dim_qm_dict:
			pointerModel.qm_fit_energy = np.array(pointerModel.subtract_hard_constraint_energy())
			qm_energies[dim][fitting_component] = pointerModel.qm_fit_energy
		else:
			pointerModel.qm_fit_energy = qm_energies[dim][fitting_component]

		#set fitting component and B coeff fitting (if exchange component)
		if fitting_component == 0 and generalizer_settings.fit_b:
			pointerModel.fit_bii = True
			pointerModel.n_isotropic_params = 2

			#set POInter instance variables required for fitting B param

			i_bparams = []

			param_i = 0
			for atomtype in dimer_atomtypes[dim]:
				if atomtype in anisoAtomtypes:
					num_aniso = len(atomtypeSphericalHarmonics[atomtype])
				else:
					num_aniso = 0
				num_A_aniso = 1 + num_aniso
				b_pos = param_i + num_A_aniso
				i_bparams.append(b_pos)
				param_i += num_A_aniso + 1

			pointerModel.i_bparams = i_bparams

		else:
			pointerModel.fit_bii = False
			pointerModel.n_isotropic_params = 1
			n_aiso = pointerModel.n_isotropic_params if not pointerModel.fit_bii else pointerModel.n_isotropic_params - 1
			n_aaniso = n_aiso

		#calc lsq_error, dlsq_error
		currentParams = tuple(get_fit_parameters(parameterList, dim))
		print(("Fitting component", fitting_component, "with current parameters", currentParams))
		pointerModel.get_num_eij = pointerModel.generate_num_eij(currentParams)
		pointerOutput = pointerModel.calc_leastsq_ff_fit(currentParams)


		print()
		print(("leastsq error: ", pointerOutput[0]))
		print()

		#add lsq_error, dlsq_error to totals

		total_lsq_error += pointerOutput[0]
		dlsq_error_components = pointerOutput[1]

		dlsq_error = add_dlsq_error(dlsq_error, dlsq_error_components, dim)	

	os.chdir(scriptDir)

	#store total_lsq_error
	component_lsq_error[energyComponents[fitting_component]] = total_lsq_error

	return total_lsq_error, np.asarray(dlsq_error)

def optimzeGeneralParameters():
	"""
	Primary function of parameterGeneralizer module. Will call scipy minimize to find optimal parameters for library.

	Parameters
	----------

	None

	Returns
	-------
	Global lsq_error for library
	"""
	#initialize qm energies dictionaries for each dimer
	global qm_energies
	qm_energies = {dim : {} for dim in list(dimers.keys())}

	#set parameters for minimize function
	maxiter = 5000
	pgtol = -1e-17
	ftol = 1e-17

	#get parameter bounds for minimization
	all_bnds = make_bounds_list(initialValuesList, True)
	a_bnds = make_bounds_list(initialValuesList, False)

	#optimize A parameters for each energy component, store in final_params list
	parameterList = initialValuesList
	metaPOInter_initial_list = remove_C_params(parameterList)
	global fitting_component
	while fitting_component < 4:

		#make list of bounds for component to pass down to POInter. Only include B bound if fitting exchange
		if fitting_component == 0 and generalizer_settings.fit_b:
			bnds = all_bnds
		else:
			bnds = a_bnds

		#make list of parameters for component to pass down to POInter
		#get only parameters for the component we are fitting

		#make sure all dimers get right exponent file
		#TODO: make sure this line is in the right place
		write_global_exponents_file()

		metaPOInter_input = get_component_parameters(metaPOInter_initial_list)

		res = minimize(metaPOInter,metaPOInter_input,method='L-BFGS-B',\
							jac=True,\
							options={'disp':True,'gtol':pgtol,'ftol':ftol,'maxiter':maxiter},\
							bounds = bnds)

		popt = res.x
		success = res.success
		message = res.message
		print(("Optimized parameters for component " + str(fitting_component) + " : " + str(popt)))

		if not res.success:
			print('Warning! Optimizer did not terminate successfully, and quit with the following error message:')
			print()
			print((res.message))
			print()
			return

		global optimized_params
		optimized_params[energyComponents[fitting_component]] = popt

		#initialize ff energy dictionary
		if fitting_component == 0:
			global ff_energies
			ff_energies = {dim : np.zeros_like(dimer_POInter_objects[dim].qm_energy[6]) for dim in list(dimers.keys())}

		#calculate ff energy for each dimer and write output files
		for dim in list(dimers.keys()):
			pointerModel = dimer_POInter_objects[dim]
			pointerModel.component = fitting_component
			pointerModel.get_num_eij = pointerModel.generate_num_eij(popt)
			pointerModel.output_params(popt)
			ff_fit_energy = pointerModel.calc_ff_energy(popt)[0]
			ff_energy = np.array(pointerModel.qm_energy[fitting_component]) - np.array(qm_energies[dim][fitting_component]) + ff_fit_energy #adding back in hard constraint energy
			pointerModel.rms_error = pointerModel.calc_rmse(ff_energy)
			pointerModel.weighted_absolute_error = pointerModel.calc_mse(ff_energy, cutoff=pointerModel.weighted_rmse_cutoff)
			pointerModel.weighted_rms_error = pointerModel.calc_rmse(ff_energy, cutoff=pointerModel.weighted_rmse_cutoff)
			pointerModel.lsq_error = pointerModel.calc_leastsq_ff_fit(popt)[0]
			#write output files
			print(("writing output for component " + str(fitting_component) + " for " + dim))
			global scriptDir
			os.chdir(scriptDir + "/dim_fits/" + dim)
			pointerModel.write_output_file(success,message)
			pointerModel.write_energy_file(ff_energy)

			total_ff_energy = ff_energies[dim]
			total_ff_energy += ff_energy

			os.chdir(scriptDir)

		#reset new B params in POInter objects after fitting exchange component
		if fitting_component == 0 and generalizer_settings.fit_b:
			updated_B_params = get_updated_B_from_exchange(popt)
			for dim in list(dimers.keys()):
				i = 0
				for atomtype in allAtomtypes:
					try:
						pointerModel.params[atomtype][0]['B'] = updated_B_params[i]
					except KeyError: #catch if atomtype not in current dimer
						pass
					i += 1
		
		fitting_component += 1

	#write output JSON file with optimized parameters
	optimized_params["Dispersion"] = [1.0 for atomtype in allAtomtypes]
	final_params_list = write_output_params_to_list(optimized_params)
	map_params(final_params_list, "optimized_params")

	#calc dispersion energy and total energy
	for dim in list(dimers.keys()):
		os.chdir(scriptDir + "/dim_fits/" + dim)
		#dispersion
		fitting_component = 4
		pointerModel = dimer_POInter_objects[dim]
		pointerModel.component = fitting_component
		if not scale_iso_disp:
			pointerModel.default_n_isotropic_params -= 1
			ff_energy = pointerModel.fit_component_parameters()
			pointerModel.write_energy_file(ff_energy)
			print((dim + " energy file for component " + str(fitting_component) + " successfully written"))
			total_ff_energy = ff_energies[dim]
			total_ff_energy += ff_energy
		else:
			raise NotImplementedError

		#total energy
		pointerModel.component = 6
		ff_energy = ff_energies[dim]
		qm_energy = pointerModel.qm_energy[pointerModel.component]
		weight = Pointer.functional_forms.weight(qm_energy, pointerModel.eff_mu, pointerModel.eff_kt)
		pointerModel.lsq_error =  np.sum(weight*(ff_energy - qm_energy)**2)
		pointerModel.write_energy_file(ff_energy)
		pointerModel.rms_error = pointerModel.calc_rmse(ff_energy)
		pointerModel.weighted_rms_error = pointerModel.calc_rmse(ff_energy, cutoff=pointerModel.weighted_rmse_cutoff)
		pointerModel.weighted_absolute_error = pointerModel.calc_mse(ff_energy, cutoff=pointerModel.weighted_rmse_cutoff)
		pointerModel.write_output_file()
		print(("Total energy file for " + dim + " successfully written"))

		os.chdir(scriptDir)

	#calc dispersion lsq_error and store
	disp_lsq_error = calc_disp_lsq_error()
	component_lsq_error[energyComponents[fitting_component]] = disp_lsq_error

	#write lsq_error output file
	with open("component_lsq_error.txt", 'w') as f:
		for component in energyComponents:
			f.write(component + ":" + '\t')
			f.write(str(component_lsq_error[component]) + '\n')

	print("===========================================================================")
	print(" Optimized parameters successfully written to optimized_params.constraints ")
	print("===========================================================================")

	return res

def main():
	try:
		useInitialParamsFile(sys.argv[1])
		print("Importing given parameters")
		print("Converting given initial parameter dictionary to list")

		allAtomtypes.sort() #alphabetize atomtypes list
		minimizeDictionaryToList()
		print("Initial parameter list successfully created")
		#TODO: get names of all dimers in library
		#TODO: perform opt.

	except IndexError:
		print("Importing all dimer parameters from library")
		importdimecularParameters()
		print("Calculating average parameter values")
		getInitalAverageParameters()

		allAtomtypes.sort() #alphabetize atomtypes list so everything iterizes over it the same way
		print("Converting initial parameter dictionary to list") #must be list for minimize function
		averageParamsDictToList()
		print("Initial parameter list successfully created")
		print("Making output directories")
		make_output_directories()
		print("Beginning parameter generalization")
		optimzeGeneralParameters()

main()
#######################################
# Known Bugs:

#######################################

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

#import settings file
import settings as generalizer_settings

# Global variables

scriptDir = os.getcwd()

#set variable to fit bii in metaPOInter
fit_b = generalizer_settings.fit_b

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

#list of all atomtypes from all imported molecules
allAtomtypes = []

#stores atomtypes held by each molecule
molecular_atomtypes = {}

#list of aniso atomtypes
anisoAtomtypes = []

#contains all imported molecular parameters
molecules = {}

#store a list of spherical harmonics for different atomtypes
atomtypeSphericalHarmonics = {}

#track how many times metaPOInter gets called
minimizeIterCount = 1

#keeps track of what component is being fit
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

#stores POInter objects for each molecule
molecule_POInter_objects = {}

#stores final parameters by component
optimized_params = {}

#store qm energies for each molecule
qm_energies = {}

#stores initial B parameters for each atomtype
atomtype_B_params = {}

#stores C parameters for each atomtype
atomtype_C_params = {}

#stores drude charge for each atomtype
atomtype_drude_charges = {}

#flag for if drude calculation has already been performed for a given molecule
drude_flags = {}

#collects total lsq_error for each component
component_lsq_error = {component : 0.0 for component in energyComponents}

##############################################################################
# Class for holding the parameters from a single ab initio FF from a JSON file
##############################################################################

class MolecularParameters:
	#creates object to hold parameters for an individual molecular ab initio FF

	def __init__(self, jsonInput, molecule):
		self.jsonInput = jsonInput
		self.moleculeName = molecule
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
		self.atomtypesList = self.jsonInput.keys()
		#add list of atomtypes to entry in molecular_atomtypes and sort
		molecular_atomtypes[self.moleculeName] = [atomtype for atomtype in self.atomtypesList]
		molecular_atomtypes[self.moleculeName].sort()
		#populate molecule with atomtype objects
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
		#creates object to store parameters for an atomtype within a molecule
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
	"""import all parameters from JSON files. Finds all .constraints files in directory of initial constrait files, converts JSON to dictionary and stores in global dictionary of initial parameters with molecule name as key

	Parameters
	----------
	None

	Returns 
	-------
	None

	"""
	try:
		os.chdir(scriptDir + '/abInitioConstraints')
		#get list of all constraint (ab initio) files
		file_list = glob.glob("*.constraints")
		file_list.sort()
		
		global inputJsons
		inputJsons = {}
		for molecule in file_list:
			molJson = readJson(molecule)
			#get name of molecule
			i = 8
			while molecule[i] != "_":
				i += 1
			molName = molecule[8:i] 
			inputJsons[molName] = molJson

		os.chdir(scriptDir)

	except OSError:
		print "need to create ab initio constraints directory!"

##################################################################################################################
# functions for use with imported inital parameters
# to be used when we already have parameters we want to start minimize() with, do not need to average from library
##################################################################################################################

def useInitialParamsFile(molName):
	"""check to see if .params file already exists, import initial values if so

	Parameters
	----------
	molName: name of molecule to find .params file for

	Returns 
	-------
	None, but stores given parameters in initialMinimizeValues

	"""
	if os.path.isfile(molName + ".params"):
		global initialMinimizeValues
		initialMinimizeValues = MolecularParameters(readJson(molName + ".params"), molName) #need to make this a list to interface with minimize

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
#fuctions to use with library of initial parameters
######################################################

def importMolecularParameters():
	"""import all ab initio parameters

	Parameters
	----------
	None

	Returns 
	-------
	None, but stores atomtypes of imported molecules in global list of atomtypes

	"""
	importParameters()
	for mol in inputJsons.keys():
		print "Importing molecular parameters: " + mol
		molecules[mol] = MolecularParameters(inputJsons[mol], mol)
		molAtomtypes = molecules[mol].atomtypesList
		print mol + " parameters successfully imported!"

def averageAvalues():
	"""return averages of all A values in dictionary

	Parameters
	----------
	None

	Returns 
	-------
	Dictionary of average A values for all atomtypes, averaged from all previously imported molecules. Nested dictionary of {atomtype {energy component : value}} 

	"""
	allAvalsRaw = {}
	allAvalsOut = {}
	for molecule in molecules:
		currentMolecule = molecules[molecule]
		for atomtype in allAtomtypes:
			if atomtype in allAvalsRaw: #if we already have this atomtype, add parameters to existing list
				for component in energyComponents:
					try: #try/except allows for multiple atomtypes in different input files
						allAvalsRaw[atomtype][component].append(currentMolecule[atomtype]["A"][component])
					except KeyError:
						pass
			else:
				componentsDict = {}
				for component in energyComponents:
					try:
						componentsDict[component] = [currentMolecule[atomtype]["A"][component]]
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
	Dictionary of average B values for all atomtypes, averaged from all previously imported molecules.

	"""
	allBvalsRaw = {}
	allBvalsOut = {}
	for molecule in molecules:
		currentMolecule = molecules[molecule]
		for atomtype in allAtomtypes:
			if atomtype in allBvalsRaw: #if we already have this atomtype, add parameters to existing list
				try:
					allBvalsRaw[atomtype].append(currentMolecule[atomtype]["B"])
				except KeyError:
					pass
			else:
				try:
					allBvalsRaw[atomtype] = [currentMolecule[atomtype]["B"]]
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
	Dictionary of average C values for all atomtypes, averaged from all previously imported molecules. Nested dictionary of {atomtype {order : value}} 

	"""
	allCvalsRaw = {}
	allCvalsOut = {}
	for molecule in molecules:
		currentMolecule = molecules[molecule]
		for atomtype in allAtomtypes:
			if atomtype in allCvalsRaw: #if we already have this atomtype, add parameters to existing list
				for term in dispersionTerms:
					try:
						allCvalsRaw[atomtype][term].append(currentMolecule[atomtype]["C"][term])
					except KeyError:
						pass
			else: #if we have not seen this atomtype yet
				termsDict = {}
				for term in dispersionTerms:
					try:
						termsDict[term] = [currentMolecule[atomtype]["C"][term]]
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
	Dictionary of average Aniso values for all atomtypes, averaged from all previously imported molecules. Nested dictionary of {atomtype {energy component {spherical harmonic : value}}} 

	"""
	#make dictionaries to store parameters
	allAnisoValsRaw = {}
	allAnisoValsOut = {}

	#get parameters for each ab initio molecule
	for molecule in molecules:
		currentMolecule = molecules[molecule]
		#for each atomtype, check if already in dictionary
		for atomtype in allAtomtypes:
			try:
				currentAtomtype = currentMolecule[atomtype]
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

##################################################################################
# metaPOInter and helper functions
# metaPOInter interfaces with POInter to get lsq_error and dlsq_error for each parameter set
##################################################################################

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
	i = 0 #keep track of position in parameterList

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
				bounds_list.append((-1e0,1e0))
				i += 1

	return bounds_list

def get_fit_parameters(all_params, fit_molecule):
	"""
	Fetches the parameters to be fit for a given component.

	Parameters
	----------

	all_params: list, parameters for all atomtypes for the given component, in form [atomtype1 A, B, aniso; atomtype2 A, B, aniso...]
	fit_molecule: string, name of molecule being worked with

	Returns
	-------

	fit_params: list, params (for the current component) to be input to POInter in format [atomtype1 A, aniso, B; atomtype2 A, aniso, B], only including atomtypes for that molecule

	"""

	molecule_atomtypes = molecular_atomtypes[fit_molecule]

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
		if fitting_component == 0:
			n_params += 1 #if we are including the B param (exchange)

		#if our molecule has this atomtype, add these parameters to the output list	
		if pot_atomtype in molecule_atomtypes:
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
	print all_params

	for atomtype in allAtomtypes:
		n_params = atomtype_nparams[atomtype] - 4 #subtract out C params
		#get A param
		output_params.append(all_params[input_index + fitting_component])
		input_index += 5
		#get aniso params
		if atomtype in anisoAtomtypes:
			init_input_index = input_index
			n_spherical_harmonics = atomtypeSphericalHarmonics[atomtype]
			for i in range(fitting_component):
				input_index += n_spherical_harmonics
			output_params.append(all_params[input_index : input_index + n_spherical_harmonics])
			input_index = init_input_index
		#get B if required (exchange)
		if fitting_component == 0:
			output_params.append(1.0)

		#move index counter to start of next atomtype
		input_index -= 5
		input_index += n_params

	return output_params

def add_dlsq_error(all_dlsq_error, new_dlsq_error, moleculeName):
	""" add dlsq_error returned by POInter to list of all dlsq_error

	parameters
	----------
	all_dlsq_error: list, all dlsq_error for all atomtypes and parameters in format corresponding to [atomtype1 A, aniso, B; atomtype2 A, aniso, B]
	new_dlsq_error: list, dlsq_error returned from POInter to be added to all_dlsq_error (in same format at all_dlsq_error, just without some atomtypes not found in that molecule)
	moleculeName: name of molecule

	returns
	-------
	all_dlsq_error updated with new_dlsq_error added
	"""

	global fitting_component
	mol_atomtypes = molecular_atomtypes[moleculeName]
	all_dlsq_error_index = 0

	#add new dlsq_error to all_dlsq_error
	new_dlsq_error_index = 0

	for pot_atomtype in allAtomtypes:
		n_params = atomtype_nparams[pot_atomtype] - 4 # excluding C params which are not fit
		if fitting_component != 0:
			n_params -= 1 # subtract B param if not fitting exchange
		n_params = n_params / 5 # to find numer of parameters for a single component energy
		if pot_atomtype in mol_atomtypes:
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
			n_spherical_harmonics = atomtypeSphericalHarmonics[atomtype]
		else:
			n_spherical_harmonics = 0

		shift = 1 + n_spherical_harmonics
		#get A params	
		for component in energyComponents:
			if component == "Exchange":
				b_scale = parameterDict[component][input_atomtype_i + n_spherical_harmonics + exchange_shift]
				exchange_shift += 1
			a_param = parameterDict[component][input_atomtype_i]
			outputList.append(a_param)
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
			n_spherical_harmonics = atomtypeSphericalHarmonics[atomtype]
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
	total_lsq_error: float, sum of all lsq_error for each individual molecular dispersion component
	"""
	#set fitting component for dispersion
	global fitting_component
	fitting_component = 4

	total_lsq_error = 0.0

	for molecule in molecules.keys():
		#get POInter object and setting fitting_component
		pointerModel = molecule_POInter_objects[molecule]
		pointerModel.component = fitting_component
		#get qm energies
		pointerModel.qm_fit_energy = np.array(pointerModel.subtract_hard_constraint_energy())
		#set POInter instance variables
		pointerModel.fit_bii = False
		pointerModel.n_isotropic_params = 1
		n_aiso = pointerModel.n_isotropic_params if not pointerModel.fit_bii else pointerModel.n_isotropic_params - 1
		n_aaniso = n_aiso
		#calc lsq_error, dlsq_error
		molecule_atomtypes = molecular_atomtypes[molecule]
		currentParams = tuple(1.0 for atomtype in molecule_atomtypes)
		print "Fitting component", fitting_component, "with current parameters", currentParams
		pointerModel.get_num_eij = pointerModel.generate_num_eij(currentParams)
		pointerOutput = pointerModel.calc_leastsq_ff_fit(currentParams)
		total_lsq_error += pointerOutput[0]

	return total_lsq_error

def metaPOInter(parameterList):
	"""
	Interface function with POInter code. Called by minimize(), will get lsq_error and dlsq_error for each molecule in library, returns total lsq_error and total dlsq_error

	Parameters
	----------

	None

	Returns
	-------

	total_lsq_error: float, sum lsq_error for entire library of molecules

	dlsq_error: list, list of total dlsq_error's for each parameter being optimized

	"""

	total_lsq_error = 0
	dlsq_error = [0.0 for i in range(len(parameterList))]
	global fitting_component

	#get lsq_error and dlsq_error for each molecule from POInter, add to running totals
	os.chdir(scriptDir + "/abInitioSAPT")
	for molecule in molecules.keys():
		if molecule not in molecule_POInter_objects:
			saptFile = molecule + "_" + molecule + ".sapt"
			#make POInter object
			pointerModel = Pointer.FitFFParameters(fit = False, energy_file = saptFile)

			#call required POInter start-up functions
			kwargs = {"mon1": molecule, "mon2": molecule}
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
			molecule_POInter_objects[molecule] = pointerModel
			#set molecule drude flag to false
			drude_flags[molecule] = False
		else:
			pointerModel = molecule_POInter_objects[molecule]
			pointerModel.component = fitting_component

		#get drude energy if component 1
		if fitting_component == 1 and not drude_flags[molecule]:
			pointerModel.get_drude_oscillator_energy()
			drude_flags[molecule] = True
		#subtract hard constraints energy from energy to be fit
		mol_qm_dict = qm_energies[molecule]
		if fitting_component not in mol_qm_dict:
			pointerModel.qm_fit_energy = np.array(pointerModel.subtract_hard_constraint_energy())
			qm_energies[molecule][fitting_component] = pointerModel.qm_fit_energy
		else:
			pointerModel.qm_fit_energy = qm_energies[molecule][fitting_component]

		#set fitting component and B coeff fitting (if exchange component)
		if fitting_component == 0:
			pointerModel.fit_bii = True
			pointerModel.n_isotropic_params = 2

			#set POInter instance variables required for fitting B param

			i_bparams = []

			param_i = 0
			for atomtype in molecular_atomtypes[molecule]:
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
		currentParams = tuple(get_fit_parameters(parameterList, molecule))
		print "Fitting component", fitting_component, "with current parameters", currentParams
		pointerModel.get_num_eij = pointerModel.generate_num_eij(currentParams)
		pointerOutput = pointerModel.calc_leastsq_ff_fit(currentParams)

		print
		print "leastsq error: ", pointerOutput[0]
		print

		#add lsq_error, dlsq_error to totals

		total_lsq_error += pointerOutput[0]
		dlsq_error_components = pointerOutput[1]

		dlsq_error = add_dlsq_error(dlsq_error, dlsq_error_components, molecule)	

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
	#initialize qm energies dictionaries for each molecule
	global qm_energies
	qm_energies = {mol : {} for mol in molecules}

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
		if fitting_component == 0:
			bnds = all_bnds
		else:
			bnds = a_bnds

		#make list of parameters for component to pass down to POInter
		#get only parameters for the component we are fitting

		#TODO: make so user can choose whether or not to fit B params
		#TODO: make sure input B params are the same for every molecule
			#TODO: use same .exp file ()
		metaPOInter_input = get_component_parameters(metaPOInter_initial_list)

		res = minimize(metaPOInter,metaPOInter_input,method='L-BFGS-B',\
							jac=True,\
							options={'disp':True,'gtol':pgtol,'ftol':ftol,'maxiter':maxiter},\
							bounds = bnds)

		popt = res.x
		print "Optimized parameters for component " + str(fitting_component) + " : " + str(popt)

		if not res.success:
			print 'Warning! Optimizer did not terminate successfully, and quit with the following error message:'
			print
			print res.message
			print
			return

		global optimized_params
		optimized_params[energyComponents[fitting_component]] = popt

		#reset new B params in POInter objects after fitting exchange component
		if fitting_component == 0:
			updated_B_params = get_updated_B_from_exchange(popt)
			for mol in molecules.keys():
				i = 0
				for atomtype in allAtomtypes:
					pointerModel = molecule_POInter_objects[mol]
					pointerModel.params[atomtype][0]['B'] = updated_B_params[i]
					i += 1
		
		fitting_component += 1

	#write output JSON file with optimized parameters
	optimized_params["Dispersion"] = [1.0 for atomtype in allAtomtypes]
	final_params_list = write_output_params_to_list(optimized_params)
	map_params(final_params_list, "optimized_params")

	#calc dispersion lsq_error and store
	disp_lsq_error = calc_disp_lsq_error()
	component_lsq_error[energyComponents[fitting_component]] = disp_lsq_error

	#write lsq_error output file
	with open("component_lsq_error.txt", 'w') as f:
		for component in energyComponents:
			f.write(component + ":" + '\t')
			f.write(str(component_lsq_error[component]) + '\n')

	print "==========================================================================="
	print " Optimized parameters successfully written to optimized_params.constraints "
	print "==========================================================================="

	return res

def main():
	try:
		useInitialParamsFile(sys.argv[1])
		print "Importing given parameters"
		print "Converting given initial parameter dictionary to list"

		allAtomtypes.sort() #alphabetize atomtypes list
		minimizeDictionaryToList()
		print "Initial parameter list successfully created"
		#TODO: get names of all molecules in library
		#TODO: perform opt.

	except IndexError:
		print "Importing all molecular parameters from library"
		importMolecularParameters()
		print "Calculating average parameter values"
		getInitalAverageParameters()

		allAtomtypes.sort() #alphabetize atomtypes list
		print "Converting initial parameter dictionary to list"
		averageParamsDictToList()
		print "Initial parameter list successfully created"
		optimzeGeneralParameters()

main()
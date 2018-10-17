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

# Local modules
import functional_forms 
from multipoles import Multipoles
import fit_ff_parameters as Pointer

# Global variables

scriptDir = os.getcwd()

#variable used so store imported JSON dictionaries when working with library of ab initio FF's
inputJsons = {}

#list of different types of parameters we care about
parameterTypes = ['A', 'B', 'C', 'aniso']

#names of energy components
energyTerms = ['Exchange', 'Electrostatics', 'Induction', 'Dhf', 'Dispersion']

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

#list of aniso atomtypes
anisoAtomtypes = []

#contains all imported molecular parameters
molecules = {}

#store a list of spherical harmonics for different atomtypes
atomtypeSphericalHarmonics = {}

#track how many times metaPOInter gets called
minimizeIterCount = 1


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
		for atomtype in self.atomtypesList:
			self.atomtypes[str(atomtype)] = self.Atomtype(self.jsonInput[atomtype], self.setAnisoFlag(self.jsonInput[atomtype]))
			global atomtypeSphericalHarmonics
			atomtypeSphericalHarmonics[str(atomtype)] = self.atomtypes[str(atomtype)].getSphericalHarmonics()
			if str(atomtype) not in allAtomtypes: #add atomtype to list of all atomtypes if not already there
				allAtomtypes.append(str(atomtype))
			if self.anisoFlag and atomtype not in anisoAtomtypes:
				anisoAtomtypes.append(str(atomtype))

	class Atomtype:
		#creates object to store parameters for an atomtype within a molecule
		def __init__(self, inputDict, anisoFlag = False):
			self.parameters = {}
			self.parameters["A"] = self.Parameters(inputDict['A'], 'A')
			self.parameters["B"] = inputDict['B'] #since only one B value for each atomtype
			self.parameters["C"] = self.Parameters(inputDict['C'], 'C')
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
				for i in range(len(energyTerms)):
					innerOutput = {}
					for j in range(len(sphericalHarmonics)):
						innerOutput[str(sphericalHarmonics[j])] = inputList[i][j]
					output[energyTerms[i]] = innerOutput
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

	allAtomtypes.sort() #alphabetize atomtypes list

	#loop over all atomtypes
	for atomtype in allAtomtypes:
		currentAtomtype = initialMinimizeValues[atomtype]
		#loop over A energy components
		aParams = currentAtomtype["A"]
		for component in energyTerms:
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
		for component in energyTerms:
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
				for term in energyTerms:
					try: #try/except allows for multiple atomtypes in different input files
						allAvalsRaw[atomtype][term].append(currentMolecule[atomtype]["A"][term])
					except KeyError:
						pass
			else:
				termsDict = {}
				for term in energyTerms:
					try:
						termsDict[term] = [currentMolecule[atomtype]["A"][term]]
						allAvalsRaw[atomtype] = termsDict
					except KeyError:
						pass
	for atomtype in allAvalsRaw: #calculate averages from list for each atomtyp for each energy type
		termAvgs = {}
		for term in energyTerms:
			termAvgs[term] = sum(allAvalsRaw[atomtype][term])/len(allAvalsRaw[atomtype][term])
		allAvalsOut[atomtype] = termAvgs
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
						for energyTerm in energyTerms:
							currentValues = currentAnisoParams[energyTerm]
							for sh in currentSphericalHarmonics:
								allAnisoValsRaw[atomtype][energyTerm][sh].append(currentValues[sh])
					else:
						energyTermValues = {}
						for energyTerm in energyTerms:
							currentValues = currentAnisoParams[energyTerm]
							shValues = {} #store values for spherical harmonics for a given energy term
							for sh in currentSphericalHarmonics:
								shValues[sh] = [currentValues[sh]]
							energyTermValues[energyTerm] = shValues
						allAnisoValsRaw[atomtype] = energyTermValues
			except KeyError:
				pass
	#perform averaging
	for atomtype in allAnisoValsRaw:
		termsAvg = {}
		for term in energyTerms:
			energyTermsAvgs = {}
			currentValues = allAnisoValsRaw[atomtype][term]
			for sh in currentValues:
				energyTermsAvgs[str(sh)] = sum(currentValues[sh])/len(currentValues[sh])
			termsAvg[str(term)] = energyTermsAvgs
		allAnisoValsOut[atomtype] = termsAvg
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
		#loop over A components
		aParams = initialAverageValues["A"][atomtype]
		for component in energyTerms:
			output.append(aParams[component])
		#B
		output.append(initialAverageValues["B"][atomtype])
		#loop over C components
		cParams = initialAverageValues["C"][atomtype]
		for component in dispersionTerms:
			output.append(cParams[component])
		#loop over aniso components
		try:
			anisoParams = initialAverageValues["aniso"][atomtype]
			sphericalHarmonics = atomtypeSphericalHarmonics[atomtype]
			for component in energyTerms:
				#loop over spherical harmonics
				for sh in sphericalHarmonics:
					output.append(anisoParams[component][sh])
		except KeyError: #catch if atomtype is not anisotropic
			pass

	global initialValuesList
	initialValuesList = initialValuesList + output

##################################################################################
# metaPOInter and helper functions
# metaPOInter interfaces with POInter to get RMSE and dRMSE for each parameter set
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

	allAtomtypes.sort() #make sure atomtypes are alphabetical order

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
		#check if atomtype is anisotropic. If so, fetch spherical harmonicsand save in atomtypeParams. Also save aniso A values, with empty lists if not anisotropic
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
		if atomtype == "O(0)":
			atomtypeParams["drude_charge"] = -sqrt(17.174511/0.1) #TODO: change to actual values
		else:
			atomtypeParams["drude_charge"] = 0.0 # TODO: change to actual values
		output[atomtype] = atomtypeParams

	#write .constraints file with parameters from list
	with open(fileName + ".constraints", 'w') as f:
		json.dump(output, f, indent = 4)

def calc_harmonic_constraint_error(params, k = 1e-2):
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
        
    dharmonic_error : list
            Derivatives of the harmonic error with respect to each
            parameter.
    """
	harmonic_error = 0
	dharmonic_error = []
	dharmonic_error = [ 0 for _ in params]
	for i in range(len(params)):
		b = params[i]
		b0 = initialValuesList[i]
		harmonic_error += k*(b - b0)**2
		dharmonic_error[i] = 2*k*(b - b0)*b0

	return harmonic_error, dharmonic_error


def metaPOInter(parameterList):
	"""
	Interface function with POInter code. Called by minimize(), will get RMSE and dRMSE for each molecule in library, returns total RMSE and total dRMSE

	Parameters
	----------

	None

	Returns
	-------

	totalRMSE: float, sum RMSE for entire library of molecules

	dRMSE: list, list of total dRMSE's for each parameter being optimized

	"""
	global minimizeIterCount
	print
	print "Minimization iteration: " + str(minimizeIterCount) 
	print

	totalRMSE = 0
	dRMSE = []

	#make temp.constraints with current parameters
	map_params(parameterList, "temp")

	#make tuple of parameters
	parameterTuple = tuple(parameterList)

	#get RMSE and dRMSE for each molecule from POInter, add to running totals
	os.chdir(scriptDir + "/abInitioSAPT")
	for molecule in molecules.keys():
		saptFile = molecule + "_" + molecule + ".sapt"
		#make POInter object with temp.constraints from map_params as input with .sapt file
		pointerModel = Pointer.FitFFParameters(fit = False, energy_file = saptFile, param_file = "temp.constraints")
		#call required POInter start-up functions
		kwargs = {"mon1": molecule, "mon2": molecule}
		pointerModel.read_settings(['default'], kwargs)
		pointerModel.read_energies()
		pointerModel.read_params()
		pointerModel.initialize_parameters()
		#set required POInter instance variables
		pointerModel.component = 0
		pointerModel.n_isotropic_params = pointerModel.default_n_isotropic_params
		pointerModel.get_num_eij = pointerModel.generate_num_eij(parameterTuple)
		pointerModel.final_energy_call = True
		pointerModel.qm_fit_energy = np.array(pointerModel.subtract_hard_constraint_energy())
		#get RMSE and dRMSE, including harmonic penalty
		pointerOutput = pointerModel.calc_leastsq_ff_fit(parameterTuple)
		harmonic_errors = calc_harmonic_constraint_error(parameterTuple)
		totalRMSE += pointerOutput[0]
		totalRMSE += harmonic_errors[0]
		if dRMSE == []:
			dRMSE = pointerOutput[1].tolist()
			for i in range(len(harmonic_errors[1])):
				dRMSE[i] += harmonic_errors[1][i]
		else:
			for i in range(len(dRMSE)):
				dRMSE[i] += pointerOutput[1].tolist()[i]
				dRMSE[i] += harmonic_errors[1][i]

	minimizeIterCount += 1
	os.chdir(scriptDir)

	return totalRMSE, np.asarray(dRMSE)

def optimzeGeneralParameters():
	"""
	Primary function of parameterGeneralizer module. Will call scipy minimize to find optimal parameters for library.

	Parameters
	----------

	None

	Returns
	-------
	Global RMSE for library
	"""

	#set parameters for minimize function
	maxiter = 5000
	pgtol = -1e-17
	ftol = 1e-17

	map_params(initialValuesList,"initial_params_test")

	res = minimize(metaPOInter,initialValuesList,method='L-BFGS-B',\
						jac=True,\
						options={'disp':True,'gtol':pgtol,'ftol':ftol,'maxiter':maxiter})
	popt = res.x
	success = res.success
	message = res.message

	if not res.success:
		print 'Warning! Optimizer did not terminate successfully, and quit with the following error message:'
		print
		print res.message
		print
		return

	#write output JSON file with optimized parameters
	map_params(popt, "optimized_params")

	#remove temp.constraints file
	os.remove("temp.constraints")

	print "==========================================================================="
	print " Optimized parameters successfully written to optimized_params.constraints "
	print "==========================================================================="

	return res

def main():
	try:
		useInitialParamsFile(sys.argv[1])
		print "Importing given parameters"
		print "Converting given initial parameter dictionary to list"
		minimizeDictionaryToList()
		print "Initial parameter list successfully created"
		#TODO: get names of all molecules in library
		#TODO: perform opt.

	except IndexError:
		print "Importing all molecular parameters from library"
		importMolecularParameters()
		print "Calculating average parameter values"
		getInitalAverageParameters()
		print "Converting initial parameter dictionary to list"
		averageParamsDictToList()
		print "Initial parameter list successfully created"
		optimzeGeneralParameters()

main()
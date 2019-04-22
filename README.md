# Parameter Generalizer for the MASTIFF Intermolecular Force Field
Written by Althea Hansel, Harvey Mudd College '19


<img src="http://www.clipartpal.com/_thumbs/pd/animal/dog/pointer_3.png" width="300" >

Introduction
------------
The MASTIFF parameter generalizer is an open-source python script for optimizing general force field parameters for the MASTIFF intermolecular force field. 
It has primarily been designed to optimize parameters by minimizing least-squared-error summed across a library of input DFT-SAPT dimer calculations.

The parameter generalizer can be called directly from the command line.


Getting Started
---------------
The parameter generalizer is a very new code. 
Email thansel at hmc dot edu
with any questions.

Documentation can be found in this repository's wiki.

The examples directory gives an example of the inputs for making general parameters fit to ethane and propane, including the required settings.py file.

Issues with the code (bugs, unclear documentation, suggestions for improvement)
should be reported directly to Althea Hansel at the email listed above.


Dependencies
------------
To run properly, POInter requires the following python packages:

* Numpy
* Scipy
* Sympy
* [numexpr](https://github.com/pydata/numexpr)
* [dill](https://github.com/uqfoundation/dill) (optional, but speeds up multipole calculations)
* [POInter] (https://git.chem.wisc.edu/schmidt/force_fields)

Downloads
---------
The MASTIFF parameter generalizer is free software.


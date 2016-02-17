#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file main.py
# This module defines the function which communicates the GAMESS output data
# to Libra and vice versa.
# It outputs the files needed for excited electron dynamics simulation.


from gamess_to_libra import *
from nve import *
from create_gamess_input import *
from exe_gamess import *

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
#cwd = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code"
#print "Using the Libra installation at", cwd
#sys.path.insert(1,cwd+"/_build/src/mmath")
#sys.path.insert(1,cwd+"/_build/src/qchem")

#print "\nTest 1: Importing the library and its content"
#from libmmath import *
#from libqchem import *

def initial_gamess_exe(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : the input data from "submit_templ.slm", in the form of dictionary
    # This function executes GAMESS program and then extracts data from GAMESS output files
    # (mainly atomic orbitals, molecular energies, molecular coefficients, gradients).
    #
    # Used in:  main.py/main

    GMS_DIR = params["GMS_DIR"]
    GMS_JOB = params["GMS_JOB"]
    params["l_gam_for"] = keep_GMS_INP_format(GMS_DIR,GMS_JOB)

    exe_gamess(params,GMS_JOB)

    # 1-st file - time "t"  old
    ao, E, C, Grad, data = unpack_file(GMS_DIR, GMS_JOB)

    return ao, E, C, Grad, data

def main(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : the input data from "submit_templ.slm", in the form of dictionary
    # This function prepares initial parameters from GAMESS output file
    # and executes classical MD in Libra and Electronic Structure Calculation in GAMESS 
    # iteratively.
    #
    # Used in:  main.py

    ao1, E1, C1, Grad, data = initial_gamess_exe(params)

    print "data= ",data

    # main routine

    nve(data,params,ao1,E1,C1,Grad)

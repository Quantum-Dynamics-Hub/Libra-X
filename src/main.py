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
import sys

# Path the the source code
sys.path.insert(1,"/user/alexeyak/Programming/libra-gamess_interface/src")

cwd = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code"
print "Using the Libra installation at", cwd
sys.path.insert(1,cwd+"/_build/src/mmath")
sys.path.insert(1,cwd+"/_build/src/qchem")
sys.path.insert(1,cwd+"/_build/src/dyn")
sys.path.insert(1,cwd+"/_build/src/chemobjects")
sys.path.insert(1,cwd+"/_build/src/hamiltonian")


import os
import sys
import math
from gamess_to_libra import *
from nve import *
from create_gamess_input import *


def main(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : the input data from "submit_templ.slm", in the form of dictionary
    # This function prepares initial parameters from GAMESS output file
    # and executes classical MD in Libra and Electronic Structure Calculation in GAMESS 
    # iteratively.
    #
    # Used in:  main.py

    ################# Step 0: Use the initial file to create a working input file ###############
 
    os.system("cp %s %s" %(params["gms_inp0"], params["gms_inp"]))

    ################# Step 1: Read initial input and run first GMS calculation ##################    
    
    params["gms_inp_templ"] = read_gms_inp_templ(params["gms_inp"])

    exe_gamess(params)

    ao, E, C, Grad, data = unpack_file(params["gms_out"])

    print data

    ################## Step 2: Initialize molecular system and run MD ###########################

    print "Initializing system..."
    syst = init_system(data, Grad)

    print "Starting MD..."
    run_MD(syst,ao,E,C,data,params)

  

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

import os
import sys
import math


# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["src_path"]) # Path the the source code
sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_qchem_path"])
sys.path.insert(1,os.environ["libra_dyn_path"])
sys.path.insert(1,os.environ["libra_chemobjects_path"])
sys.path.insert(1,os.environ["libra_hamiltonian_path"])

from gamess_to_libra import *
from md import *
from create_gamess_input import *


def main(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params  the input data from "submit_templ.slm", in the form of dictionary
    # \param[out] test_data  the output data for debugging, in the form of dictionary
    # \param[out] data  the data extracted from gamess output file, in the form of dictionary
    # This function prepares initial parameters from GAMESS output file
    # and executes classical MD in Libra and Electronic Structure Calculation in GAMESS 
    # iteratively.
    # Parallelly, it executes TD-SE calculation for simulating excited eletronic dynamics.
    #
    # Used in:  main.py

    ################# Step 0: Use the initial file to create a working input file ###############
 
    os.system("cp %s %s" %(params["gms_inp0"], params["gms_inp"]))

    ################# Step 1: Read initial input and run first GMS calculation ##################    
    
    params["gms_inp_templ"] = read_gms_inp_templ(params["gms_inp"])

    exe_gamess(params)

    ao, E, C, Grad, data = unpack_file(params["gms_out"],params["debug_gms_unpack"])

    ################## Step 2: Initialize molecular system and run MD with TD-SE ####

    print "Initializing system..."
    syst = []
    # store several initial nuclei systems with different momenta
    for i in xrange(params["nconfig"]):
        syst.append(init_system(data, Grad,params["Temperature"]))
    
    #print "Initializing electronic variables"    
    #el = []
    #nstates = len(params["excitations"])
    #for i_ex in xrange(nstates):  # loop over all initial excitations
    #    eltmp = Electronic(nstates,i_ex)
    #    el.append(eltmp)

    print "Starting MD..."
    test_data = run_MD(syst,ao,E,C,data,params)

    return data, test_data

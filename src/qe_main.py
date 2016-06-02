#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file main.py
# This module defines the function which communicates QUANTUM ESPRESSO output data
# to Libra and vice versa.
# It outputs the files needed for excited electron dynamics simulation.
import os
import sys
import math


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *



#Import libraries
from read_qe_inp_templ import*
from exe_espresso import*
from create_qe_input import*
from unpack_file import*
from md import *
from export_wfc import *


def main(params):
##
# Finds the keywords and their patterns and extracts the parameters
# \param[in] params : the input data from "submit_templ.slm", in the form of dictionary
# This function prepares initial parameters from QUANTUM ESPRESSO output file
# and executes classical MD in Libra and Electronic Structure Calculation in QUANTUM ESPRESSO
# iteratively.
#
# Used in:  main.py

    ################# Step 0: Use the initial file to create a working input file ###############
    os.system("cp %s %s" %(params["qe_inp00"], params["qe_inp0"]))

    ################# Step 1: Read initial input and run first QS calculation ##################    

    params["qe_inp_templ"] = read_qe_inp_templ(params["qe_inp0"])

    exe_espresso(params["qe_inp0"], params["qe_out0"])
    tot_ene, label, R, grad,params["norb"],params["nel"],params["nat"],params["alat"] = unpack_file(params["qe_out0"], params["qe_debug_print"],1)

    ################## Step 2: Initialize molecular system and run MD ###########################

    print "Initializing system..."
    df = 0 # debug flag
    #Generate random number
    rnd = Random()

    # Here we use libra_py module!
    syst = init_system.init_system(label, R, grad, rnd, params["Temperature"], params["sigma_pos"], df, "elements.txt")      

    # Create a variable that will contain propagated nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    
    n_el = params["nel"]
    n_mo = params["num_MO"]
    wfc = {}  # wavefunction dictionary, where all the coefficients of the MO basis will be saved
    # Running SCF calculation for different excited states, extracting their Energies, Forces and wavefucntion coefficients
    # savings coefficients as coeff_old
    for i in xrange(len(params["excitations"])):
        write_qe_input(params["qe_inp%i" %i],label,mol,params["excitations"][i],params)
        exe_espresso(params["qe_inp%i" % i], params["qe_out%i" % i] )
        wfc["coeff_old_%i"%i] = read_qe_wfc("x%i.export/wfc.1"%i, "Kpoint.1", n_el, n_mo)



    # starting MD calculation
    test_data = run_MD(label,syst,params,wfc)
    return test_data


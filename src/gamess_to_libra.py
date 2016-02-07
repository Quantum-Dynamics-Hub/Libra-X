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

## \file gamess_to_libra.py 
# This module implements the functions that extract parameters from the gamess output file:
# atomic forces , eigenenergies, eigenfunctions, and atomic basis information.
# The forces are used for simulating Classical MD on Libra.
#

# Need to move the descriptions below somewhere else:
#* The eigenfunctions and the atomic basis sets are communicated to Libra modules
#* to calculate the overlap matrix of the atomic orbitals and eigenfunctions.
#* From the overlap matrix of eigenfunctions, We get Non-Adiabatic Couplings(NAC).
#* Eigenenergies and NAC are used for simulating excited electron dynamics.
#**********************************************************************************

from detect import *
from extract import *
from ao_basis import *
from overlap import *
from Ene_NAC import *

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

def unpack_file(filename,runtype,basis_opt):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] filename The name of the (GAMESS output) file from which we will be getting data
    # \param[in] runtype The option controlling how to interpret the file(single point or optimization)
    # \param[in] basis_op The option controlling assumed orthogonality of basis (as in semiempirics)
    # This function returns the data extracted from the file, in the form of dictionary
    #
    # Used in:  gamess_to_libra.py/gamess_to_libra

    f_gam = open(filename,"r")
    l_gam = f_gam.readlines()
    f_gam.close()


    data = {}

    # detect the columns showing parameters
    detect(l_gam,data,runtype)    

    # extract the parameters from the columns detected
    extract(l_gam,data)

    # Construct the AO basis
    data["ao_basis"] = ao_basis(data,basis_opt) # the construction of the AO basis should not
                                                # depend on basis_opt

    return data["ao_basis"], data["E"], data["C"], data["gradient"], data

    

def gamess_to_libra(par):
    ## 
    # Extracts data (coordinates, forces, MOs, and basis info) GAMESS output file
    # and create AO basis
    # \param[in] par The dictionary of all necessary parameters
    # 
    # Used in: main.py/


    # 1-st file - time "t" 
    ao1, E1, C1, Grad1, data1 = unpack_file(par["gamess_out1"],par["runtype"],par["basis_option"])

    # 2-nd file - time "t+dt"
    ao2, E2, C2, Grad2, data2 = unpack_file(par["gamess_out2"],par["runtype"],par["basis_option"])



    # calculate overlap matrix of atomic orbitals and eigenfunctions
    P11, P22, P12, P21 = overlap(ao1,ao2,C1,C2,par["basis_option"])

    print "P11 and P22 matrixes should show orthogonality"
    print "P11 is";    P11.show_matrix()
    print "P22 is";    P22.show_matrix()

    print "P12 and P21 matrixes show overlap of MOs for different molecular geometries "
    print "P12 is";    P12.show_matrix()
    print "P21 is";    P21.show_matrix()

    # calculating energies and Non-Adiabatic Coupling
    E, D = Ene_NAC(E1,E2,P12,P21,par["dt_nuc"])

    print "E matrix is";  E.show_matrix()
    print "D matrix is";  D.show_matrix()


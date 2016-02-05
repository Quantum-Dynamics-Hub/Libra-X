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

#**********************************************************************************
#* This program extracts from the gamess output file the parameters
#* like forces , eigenenergies, eigenfunctions ,and atomic basis sets .
#* The forces are used for simulating Classical MD on Libra.
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

def detect_and_extract(l_gam,runtype,basis_sets):

    params = {}
    # detect the columns showing parameters
    detect(l_gam,params,runtype)    

    # extract the parameters from the columns detected
    extract(l_gam,params)

    # ****** construct atomic orbitals *****************************

    ao = ao_basis(params,basis_sets)

    return ao, params["E"], params["C"], params["gradient"]

def gamess_to_libra(inputs):

    # *************** open the GAMESS output file  ******************
    # 1st
    f_gam1 = open(inputs["gamess_out1"],"r")
    l_gam1 = f_gam1.readlines()
    f_gam1.close()

    # 2nd
    f_gam2 = open(inputs["gamess_out2"],"r")
    l_gam2 = f_gam2.readlines()
    f_gam2.close()

    # detect and extract parameters from GAMESS output file and create AO constructor
    ao1, E1, C1 ,gradient = detect_and_extract(l_gam1,inputs["runtype"],inputs["basis_sets"])
    ao2, E2, C2 ,gradient = detect_and_extract(l_gam2,inputs["runtype"],inputs["basis_sets"])

    # calculate overlap matrix of atomic orbitals and eigenfunctions
    P11, P22, P12 , P21 = overlap(ao1,ao2,C1,C2,inputs["basis_sets"])

    print "P11 and P22 matrixes should show orthogonality"
    print "P11 is"
    P11.show_matrix()
    print "P22 is"
    P22.show_matrix()

    print "P12 and P21 matrixes show overlap of different eigenfunctions "
    print "P12 is"
    P12.show_matrix()
    print "P21 is"
    P21.show_matrix()

    # calculating energies and Non-Adiabatic Coupling
    E, D = Ene_NAC(E1,E2,P12,P21,inputs["dt_nuc"])

    print "E matrix is"
    print E.show_matrix()
    print "D matrix is"
    print D.show_matrix()

    return

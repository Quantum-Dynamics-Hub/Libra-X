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

#**********************************************************
#* This program extracts parameters from the GAMESS output
#* and communicate them to Libra modules.
#* Inside the modules, gradients acting on atoms are used
#* for classical molecular dynamics
#* and eigenenergies and eigenfunctions are used
#* for simulating excited electron dynamics.
#**********************************************************


from detect import *
from detect1 import *
from extract import *
from ao_basis import *
from overlap import *
#from test_AO import *

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

def detect_and_extract(l_gam):

    params = {}

    # ****** detect columns containing parameters*******************

    detect1(l_gam,params)           # for optimize
    #detect(l_gam,params)             # for single point

    # ****** extract parameters from the columns detected***********

    extract(l_gam,params)

    # ****** construct atomic orbitals *****************************

    ao = ao_basis(params)

    return ao, params["E"], params["C"], params["gradient"]

def gamess_to_libra():

    # *************** open the GAMESS output file  ******************
    f_gam1 = open("../input/exam01.out","r")
    #f_gam1 = open("../input/exam03_STO_single.out","r")
    #f_gam1 = open("../input/exam03_AM1_single.out","r")
    l_gam1 = f_gam1.readlines()
    f_gam1.close()

    # extracte parameters and create AO constructor
    ao1, E1, C1 ,gradient = detect_and_extract(l_gam1)

    # calculate overlap matrix of AO and eigenfunctions
    P11 , P22 = overlap(ao1,ao1,C1,C1)

    print "P11 matrix should show orthogonality"
    P11.show_matrix()
    
    return

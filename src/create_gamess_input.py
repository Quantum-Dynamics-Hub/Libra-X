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

## \create_gamess_input.py
#  This module defines the functions to create a GAMESS input files as JOB.inp.


import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd= "/projects/academic/alexeyak/kosukesa/libracode-code/"
print "Current working directory", cwd
sys.path.insert(1,cwd+"_build/src/mmath")
sys.path.insert(1,cwd+"_build/src/chemobjects")
sys.path.insert(1,cwd+"_build/src/hamiltonian")
sys.path.insert(1,cwd+"_build/src/dyn")
#sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/hamiltonian_atomistic")

print "\nTest 1: Importing the library and its content"
from libmmath import *
from libchemobjects import *
from libhamiltonian import *
from libdyn import *
#from cyghamiltonian_atomistic import *


def keep_GMS_INP_format(GMS_DIR,GMS_JOB):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] GMS_DIR : The directory  where GAMESS is executed 
    # \param[in] GMS_JOB : The name of GAMESS input file
    # This function returns the format of GAMESS information before 
    # coordinates of atoms.
    #
    # Used in:  main.py/main/initial_gamess_exe    

    filename = GMS_DIR + GMS_JOB + ".inp"
    f = open(filename,"r")
    l_gam_for = f.readlines()
    f.close()

    N = len(l_gam_for)
    for i in range(0,N):
        s = l_gam_for[i].split()
        if len(s) > 0 and s[0] == "$DATA":
            ikeep = i + 2
            break

    l_gam_for[ikeep+1:N] = []

    return l_gam_for
    
def write_GMS_INP(l_atoms,l_charges,l_gam_for,mol,GMS_JOB,GMS_DIR):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_atoms   : a list of atoms
    # \param[in] l_charges : a list of atomic charges
    # \param[in] l_gam_for : includes format of GAMESS input file 
    # \param[in] GMS_DIR : The directory  where GAMESS is executed
    # \param[in] GMS_JOB : The name of GAMESS input file
    # This function returns the GAMESS input file
    #
    # Used in:  main.py/main/nve/nve_MD

    B_to_A = 0.529177208 # Bohr to Angstrom
    filename = GMS_DIR + GMS_JOB + ".inp"
    g = open(filename,"w")

    for i in range(0,len(l_gam_for)):
        g.write(l_gam_for[i])

    g.write("\n")
    Natoms = len(l_atoms)
    for i in range(0,Natoms):
        g.write("%s   %2.1f    %12.7f    %12.7f    %12.7f  \n" \
                % (l_atoms[i], l_charges[i], B_to_A*mol.q[3*i], B_to_A*mol.q[3*i+1], B_to_A*mol.q[3*i+2]) )

    g.write(" $END \n")

    

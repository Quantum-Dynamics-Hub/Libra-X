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

#*********************************************************************************
# This module calculates the overlap matrixes of atomic and molecular orbitals.
# When ab initio calculation is done, 
# the overlap matrix of atomic orbitals is necessary for calculating 
# the overlap matrix of eigenfunctions.
# However, When semi-empirical calculation is done, it is not.
#********************************************************************************

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
cwd = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code"
print "Using the Libra installation at", cwd
sys.path.insert(1,cwd+"/_build/src/mmath")
sys.path.insert(1,cwd+"/_build/src/qchem")

from libmmath import *
from libqchem import *

def AO_overlap(ao_i, ao_j):
    # calculate overlap matrix of atomic orbitals.

    Norb = len(ao_i)

    S = MATRIX(Norb,Norb)

    # overlap matrix of S
    for i in range(0,Norb): # all orbitals
        for j in range(0,Norb):
            S.set(i,j,gaussian_overlap(ao_i[i],ao_j[j]))

    return S

def MO_overlap(ao_i, ao_j, Ci, Cj, basis_option):
    # calculate overlap matrix of MOs:
    # <MO(i)|MO(j)> 

    Norb = len(ao_i)
    P = MATRIX(Norb, Norb)

    if basis_option == 1: # ab initio calculation
        S = AO_overlap(ao_i, ao_j)
        P = Ci.T() * S * Cj
    elif basis_option == 2: # semi empirical calculation
        P = Ci.T() * Cj
    else:
        print "basis_sets has an illegal value, so stop"
        sys.exit()
    
    return P


def overlap(ao1,ao2,C1,C2,basis_sets):
    # this is mostly a test function

    S11 = AO_overlap(ao1,ao1)
    S22 = AO_overlap(ao2,ao2)
    S12 = AO_overlap(ao1,ao2)
    S21 = AO_overlap(ao2,ao1)

    P11 = MO_overlap(ao1,ao1,C1,C1,basis_sets)
    P22 = MO_overlap(ao2,ao2,C2,C2,basis_sets)
    P12 = MO_overlap(ao1,ao2,C1,C2,basis_sets)
    P21 = MO_overlap(ao2,ao1,C2,C1,basis_sets)
    
    return P11, P22, P12, P21

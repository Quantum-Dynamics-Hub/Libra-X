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

def AO_overlap(ao1,ao2):
    # This function takes atomic orbitals from AO
    # and calculate overlap matrixes of AO and eigenfunctions.
    # Also, when atomic orbitals are at the same time
    # ,overlap matrix of eigenfunctions should show orthogonality.
    Norb = len(ao1)

    S = MATRIX(Norb,Norb)

    # overlap matrix of S
    for i in range(0,Norb): # all orbitals
        for j in range(0,Norb):
            S.set(i,j,gaussian_overlap(ao1[i],ao2[j]))

    print "The overlap matrix of S is"
    S.show_matrix()
    
    return S

def eigenfunction_overlap(S,Ci,Cj):

    # overlap matrix of < psi_i | psi_j > which should be Kronecker's delta
    #print "Orthogonality check!!!"
    #(C.T()*S*C).show_matrix()

    return (Ci.T()*S*Cj)

def overlap(ao1,ao2,C1,C2):

#    print len(aoi)
    S12 = AO_overlap(ao1,ao2)
    S21 = AO_overlap(ao2,ao1)

    P12 = eigenfunction_overlap(S12,C1,C2)
    P21 = eigenfunction_overlap(S21,C2,C1)
    

    return P12, P21

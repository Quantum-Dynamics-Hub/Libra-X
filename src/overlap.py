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

## \file overlap.py
# This program implements the module that calculates
# the overlap matrixes of atomic and molecular orbitals with different time steps.
# This returns the overlap matrix of molecular orbitals like  <MO(t)|MO(t+dt)>.


import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def AO_overlap(ao_i, ao_j):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] ao_i, ao_j : atomic orbital basis at different time step.
    # This function returns overlap matrix of atomic orbitals with different time step
    # like <AO(t)|AO(t+dt)>.
    #
    # Used in: overlap.py/overlap

    Norb = len(ao_i)

    S = MATRIX(Norb,Norb)

    # overlap matrix of S
    for i in range(0,Norb): # all orbitals
        for j in range(0,Norb):
            S.set(i,j,gaussian_overlap(ao_i[i],ao_j[j]))

    return S

def MO_overlap(S,ao_i, ao_j, Ci, Cj, basis_option):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] S : overlap matrix of atomic orbitals
    # \param[in] ao_i, ao_j : atomic orbital basis at different time step.
    # \param[in] Ci, Cj : molecular coefficients at different time step.
    # \param[in] basis_option : "= 1" means ab initio and "= 2" semi empirical .
    # This function returns overlap matrix of molecular orbitals with different time step
    # like <MO(t)|MO(t+dt)>.
    #
    # Used in: overlap.py/overlap

    Norb = len(ao_i)
    P = MATRIX(Norb, Norb)

    if basis_option == 1: # ab initio calculation
#        S = AO_overlap(ao_i, ao_j)
        P = Ci.T() * S * Cj
    elif basis_option == 2: # semi empirical calculation
        P = Ci.T() * Cj
    else:
        print "basis_sets has an illegal value, so stop"
        sys.exit()
    
    return P


def overlap(ao1,ao2,C1,C2,basis_sets):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] ao1, ao2 : atomic orbital basis at different time step.
    # \param[in] C1, C2 : molecular coefficients at different time step.
    # \param[in] basis_option : "= 1" means ab initio and "= 2" semi empirical .
    # This function returns overlap matrix of atomic orbitals with different time step
    # like <MO(t)|MO(t+dt)>.
    #
    # Used in: gamess_to_libra.py/gamess_to_libra
    # this is mostly a test function

    S11 = AO_overlap(ao1,ao1)
    S22 = AO_overlap(ao2,ao2)
    S12 = AO_overlap(ao1,ao2)
    S21 = AO_overlap(ao2,ao1)

    P11 = MO_overlap(S11,ao1,ao1,C1,C1,basis_sets)
    P22 = MO_overlap(S22,ao2,ao2,C2,C2,basis_sets)
    P12 = MO_overlap(S12,ao1,ao2,C1,C2,basis_sets)
    P21 = MO_overlap(S21,ao2,ao1,C2,C1,basis_sets)
    
    return P11, P22, P12, P21

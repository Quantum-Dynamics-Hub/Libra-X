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

## \file moment.py
# This program implements the module that calculates and returns
# the dipole moment matrixes at given space coordinates r like  <MO(t)| r |MO(t+dt)>.
#
# Used in: main.py/main/nve_MD/gamess_to_libra


import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_qchem_path"])

from libmmath import *
from libqchem import *

def AO_moment(ao_i, ao_j,g):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] ao_i, ao_j : atomic orbital basis at different time step.
    # \param[in] g : PrimitiveG object at given space
    # This function returns dipole moment matrix of atomic orbitals with different time step
    # like <AO(t)| r |AO(t+dt)>.
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/moment

    Norb = len(ao_i)
    D = MATRIX(Norb,Norb)

    
    # overlap matrix of D
    for i in range(0,Norb): # all orbitals
        for j in range(0,Norb):
            D.set(i,j,gaussian_moment(ao_i[i],g,ao_j[j]))

    return D

def MO_moment(D,ao_i,ao_j,Ci,Cj):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] D : dipole moment matrix of atomic orbitals
    # \param[in] ao_i, ao_j : atomic orbital basis at different time step.
    # \param[in] Ci, Cj : molecular coefficients at different time step.
    # This function returns overlap matrix of molecular orbitals with different time step
    # like <MO(t)| r |MO(t+dt)>.
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/overlap

    Norb = len(ao_i)
    DM = MATRIX(Norb, Norb)

    #S = AO_overlap(ao_i, ao_j)
    DM = Ci.T() * D * Cj

    return DM


def moment(ao1,ao2,C1,C2,g):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] ao_i, ao_j : atomic orbital basis at different time step.
    # \param[in] Ci, Cj : molecular coefficients at different time step.
    # \param[in] g : PrimitiveG vector at given space
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra

    # define PrimitiveG object to calculate dipole moment
    D11 = AO_moment(ao1,ao1,g)
    D22 = AO_moment(ao2,ao2,g)
    D12 = AO_moment(ao1,ao2,g)
    D21 = AO_moment(ao2,ao1,g)

    M11 = MO_moment(D11,ao1,ao1,C1,C1)
    M22 = MO_moment(D22,ao2,ao2,C2,C2)
    M12 = MO_moment(D12,ao1,ao2,C1,C2)
    M21 = MO_moment(D21,ao2,ao1,C2,C1)
    
    return M11, M22, M12, M21

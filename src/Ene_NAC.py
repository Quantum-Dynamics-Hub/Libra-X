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
## \file Ene_NAC.py
# This module defines the functions that return time-averaged energy and
# the Non-Adiabatic couplings (NACs).

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def NAC(P12,P21,dt_nucl):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] P12, P21 : overlap matrix of molecular orbitals at different time step.
    # \param[in] dt_nucl  : time step width of nuclear motion
    # This function returns Non-Adiabatic Couplings(NACs)
    #
    # Used in: gamess_to_libra.py/gamess_to_libra

    Norb = P12.num_of_rows
    D = MATRIX(Norb,Norb)

    D = 0.50/dt_nucl * ( P12 - P21 )

    return D

def average_E(E1,E2):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] E1, E2 : molecular energies at different time step.
    # This function returns the time-averaged molecular energies.
    #
    # Used in: gamess_to_libra.py/gamess_to_libra

    Norb = E1.num_of_rows
    E = MATRIX(Norb,Norb)

    E = 0.50 * (E1 + E2)

    return E

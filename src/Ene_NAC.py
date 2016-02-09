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

from libmmath import *
from libqchem import *


# This program calculates the averaged molecular energies and 
# NACs necessary for simulating excited electron dynamics.

def NAC(P12,P21,dt_nuc):
    Norb = P12.num_of_rows
    D = MATRIX(Norb,Norb)

    D = 0.50/dt_nuc * ( P12 - P21 )
    return D

def average_E(E1,E2):
    Norb = E1.num_of_rows
    E = MATRIX(Norb,Norb)

    E = 0.50 * (E1 + E2)
    return E

def Ene_NAC(E1,E2,P12,P21,dt_nuc):
    E = average_E(E1,E2)
    D = NAC(P12,P21,dt_nuc)

    return E, D

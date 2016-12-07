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

## \file include_mm.py
# This module defines a function initializing "Hamiltonian" objects for MM part.

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def init_hamiltonian_mm(syst, ff):
    #\param[in,out]     syst System object including atomic coordinates and connectivity information
    #\param[in]          ff  Forcefield object defining the MM interactions
    #
    # returned variable -> ham_mm  A list of MM Hamiltonians, linked to the corresponding system replicas

    ham_mm = []

    for i in xrange(len(syst)):

        ham = Hamiltonian_Atomistic(1, 3*syst[i].Number_of_atoms)
        ham.set_Hamiltonian_type("MM")
        atlst1 = range(1,syst[i].Number_of_atoms+1)
        ham.set_interactions_for_atoms(syst[i], atlst1, atlst1, uff, 1, 0) # 0 - verb, 0 - assign_rings
        ham.set_system(syst[i])
        ham.show_interactions_statistics()
        ham.compute()
        ham_mm.append(ham)

    return ham_mm

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

## \file create_states.py
# This module defines the function which creates a list of ground and excited states.
# It outputs the key parameter named "excitation"

# ************************ CAUTION ***************************************
# This function creates only singlet type excitation states.
# ************************************************************************

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["libra_hamiltonian_path"] + "/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters")
sys.path.insert(1,os.environ["libra_mmath_path"])

#from libhamiltonian import *
from libcontrol_parameters import *
from libmmath import *

def create_states(Nmin,HOMO,Nmax,spin,flip):
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in]  Nmin  lowest molecular orbital taken for TD-SE calculation
    # \param[in]  HOMO  Highest Occupied Molecular Orbital
    # \param[in]  Nmin  Highest molecular orbital taken for TD-SE calculation
    # \param[in]  spin  spin is considered : option 0 -> no, 1 -> yes
    # \param[in]  flip  spin flip is considered if spin = 1: option 0 -> no, 1 -> yes
    # Used in:  run.py

    LUMO = HOMO + 1
    
    excitations = []

    if spin == 0:
        sp_st = [1]
    elif spin == 1:
        sp_st = [1,-1]
    else:
        print "Error in create_states: spin must be 0 or 1"
        print "Value given = ", spin
        print "Exiting..."
        sys.exit(0)

    # note "excitation" object.
    # "excitation" has 4 index : (from orbital, from spin, to orbital, to spin)
    # \orbital index  0 -> HOMO, 1 -> LUMO, 2 -> LUMO+1, -1 -> HOMO-1, etc....
    # \spin index 1 -> alpha -1 -> beta
    # e.g. (0,1,0,1) means ground state (no excitation)
    #      (0,1,1,1)       alpha electron in HOMO is excited to LUMO without spin flip
    #      (0,1,1,-1)      alpha electron in HOMO is excited to LUMO with spin flip (in this case, spin-orbital coupling should be included)

    excitations.append(excitation(0,1,0,1)) # Add a Ground State

    count=1
    if flip == 0: # Excited States without spin flip
        for sp in sp_st:
            for om in range(Nmin,HOMO+1):
                for uom in range(LUMO,Nmax+1):

                    excitations.append(excitation(om-HOMO,sp,uom-HOMO,sp))

    elif flip == 1 and spin == 1: # Excited states with spin flip
        for sp1 in sp_st:
            for sp2 in sp_st:
                for om in range(Nmin,HOMO+1):
                    for uom in range(LUMO,Nmax+1):

                        excitations.append(excitation(om-HOMO,sp1,uom-HOMO,sp2))

    else:
        print "Error in create_states: flip must be 0 or 1 and spin 0 or 1"
        print "And,when spin = 0, flip should not be 1"
        print "Value flip and spin given = ", flip, spin
        print "Exiting..."
        sys.exit(0)

    print "excitations length is ",len(excitations)
    #print excitations

    return excitations

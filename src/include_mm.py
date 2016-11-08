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
# This module defines functions initializing "System" and "Hamiltonian" objects for MM part.

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def init_system(mb_functional, R_vdw_on, R_vdw_off, elem_file, uff_file, ent_file):
    #\param[in]       mb_functional vdw functional type
    #\param[in] R_vdw_on, R_vdw_off vdw interaction is included in the range of R_vdw_off and R_vdw_on  
    #\param[in]           elem_file containing elements information
    #\param[in]            uff_file containing UFF parameters
    #\param[in]            ent_file containing molecular information
    #
    # returned variable -> syst_mm System object for MM part

    syst_mm = System()

    # Create Universe and populate it
    U = Universe(); LoadPT.Load_PT(U, elem_file)
    # Create force field
    uff = ForceField({"mb_functional":mb_functional,"R_vdw_on":R_vdw_on,"R_vdw_off":R_vdw_off })
    LoadUFF.Load_UFF(uff, uff_file)

    LoadMolecule.Load_Molecule(U, syst_mm, ent_file,"pdb")

    return syst_mm

def init_hamiltonian_ele(syst_mm, mb_functional, R_vdw_on, R_vdw_off, elem_file, uff_file):
    #\param[in,out]         syst_mm system object for MM part.
    #\param[in]       mb_functional vdw functional type
    #\param[in] R_vdw_on, R_vdw_off vdw interaction is included in the range of R_vdw_off R_vdw_off and R_vdw_on
    #\param[in]           elem_file containing elements information
    #\param[in]            uff_file containing UFF parameters 
    #
    # returned variable -> ham_mm  MM hamiltonian linked to syst_mm
    #                       el_mm  list of electronic(1,0)
    ham_mm = []
    el_mm = []

    for i in xrange(len(syst_mm)):

        # Create Universe and populate it                                                                                                                    
        U = Universe(); LoadPT.Load_PT(U, elem_file)
        # Create force field                                                                                                                                 
        uff = ForceField({"mb_functional":mb_functional,"R_vdw_on":R_vdw_on,"R_vdw_off":R_vdw_off })
        LoadUFF.Load_UFF(uff, uff_file)

        ham = Hamiltonian_Atomistic(1, 3*syst_mm[i].Number_of_atoms)
        ham.set_Hamiltonian_type("MM")
        atlst1=range(1,syst_mm[i].Number_of_atoms+1)
        ham.set_interactions_for_atoms(syst_mm[i], atlst1, atlst1, uff, 1, 0) # 0 - verb, 0 - assign_rings
        ham.set_system(syst_mm[i])
        ham.show_interactions_statistics()
        ham.compute()
        ham_mm.append(ham)

        el_mm.append(Electronic(1,0))

    return ham_mm, el_mm

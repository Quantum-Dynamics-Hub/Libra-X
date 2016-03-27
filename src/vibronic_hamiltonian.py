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
## \file vibronic_hamiltonian.py
# This module implements the function which creates vibronic hamiltonian object.

import os
import sys
import math

sys.path.insert(1,os.environ["libra_qchem_path"])
sys.path.insert(1,os.environ["libra_hamiltonian_path"] + "/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters")
sys.path.insert(1,os.environ["libra_mmath_path"])

from libqchem import *
from libcontrol_parameters import *
from libmmath import *

def vibronic_hamiltonian(params,E_mol_red,D_mol_red,suffix):
    ##
    # This function transforms the 1-electron orbital energies matrix and the matrix of 
    # nonadiabatic couplings into the matrix of excitation energies and corresponding couplings
    # between the considered excited Slater determinants (SD). In doing this, it will return 
    # the vibronic and electronic (in SD basis) Hamiltonians.
    # 
    # \param[in] params  contains the dictionary of the input parameters
    # \param[in] E_mol_red   the matrix of the 1-electron MO energies, in reduced space
    # \param[in] D_mol_red   the matrix of the NACs computed with the 1-electon MOs, in reduced space
    # \param[in] suffix the suffix to add to the file names for the files created in this function
    #
    # This function returns the electronic and vibronic Hamiltonians, H_el, and Hvib
    #
    # Used in: main.py/main/run_MD

    HOMO = params["HOMO"]
    min_shift = params["min_shift"]
    max_shift = params["max_shift"]

    states = params["excitations"]
    nstates = len(states)  
    H_el = MATRIX(nstates,nstates)  # electronic Hamiltonian
    D_el = MATRIX(nstates,nstates)  # nonadiabatic couplings
    flag = params["print_sd_ham"]

    # Excitation energy : 
    # ex) GS = (0,1,0,1) -> E_GS = 0
    #     SE0 = (0,1,1,1) -> E_SE0 = E(LUMO) - E(HOMO)
    #     SE1 = (0,1,2,1) -> E_SE1 = E(LUMO+1) -E(HOMO) 
    #     etc ........

    # E.g. 
    # HOMO = 5, min_shift = -2
    # [0,1,2,3,4,5,  6,7,8,9]
    #     [0,1,2,3,  4,5]
    #

    #============ EX energies ================
    for i in xrange(nstates):

        h_indx = states[i].from_orbit[0] + min_shift  # index of the hole orbital w.r.t. the lowest included in the active space orbital
        e_indx = states[i].to_orbit[0] + min_shift  # --- same, only for the electron orbital

        EX_ene = E_mol_red.get(e_indx,e_indx) - E_mol_red.get(h_indx,h_indx)
       
        H_el.set(i,i,EX_ene)


    #============== Couplings =================

    for I in range(0,nstates):
        h_indx_I = states[I].from_orbit[0] + min_shift
        e_indx_I = states[I].to_orbit[0] + min_shift  


        for J in range(0,nstates):
            h_indx_J = states[J].from_orbit[0] + min_shift
            e_indx_J = states[J].to_orbit[0] + min_shift

        coupled = 0
        #====== Check whether the determinants I and J are coupled =======
        if I!=J: # exclude self-coupling
            if h_indx_I == h_indx_J:  # same hole position, electronic transitions
                coupled = 1
                i = e_indx_I
                j = e_indx_J
            elif e_indx_I == e_indx_J: # same electronic position, hole transitions
                couple = 1
                i = h_indx_J  ### TO DO: Check the ordering of the indices 
                j = h_indx_I  ### END TO DO

        if coupled:
            D_el.set(I,J, -D_mol_red.get(i,j))
        else:
            D_el.set(I,J, 0.0)

      
    if params["print_sd_ham"] == 1:
        H_el.show_matrix(params["sd_ham"] + "SD_re_Ham_" + suffix)
        D_el.show_matrix(params["sd_ham"] + "SD_im_Ham_" + suffix)


    Hvib = CMATRIX(H_el, D_el) # vibronic Hamiltonian

    # Returned values:
    # H_el - electronic Hamiltonian - the real diagonal matrix with adiabatic energies of the states
    #         w.r.t the ground state energy
    # Hvib - vibronic Hamiltonian - the complex-valued matrix, also containing nonadiabatic couplings
    # on the off-diagonals
    return H_el, Hvib

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

def vibronic_hamiltonian(params,E_mol,D):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params :  contains input parameters , in the directory form
    # \param[in] E_mol : time-averaged molecular energies
    # \param[in] D : Non-Adiabatic couplings
    # This function returns the vibronic hamiltonian "Hvib".
    #
    # Used in: main.py/main/run_MD

    HOMO = params["HOMO"]
    states = params["excitations"]
    nstates = len(states)
    Hvib = CMATRIX(nstates,nstates)
    nac = MATRIX(nstates,nstates)

    Nmin = states[1].from_orbit[0] + HOMO - 1
    Nmax = states[-1].to_orbit[0] + HOMO -1

    # Excitation energy : sum of the molecular orbital energies occupied by electrons.
    # E_I = sum_Ii e_Ii
    # ex) GS = [1,-1,2,-2,3,-3] -> E_GS = 2 * (e1 + e2 + e3)
    #     SE0 = [4,-1,2,-2,3,-3] -> E_SE0 = E_GS + e4 - e1
    #     SE1 = [1,-1,4,-2,3,-3] -> E_SE1 = E_GS + e4 - e2 
    #     etc ........

    # GS energy
    E0 = 0.0
    for i in range(Nmin,HOMO+1):
        E0 += 2 * E_mol.get(i-Nmin,i-Nmin)

    # EX energy
    for i in range(1,nstates):

        h_indx = states[i].from_orbit[0] + HOMO - 1
        e_indx = states[i].to_orbit[0] + HOMO - 1
        EX_ene = E0 + E_mol.get(e_indx-Nmin,e_indx-Nmin) - E_mol.get(h_indx-Nmin,h_indx-Nmin)

        Hvib.set(i,i,EX_ene,0.0)
    #print "Ex_ene="
    #Hvib.show_matrix()

    # NACs : a value exists when 2 excitonic wavefunctions differ only one electron position.
    # D_(I,J) = d_(i,j)
    # ex)  GS = [1,-1,2,-2,3,-3] 
    #     SE0 = [4,-1,2,-2,3,-3]    SE3 = [5,-1,2,-2,3,-3]
    #     SE1 = [1,-1,4,-2,3,-3]    SE4 = [1,-1,5,-2,3,-3]
    #     SE2 = [1,-1,2,-2,4,-3]    SE5 = [1,-1,2,-2,5,-3]
    #   D_(GS,SE0) = d_(1,4)
    #   D_(GS,SE1) = d_(2,4)
    #  D_(SE0,SE1) = d_(2,1) 
    #  D_(SE0,SE2) = d_(3,1) 
    #  D_(SE0,SE3) = d_(4,5)
    #  etc.......

    for i in range(1,nstates):
        h_indx_i = states[i].from_orbit[0] + HOMO - 1
        e_indx_i = states[i].to_orbit[0] + HOMO - 1

        # EX -> GS
        Hvib.set(i,0,0.0,-D.get(e_indx_i-Nmin,h_indx_i-Nmin))
        nac.set(i,0,-D.get(e_indx_i-Nmin,h_indx_i-Nmin))
        # GS -> EX
        Hvib.set(0,i,0.0,-D.get(h_indx_i-Nmin,e_indx_i-Nmin))
        nac.set(0,i,-D.get(h_indx_i-Nmin,e_indx_i-Nmin))

        for j in range(1,nstates):
            if not i == j:
                h_indx_j = states[j].from_orbit[0] + HOMO - 1
                e_indx_j = states[j].to_orbit[0] + HOMO - 1

                # EX -> EX
                if not e_indx_i == e_indx_j and h_indx_i == h_indx_j: # difference of the orbital occupied by excited electron 
                    Hvib.set(i,j,0.0,-D.get(e_indx_i-Nmin,e_indx_j-Nmin))
                    nac.set(i,j,-D.get(e_indx_i-Nmin,e_indx_j-Nmin))
                if e_indx_i == e_indx_j and not h_indx_i == h_indx_j: # difference of the orbital occupied by left hole
                    Hvib.set(i,j,0.0,-D.get(h_indx_j-Nmin,h_indx_i-Nmin))
                    nac.set(i,j,-D.get(h_indx_j-Nmin,h_indx_i-Nmin))


    print "D="
    D.show_matrix()
    print "nac(Im(Hvib))="
    nac.show_matrix()
    print "Hvib ="
    Hvib.show_matrix()

    return Hvib

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

def vibronic_hamiltonian(params,E_mol,D_mol,ite):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params  contains input parameters , in the directory form
    # \param[in] E_mol   time-averaged molecular energies on MO basis
    # \param[in] D_mol   Non-Adiabatic couplings on MO basis
    # \param[in] ite     The number of iteration 
    # This function returns the vibronic hamiltonian "Hvib".
    #
    # Used in: main.py/main/run_MD

    HOMO = params["HOMO"]
    states = params["excitations"]
    nstates = len(states)
    Hvib = CMATRIX(nstates,nstates)
    E_SD = MATRIX(nstates,nstates)
    D_SD = MATRIX(nstates,nstates)

    Nmin = states[-1].from_orbit[0] + HOMO - 1
    Nmax = states[-1].to_orbit[0] + HOMO -1

    # Excitation energy : 
    # ex) GS = (0,1,0,1) -> E_GS = 0
    #     SE0 = (0,1,1,1) -> E_SE0 = E(LUMO) - E(HOMO)
    #     SE1 = (0,1,2,1) -> E_SE1 = E(LUMO+1) -E(HOMO) 
    #     etc ........

    # EX energy
    for i in range(1,nstates):

        h_indx = states[i].from_orbit[0] + HOMO - 1
        e_indx = states[i].to_orbit[0] + HOMO - 1
        EX_ene = E_mol.get(e_indx-Nmin,e_indx-Nmin) - E_mol.get(h_indx-Nmin,h_indx-Nmin)

        Hvib.set(i,i,EX_ene,0.0)
        E_SD.set(i,i,EX_ene)

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

    if 1==1: # for debug mode
        print "(e_indx,h_indx) : e_indx is MO index of excited electron and h_indx is that of left hole"

    for i in range(1,nstates):
        h_indx_i = states[i].from_orbit[0] + HOMO - 1
        e_indx_i = states[i].to_orbit[0] + HOMO - 1

        # EX -> GS
        Hvib.set(i,0,0.0,-D_mol.get(e_indx_i-Nmin,h_indx_i-Nmin))
        D_SD.set(i,0,-D_mol.get(e_indx_i-Nmin,h_indx_i-Nmin))
        # GS -> EX
        Hvib.set(0,i,0.0,-D_mol.get(h_indx_i-Nmin,e_indx_i-Nmin))
        D_SD.set(0,i,-D_mol.get(h_indx_i-Nmin,e_indx_i-Nmin))

        # for debug mode
        if 1==1 : # for debug mode
            print "Imaginary part of Hvib(i=%i ((%i,%i) state) ,j=0 ((0,0) state) ) is -D(%i,%i)" %(i,e_indx_i,h_indx_i,e_indx_i,h_indx_i)
            print "Imaginary part of Hvib(i=0 ((0,0) state) ,j=%i ((%i,%i) state) ) is -D(%i,%i)" %(i,e_indx_i,h_indx_i,h_indx_i,e_indx_i)
        

        for j in range(1,nstates):
            if not i == j:
                h_indx_j = states[j].from_orbit[0] + HOMO - 1
                e_indx_j = states[j].to_orbit[0] + HOMO - 1

                # EX -> EX
                if not e_indx_i == e_indx_j and h_indx_i == h_indx_j: # difference of the orbital occupied by excited electron 
                    
                    Hvib.set(i,j,0.0,-D_mol.get(e_indx_i-Nmin,e_indx_j-Nmin))
                    D_SD.set(i,j,-D_mol.get(e_indx_i-Nmin,e_indx_j-Nmin))
                    if 1==1 : # for debug mode
                        print "Imaginary part of Hvib(i=%i ((%i,%i) state) ,j=%i ((%i,%i) state) ) is -D(%i,%i)" %(i,e_indx_i,h_indx_i,j,e_indx_j,h_indx_j,e_indx_i,e_indx_j)
                if e_indx_i == e_indx_j and not h_indx_i == h_indx_j: # difference of the orbital occupied by left hole
                    
                    Hvib.set(i,j,0.0,-D_mol.get(h_indx_j-Nmin,h_indx_i-Nmin))
                    D_SD.set(i,j,-D_mol.get(h_indx_j-Nmin,h_indx_i-Nmin))
                    if 1==1 : # for debug mode
                        print "Imaginary part of Hvib(i=%i ((%i,%i) state) ,j=%i ((%i,%i) state) ) is -D(%i,%i)" %(i,e_indx_i,h_indx_i,j,e_indx_j,h_indx_j,h_indx_j,h_indx_i)

    print "D_mol="
    D_mol.show_matrix()
    print "nac(Im(Hvib))="
    D_SD.show_matrix()
    print "Hvib ="
    Hvib.show_matrix()

    ene_filename = params["sd_ham"] + "re_Ham_" + str(ite)
    nac_filename = params["sd_ham"] + "im_Ham_" + str(ite)

    E_SD.show_matrix(ene_filename)
    D_SD.show_matrix(nac_filename)

    return Hvib

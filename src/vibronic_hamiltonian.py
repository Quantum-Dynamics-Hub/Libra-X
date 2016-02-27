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

sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_qchem_path"])
from libmmath import *
from libqchem import *

def vibronic_hamiltonian(params,E_mol,D):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params :  contains input parameters , in the directory form
    # \param[in] E_mol : time-averaged molecular energies
    # \param[in] D : Non-Adiabatic couplings
    # This function returns the vibronic hamiltonian "Hvib".
    #
    # Used in: main.py/main/run_MD

    nstates = params["excitations"]
    Hvib = CMATRIX(nstates,nstates)
    nac = MATRIX(nstates,nstates)

    #states = params["states"]
    #act = params["active_space"]
    #Nmin = act[0]
    #Nmax = act[-1]
    #Hvib = CMATRIX(len(states),len(states))
    #nac = MATRIX(len(states),len(states))
    #>>>>>>> devel

    # Excitation energy : sum of the molecular orbital energies occupied by electrons.
    # E_I = sum_Ii e_Ii
    # ex) GS = [1,-1,2,-2,3,-3] -> E_GS = 2 * (e1 + e2 + e3)
    #     SE0 = [4,-1,2,-2,3,-3] -> E_SE0 = E_GS + e4 - e1
    #     SE1 = [1,-1,4,-2,3,-3] -> E_SE1 = E_GS + e4 - e2 
    #     etc ........

    for i in range(0,len(states)):

        ene = 0.0
        for ii in states[i][1]:
            ene += E_mol.get(abs(ii)-Nmin,abs(ii)-Nmin)

        Hvib.set(i,i,ene,0.0)
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

    for i in range(0,len(states)):
        for j in range(0,len(states)):
            if not i == j:
                dif_i = []
                dif_j = []
                for k in range(0,len(states[i][1])):
                    if not states[i][1][k] == states[j][1][k]:
                        dif_i.append(states[i][1][k])
                        dif_j.append(states[j][1][k])
                if len(dif_i) == 1: # difference of the orbital occupied by excited electron 
                    Hvib.set(i,j,0.0,-D.get(abs(dif_i[0])-Nmin,abs(dif_j[0])-Nmin))
                    nac.set(i,j,-D.get(abs(dif_i[0])-Nmin,abs(dif_j[0])-Nmin))
                elif len(dif_i) == 2: # difference of the orbital occupied by left hole
                    if dif_i[0] == dif_j[1]:
                        Hvib.set(i,j,0.0,-D.get(abs(dif_i[1])-Nmin,abs(dif_j[0])-Nmin))
                        nac.set(i,j,-D.get(abs(dif_i[1])-Nmin,abs(dif_j[0])-Nmin))
                    elif dif_i[1] == dif_j[0]:
                        Hvib.set(i,j,0.0,-D.get(abs(dif_i[0])-Nmin,abs(dif_j[1])-Nmin))
                        nac.set(i,j,-D.get(abs(dif_i[0])-Nmin,abs(dif_j[1])-Nmin))

    #Hvib.set(3,5,0.0,0.01)
    #Hvib.set(5,3,0.0,0.01)
    print "D="
    D.show_matrix()
    print "nac(Im(Hvib))="
    nac.show_matrix()
    print "Hvib ="
    Hvib.show_matrix()

    return Hvib

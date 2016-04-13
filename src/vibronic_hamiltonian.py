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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#sys.path.insert(1,os.environ["libra_qchem_path"])
sys.path.insert(1,os.environ["libra_hamiltonian_path"] + "/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters")
#sys.path.insert(1,os.environ["libra_mmath_path"])
#sys.path.insert(1,os.environ["libra_util_path"])
#sys.path.insert(1,os.environ["libra_converters_path"])

#from libqchem import *
from libcontrol_parameters import *
#from libmmath import *
#from libconverters import *
#from libutil import *

def pyxaid_states(states, min_shift, max_shift):
    ##
    # This function converts our Libra type states into Pyxaid convention
    # E.g.
    # HOMO = 3, min_shift = -1, max_shift = 1
    # [0,1,2,3,4,5,6,7]  <- original set
    #     [0,1,2]        <- reduced set
    # nocc = 2
    # active_space_sz = 3
    # so: gs = [1,-1,2,-2]

    active_space_sz = max_shift - min_shift + 1
    nstates = len(states)

    nocc = -min_shift + 1 # the number of double occupied orbitals

    gs = [] # ground states
    for i in xrange(nocc): 
        indx = i + 1  # +1 so we can distinguish 0 and -0
        gs.append(indx)
        gs.append(-indx)

    # Now lets create all excitations from the reference state
    res = []
    print "All generated configurations (pyxaid indexing)\n"
    for i in xrange(nstates):        
        es = list(gs)
        print "Configuration ", i, es
     
        # Here +1 is needed to be consistent with Pyxiad indexing convention
        a = states[i].from_orbit[0] - min_shift + 1
        a_s = states[i].from_spin[0]  # +1 for alp, -1 for bet
        b = states[i].to_orbit[0] - min_shift + 1
        b_s = states[i].to_spin[0]  # +1 for alp, -1 for bet

        # orbital indices with sign = spin-orbit indices
        A = a * a_s  
        B = b * b_s

        # Replace A with B
        indx = es.index(A)
        es[indx] = B

        #print "Replace orbital ",A, " at position ",indx," with orbital ",B
   
        res.append(es)

    return res
      



def update_vibronic_hamiltonian(ham_el, ham_vib, params,E_mol_red,D_mol_red,suffix):
    ##
    # This function transforms the 1-electron orbital energies matrix and the matrix of 
    # nonadiabatic couplings into the matrix of excitation energies and corresponding couplings
    # between the considered excited Slater determinants (SD). In doing this, it will return 
    # the vibronic and electronic (in SD basis) Hamiltonians.
    # 
    # \param[out] ham_el Electronic (adiabatic) Hamiltonian (MATRIX)
    # \param[out] ham_vib Vibronic Hamiltonian (CMATRIX)
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
    #     SE1 = (-1,1,1,1) -> E_SE1 = E(LUMO) - E(HOMO-1) 
    #     etc ........

    # E.g. 
    # HOMO = 3, min_shift = -1, max_shift = 1
    # [0,1,2,3,4,5,6,7]  <- original set
    #     [0,1,2]        <- reduced set
    #

    pyx_st = pyxaid_states(states, min_shift, max_shift)

    #============ EX energies ================
    for i in xrange(nstates):

        h_indx = states[i].from_orbit[0] - min_shift  # index of the hole orbital w.r.t. the lowest included in the active space orbital
        e_indx = states[i].to_orbit[0] - min_shift  # --- same, only for the electron orbital

        EX_ene = E_mol_red.get(e_indx,e_indx) - E_mol_red.get(h_indx,h_indx)
       
        H_el.set(i,i,EX_ene)
        ham_el.set(i,i,EX_ene)
        ham_vib.set(i,i,EX_ene, 0.0)

        if flag==1:
            print "excitation ", i, " h_indx = ", h_indx, " e_indx = ", e_indx, " EX_ene = ", EX_ene
            print "source orbital energy = ",E_mol_red.get(h_indx,h_indx)
            print "target orbital energy = ",E_mol_red.get(e_indx,e_indx)


    #============== Couplings =================
    print "Do couplings"
    for I in xrange(nstates):
        for J in xrange(nstates):

            if I != J:
                st_I = Py2Cpp_int(pyx_st[I])
                st_J = Py2Cpp_int(pyx_st[J])
                print "st_I= ", pyx_st[I] # show_vector(st_I)
                print "st_J= ", pyx_st[J] # show_vector(st_J)
                coupled,a,b = delta(st_I, st_J)
                print "delta returns ", coupled, a, b
            
                if coupled:
                    i = abs(a) - 1
                    j = abs(b) - 1
                    print "pair of SD (",I,",",J,") is coupled via orbitals(in reduced space) ", i,j 
                    D_el.set(I,J, -D_mol_red.get(i,j))
                    ham_vib.set(I,J,0.0, -D_mol_red.get(i,j))
                else:
                    D_el.set(I,J, 0.0)
                    ham_vib.set(I,J, 0.0, 0.0)
                print "\n"

      
    if params["print_sd_ham"] == 1:
        H_el.show_matrix(params["sd_ham"] + "SD_re_Ham_" + suffix)
        D_el.show_matrix(params["sd_ham"] + "SD_im_Ham_" + suffix)


    #Hvib = CMATRIX(H_el, D_el) # vibronic Hamiltonian


    # Returned values - actually we just update matrices ham_el and ham_vib 
    # 
    # !!! Don't create new instances (new MATRIX and CMATRIX objects) 
    #
    # ham_el - electronic Hamiltonian - the real diagonal matrix with adiabatic energies of the states
    #         w.r.t the ground state energy
    # ham_vib - vibronic Hamiltonian - the complex-valued matrix, also containing nonadiabatic couplings
    # on the off-diagonals   

#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Kosuke Sato, Alexey V. Akimov
#* 
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file gamess_to_libra.py 
# This module implements the functions that extract parameters from the gamess output file:
# atomic forces , molecular energies, molecular orbitals, and atomic basis information.
# The forces are used for simulating Classical MD on Libra 
# and the others for calculating time-averaged energies and Non-Adiabatic Couplings(NACs).

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from extract import *
from overlap import *
from Ene_NAC import *
from moment import *




def qe_to_libra(params, E, C, suff):
    ## 
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params :  contains input parameters , in the directory form
    # \param[in,out] E  :  molecular energies at "t" old, will be updated
    # \param[in,out] C  :  molecular coefficients at "t" old, will be updated
    # \param[in] suff : The suffix to add to the name of the output files
    # this suffix is now considered to be of a string type - so you can actually encode both the
    # iteration number (MD timestep), the nuclear cofiguration (e.g. trajectory), and any other
    # related information
    #
    # This function outputs the files for excited electron dynamics
    # in "res" directory.
    # It returns the forces which act on the atoms.
    # Also, it returns new atomic orbitals, molecular energies, and
    # molecular coefficients used for calculating time-averaged
    # molecular energies and Non-Adiabatic Couplings(NACs).
    #
    # Used in: md.py/run_MD

    # 2-nd file - time "t+dt"  new
    label, Q, R, Grad, E2, C2, ao2, tot_ene = extract(params["gms_out"],params["debug_gms_unpack"])

    # calculate overlap matrix of atomic and molecular orbitals
    P11, P22, P12, P21 = overlap(ao,ao2,C,C2,params["basis_option"])

    # calculate transition dipole moment matrices in the MO basis:
    # mu_x = <i|x|j>, mu_y = <i|y|j>, mu_z = <i|z|j>
    # this is done for the "current" state only    
    mu_x, mu_y, mu_z = transition_dipole_moments(ao2,C2)
    mu = [mu_x, mu_y, mu_z]

    if params["debug_mu_output"]==1:
        print "mu_x:";    mu_x.show_matrix()
        print "mu_y:";    mu_y.show_matrix()
        print "mu_z:";    mu_z.show_matrix()
 
    if params["debug_densmat_output"]==1:
        print "P11 and P22 matrixes should show orthogonality"
        print "P11 is";    P11.show_matrix()
        print "P22 is";    P22.show_matrix()
        print "P12 and P21 matrixes show overlap of MOs for different molecular geometries "
        print "P12 is";    P12.show_matrix()
        print "P21 is";    P21.show_matrix()


    ### TO DO: In the following section, we need to avoid computing NAC matrices in the full
    # basis. We will need the information on cropping, in order to avoid computations that 
    # we do not need (the results are discarded anyways)
    # calculate molecular energies and Non-Adiabatic Couplings(NACs) on MO basis
    E_mol = average_E(E,E2)
    D_mol = NAC(P12,P21,params["dt_nucl"])

    # reduce the matrix size
    E_mol_red = reduce_matrix(E_mol,params["min_shift"], params["max_shift"],params["HOMO"])
    D_mol_red = reduce_matrix(D_mol,params["min_shift"], params["max_shift"],params["HOMO"])
    ### END TO DO

    if params["print_mo_ham"]==1:
        E_mol.show_matrix(params["mo_ham"] + "full_re_Ham_" + suff)
        D_mol.show_matrix(params["mo_ham"] + "full_im_Ham_" + suff)
        E_mol_red.show_matrix(params["mo_ham"] + "reduced_re_Ham_" + suff)
        D_mol_red.show_matrix(params["mo_ham"] + "reduced_im_Ham_" + suff)

    # store "t+dt"(new) parameters on "t"(old) ones
    for i in range(0,len(ao2)):
        ao[i] = AO(ao2[i])
    E = MATRIX(E2)
    C = MATRIX(C2)

    # Returned data:
    # Grad: Grad[k][i] - i-th projection of the gradient w.r.t. to k-th nucleus (i = 0, 1, 2)
    # data: a dictionary containing transition dipole moments
    # E_mol: the matrix of the 1-el orbital energies in the full space of the orbitals
    # D_mol: the matrix of the NACs computed with 1-el orbitals. Same dimension as E_mol
    # E_mol_red: the matrix of the 1-el orbital energies in the reduced (active) space
    # D_mol_red: the matrix of the NACs computed with 1-el orbital. Same dimension as E_mol_red

    return tot_ene, Grad, mu, E_mol, D_mol, E_mol_red, D_mol_red

